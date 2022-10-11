
counts = function(full, uniques)
{
    n = rep(0, length(uniques))
    for (i in seq(along.with = uniques))
    {
        n[i] = sum(full == uniques[i])
    }
    return(n)
}

BVVA2 <- function(iter = 1000, X, y, v_0 = 0.0001, a1 = 5, a2 = 50, b1 = 1, b2 = 20, d1, d2, verbose=FALSE, scaling = TRUE){
  
  if(missing(X)) {
    stop('The X Parameter has not been supplied')
  }
  if(missing(y)) {
    stop('The y Parameter has not been supplied')
  }
  
  # if (scaling)
    Y_Star = scale(y, scale = FALSE) #larger variances in the Y 
  # else
    # Y_Star = sapply( y, as.numeric )
  # Y.mu = mean(Y)
  
  n=length(y)
  K=ncol(X)
   
  # set some initial values 
  big_gamma = rep(1, K)
  tau = c(rep(1, K))
  I = rep(1,K)
  
  w = runif(1)
  
  if (scaling)
    X = scale(X, scale = FALSE)
  # else
  #   X <- scale(X)
  
  XX = t(X) %*% X
  
  XY = t(X) %*% Y_Star
  
  sigma = rep(var(Y_Star), n)
  
  # set space to store the results
  beta.store=matrix(rep(0,iter*K),nrow=iter)
  sigma.store = matrix(rep(0,iter*n),nrow=iter)
  c.store = sigma.store
  I.store=matrix(rep(0,iter*K),nrow=iter)
  tau.store=matrix(rep(1,iter*K),nrow=iter)
  w.store = rep(0,iter)
  w2.store = matrix(nrow = iter, ncol = K)
  K.store = w.store;  alpha.store = w.store;  max.variance = w.store; max.mu = w.store
  
  #Conentration parameter for DP
  alpha = rgamma(1,1,1)
  
  #Chinese restaurant process
  M = ifelse(round(alpha*log(n)), round(alpha*log(n)), 1)
  tables = 1:M
  customer = sample(1:M, n, replace = T)
  unique_sigma = rep(sigma[1], M)
  table_counts = counts(customer, tables)
  
  c = 10 
  # Run the Gibbs sampler
  for (i in 1:iter)
  {
    if (verbose) 
      cat('\r', M, '-', alpha, '-', quantile(unique_sigma, probs = c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95)))
    if(verbose==TRUE && (i>(c/100)*iter) ){
      print(paste(c,'%'))
      c = c+10
    }
    
    #step 3
    if (n >= K) {
    E        = diag(1/sigma)
    TXE      = t(X) %*% E #efficiency 
    try(variance <- solve(TXE %*% X + diag(1/big_gamma) )) #efficiency
    mu       = variance %*%  (TXE %*% Y_Star)
    beta = MASS::mvrnorm(1, mu, variance)
    } else {
    u = rnorm(K, 0, sqrt(big_gamma))
    delta = rnorm(n)
    E_vec = sqrt(1/sigma)
    Phi = t(t(X) * outer(rep.int(1L, ncol(X)), E_vec))
    GPhi = t(Phi * outer(rep.int(1L, nrow(Phi)), big_gamma))
    aalpha = E_vec * Y_Star
    v = Phi %*% u + delta
    ww = solve(Phi %*% GPhi + diag(n), aalpha - v)
    beta = u + GPhi %*% ww
    }
    
    # for (k in  1:K) 
    loop_fun = function(beta, I){
    
      #step 4
      tau = 1/ (rgamma(n = 1, shape = a1 + 0.5 , rate = a2 + (beta^2/(2 * I) )))    

      #step 5
      w1 = (1 - w) * v_0^(-0.5) * exp(-1 * (beta^2 / (2 * v_0 * tau)) ) 
      w2 = w * exp(-1 * (beta^2 / (2 * tau)) )
      I = sample(c(v_0, 1), 1, prob=c(w1,w2))     

      w2s = w2/(w1+w2)
      return(c(I, tau, w2s))
    }
    Itw2 = mapply(loop_fun, beta, I)
    I = Itw2[1,]
    tau = Itw2[2,]
    w2.store[i,] = Itw2[3,]
    
    #step 6
    w = rbeta(1 , 1 + sum(I == 1), 1 + sum(I == v_0)  )
    
    #step 1
    square.diff = (Y_Star - (X %*% beta))^2
    prob.new = ( alpha * (1/ sqrt(2 * pi)) * ((b2^b1)/gamma(b1))) * (((square.diff / 2) + b2 )^ (-1 *( b1 + 0.5))  ) * gamma(b1 + 0.5)
    # prob.new = alpha*(1/sqrt(2*pi))*gamma(1/2)*((square.diff/2)^(-1/2))
    for (j in 1:n)
    {
      which_table = which(customer[j] == tables)
      table_counts[which_table] = table_counts[which_table] - 1
      if (table_counts[which_table] == 0)
      {
        table_counts = table_counts[-which_table]
        unique_sigma = unique_sigma[-which_table]
        tables = tables[-which_table]
        M = M - 1
      }
      dnormal = dnorm(Y_Star[j], X[j,] %*% beta, unique_sigma)
      propto = c(table_counts*dnormal, prob.new[j])
      new_j = sample(1:(M+1), 1, replace = T, prob = propto)
      if (new_j == M+1)
      {
        table_counts[M+1] = 1
        unique_sigma[M+1] = 1/ (rgamma(1, b1 + table_counts[M+1] /2, rate =  b2 + ((square.diff[j])/2)))
        tables[M+1] = max(tables) + 1
        M = M + 1
      }
      else
      {
        table_counts[new_j] = table_counts[new_j] + 1
      }
      customer[j] = tables[new_j]
    }
    
    #step 2
    for (m in 1:M) 
    {
      which_sigma = which(tables[m] == customer)
      unique_sigma[m] = 1/ (rgamma(1, b1 + table_counts[m] /2, rate =  b2 + ((sum(square.diff[which_sigma]))/2)))
      # unique_sigma[m] = 1/ (rgamma(1, table_counts[m] /2, rate = ((sum(square.diff[which_sigma]))/2)))
      sigma[which_sigma] = unique_sigma[m]
    }
        
    # Set new value of Big_ This hypervariance shrinks betas
    big_gamma = ifelse(is.infinite(I * tau), big_gamma, I*tau)
    
    #step 7
    # As alpha_x increases the scale decreases of the gamma which increases alpha value
    alpha_x = rbeta(1, alpha+1, n)
    aw1 = d1 + K + 1
    aw2 = n * (d2 - log(alpha_x))
    a_extra = rbinom(1, 1, aw2 / (aw2 + aw1))  
    alpha = rgamma(1, d1 + M + a_extra, rate = (d2 - log(alpha_x)) )
    
    # store results of this iteration
    c.store[i,] = sapply(customer, function(x) which(x == tables))
    alpha.store[i] = alpha
    beta.store[i,] = beta
    sigma.store[i,] = sigma
    K.store[i] = M
    I.store[i,]    = I
    tau.store[i,]  = tau
    w.store[i]     = w
  }
  return(list('c'=c.store, 'w'=w.store,'alpha'=alpha.store,'beta'=beta.store,'I'=I.store,'tau'=tau.store,'sigma'=sigma.store,'K'=K.store, 'w2'=w2.store))

}
