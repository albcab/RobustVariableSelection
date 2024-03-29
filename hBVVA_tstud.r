
counts = function(full, uniques)
{
    n = rep(0, length(uniques))
    for (i in seq(along.with = uniques))
    {
        n[i] = sum(full == uniques[i])
    }
    return(n)
}

hBVVA2 <- function(iter = 1000, X, y, v=2, b1 = 1, b2 = 20, d1, d2, verbose=FALSE, scaling = TRUE){
  
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
  tau = 1
  lambda = rep(1, K)
  big_gamma = tau*lambda
  
  if (scaling)
    X = scale(X, scale = FALSE)
  else
    X <- scale(X)
  
  XX = t(X) %*% X
  
  XY = t(X) %*% Y_Star
  
  sigma = rep(var(Y_Star), n)
  G = rgamma(n, v/2, v*sigma/2)
  
  # set space to store the results
  beta.store=matrix(rep(0,iter*K),nrow=iter)
  sigma.store = matrix(rep(0,iter*n),nrow=iter)
  G.store = matrix(rep(0,iter*n),nrow=iter)
  lambda.store=matrix(ncol = K, nrow = iter)
  tau.store=rep(0,iter)
  w.store = rep(0,iter)
  v.store = rep(0,iter)
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
    E        = diag(G)
    TXE      = t(X) %*% E #efficiency 
    try(variance <- solve(TXE %*% X + diag((1/big_gamma)) )) #efficiency 
    mu       = variance %*%  (TXE %*% Y_Star)
    beta = MASS::mvrnorm(1, mu, variance)
    } else {
    u = rnorm(K, 0, sqrt(big_gamma))
    delta = rnorm(n)
    E_vec = sqrt(G)
    Phi = t(t(X) * outer(rep.int(1L, ncol(X)), E_vec))
    GPhi = t(Phi * outer(rep.int(1L, nrow(Phi)), big_gamma))
    aalpha = E_vec * Y_Star
    vv = Phi %*% u + delta
    ww = solve(Phi %*% GPhi + diag(n), aalpha - vv)
    beta = u + GPhi %*% ww
    }
    
    #resample tau
    xi = 1/rgamma(1, 1, rate = 1 + 1/tau)
    tau = 1/rgamma(1, (K+1)/2, rate = 1/xi + 1/2*sum(beta^2/lambda))
    
    # Resample lambda
    # for (k in  1:K) 
    loop_fun = function(lambda, beta) {
        vj = 1/rgamma(1, 1, rate = 1 + 1/lambda)
        lambda = 1/rgamma(1, 1, rate = 1/vj + beta^2/(2*tau))
        return(lambda)
    }
    lambda = mapply(loop_fun, lambda, beta)
    
    #step 1
    prob.new = alpha * (v/2)^(v/2)/gamma(v/2) * G^(v/2-1) * b2^b1/gamma(b1) * gamma(b1+v/2)/( (b2+v*G/2)^(b1+v/2) )
    # log_prob.new = (v/2) * log(v/2) - lgamma(v/2) + (v/2-1) * log(G) + b1 * log(b2) - lgamma(b1) + lgamma(b1 + v/2) - (b1 + v/2) * log(b2 + v*G/2)
    # prob.new = alpha * exp(log_prob.new)
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
      dgama = dgamma(G[j], v/2, v*unique_sigma/2)
      propto = c(table_counts*dgama, prob.new[j])
      new_j = sample(1:(M+1), 1, replace = T, prob = propto)
      if (new_j == M+1)
      {
        table_counts[M+1] = 1
        unique_sigma[M+1] = rgamma(1, b1 + (table_counts[M+1]*v/2), rate =  b2 + (v*G[j]/2))
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
      unique_sigma[m] = rgamma(1, b1 + (table_counts[m]*v/2), rate =  b2 + (v*sum(G[which_sigma])/2))
      sigma[which_sigma] = unique_sigma[m]
    }
    
    #step 0
    square.diff = (Y_Star - (X %*% beta))^2
    G = rgamma(n, (v+1)/2, (square.diff + v*sigma)/2)

    # #Find v
    # delta = sum((Y_Star - X %*% beta) ^ 2 / sigma)
    # v_fn = function(v) {
    #   v = exp(v)
    #   w = (v + n) / (v + delta)
    #   return(1 - digamma(v/2) + log(v/2) + log(w) - w + digamma((v+n)/2) - log((v+n) / 2))
    # }
    # v = exp(uniroot(v_fn, lower=log(1e-6), upper=log(20), extendInt="yes")$root)
        
    # Set new value of Big_ This hypervariance shrinks betas
    big_gamma = tau * lambda
    
    #step 7
    # As alpha_x increases the scale decreases of the gamma which increases alpha value
    alpha_x = rbeta(1, alpha+1, n)
    aw1 = d1 + K + 1
    aw2 = n * (d2 - log(alpha_x))
    a_extra = rbinom(1, 1, aw2 / (aw2 + aw1))  
    alpha = rgamma(1, d1 + M + a_extra, rate = (d2 - log(alpha_x)) )
    
    # store results of this iteration
    alpha.store[i] = alpha
    beta.store[i,] = beta
    sigma.store[i,] = sigma
    G.store[i,] = G
    K.store[i] = M
    tau.store[i]  = tau
    lambda.store[i,] = lambda
    v.store[i] = v
  }
  
  return(list('alpha'=alpha.store,'beta'=beta.store,'lambda'=lambda.store,'tau'=tau.store,'sigma'=sigma.store,'K'=K.store, 'G' = G.store, 'v'=v.store))
}
