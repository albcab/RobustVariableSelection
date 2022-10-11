# JUST SVSS

SVSS <- function(iter = 1000, X, y, v_0 = 0.0001, a1 = 5, a2 = 50, b1 = 1, b2 = 20, verbose=FALSE, scaling = TRUE){
  if(missing(X)) {
    stop('The X Parameter has not been supplied')
  }
  if(missing(y)) {
    stop('The y Parameter has not been supplied')
  }
  
  n = length(y)
  K = ncol(X)
  
  # set some initial values 
  big_gamma = rep(1, K)
  tau = c(rep(1, K))
  I=rep(0,K)
  w = runif(1)
  
  if (scaling)
    X = scale(X, scale = FALSE)
  # else
  #   X <- scale(X)
  
  # if (scaling)
    Y_Star = scale(y, scale = FALSE) #larger variances in the Y 
  # else
    # Y_Star = sapply( y, as.numeric )
  
  # Moved this to below
  # Should it be sd or sd^2?
  sigma = sd(Y_Star)^2
  
  # set space to store the results
  beta.store=matrix(rep(0,iter*K),nrow=iter)
  I.store=matrix(rep(0,iter*K),nrow=iter)
  tau.store=matrix(rep(0,iter*K),nrow=iter)
  w.store=rep(0,iter)
  w2.store = matrix(nrow = iter, ncol = K)
  sigma.store=rep(0,iter)
  c=10
  for (i in 1:iter)
  {
    if(verbose==TRUE && (i>(c/100)*iter) ){
      print(paste(c,'%'))
      c = c+10
    }
    
    # Resample Beta 
    if (n >= K) {  
    try(variance <-  solve(t(X) %*% X  + sigma * diag(1/big_gamma))) #more efficient, took out n 
    mu = variance %*% t(X) %*% Y_Star
    beta = MASS::mvrnorm(1, mu, sigma*variance) #added sigma multiplication
    } else {
    u = rnorm(K, 0, sqrt(big_gamma))
    delta = rnorm(n)
    v = (X/sqrt(sigma)) %*% u + delta
    GX = t(X/sqrt(sigma) * outer(rep.int(1L, nrow(X)), big_gamma))
    ww = solve((X/sqrt(sigma)) %*% GX + diag(n), Y_Star/sqrt(sigma) - v)
    beta = u + GX %*% ww
    }
    # Resample I and tau
    # for (k in  1:K) 
    loop_fun = function(b, t) {
      w1 = (1 - w) * v_0^(-0.5) * exp(-1 * (b^2 / (2 * v_0 * t)) ) 
      w2 = w * exp(-1 * (b^2 / (2 * t)) )
      if(w1 + w2 == 0)
      {
        print('w')
        w1 = 1
        w2 = 1
      }
      
      I = sample(c(v_0, 1), 1, prob=c(w1,w2))
      tau = 1/ (rgamma(n = 1, shape = a1 + 0.5 , rate = a2 + (b^2/(2 * I) )))

      w2s = w2/(w1+w2)
      return(c(I, tau, w2s))
    }
    Itw2 = mapply(loop_fun, beta, tau)
    I = Itw2[1,]
    tau = Itw2[2,]
    w2.store[i,] = Itw2[3,]
    
    # Resample w
    w = rbeta(1 , 1 + sum(I == 1), 1 + sum(I == v_0)  )
    
    # Resample sigma
    sigma = 1/ (rgamma(1, shape = b1 + n/2 , rate = b2 + (1/2 ) * sum((Y_Star - (X %*% beta ))^2) ))
    #sigma = 1/ (rgamma(1, shape = b1 , rate = b2 + (1/2 ) * sum((Y_Star - (X %*% beta))^2) ))
    
    # Set new value of Big_
    big_gamma = I * tau
    
    # store results of this iteration
    beta.store[i,]= beta
    I.store[i,]   = I
    tau.store[i,] = tau
    w.store[i]    = w
    sigma.store[i]= sigma
  }
  
  return(list('beta'=beta.store, 'I'=I.store, 'tau'=tau.store, 'w'=w.store, 'sigma'=sigma.store, 'w2'=w2.store))
}
