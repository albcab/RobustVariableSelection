# JUST SVSS

hSVSS <- function(iter = 1000, X, y, v=2, b1 = 1, b2 = 20, verbose=FALSE, scaling = TRUE){
  if(missing(X)) {
    stop('The X Parameter has not been supplied')
  }
  if(missing(y)) {
    stop('The y Parameter has not been supplied')
  }
  
  n = length(y)
  K = ncol(X)
  
  # set some initial values 
  tau = 1
  lambda = rep(1, K)
  big_gamma = tau*lambda
  
  if (scaling)
    X = scale(X, scale = FALSE)
  else
    X <- scale(X)
  
  # if (scaling)
    Y_Star = scale(y, scale = FALSE) #larger variances in the Y 
  # else
    # Y_Star = sapply( y, as.numeric )
  
  # Moved this to below
  # Should it be sd or sd^2?
  sigma = sd(Y_Star)^2
  G = rgamma(1, v/2, v*sigma/2)
  
  # set space to store the results
  beta.store=matrix(rep(0,iter*K),nrow=iter)
  tau.store=rep(0,iter)
  sigma.store=rep(0,iter)
  G.store=rep(0,iter)
  v.store = G.store
  lambda.store=matrix(ncol = K, nrow = iter)
  c=10
  for (i in 1:iter)
  {
    if(verbose==TRUE && (i>(c/100)*iter) ){
      print(paste(c,'%'))
      c = c+10
    }

    #Resample Beta
    if (n >= K) {   
    try(variance <-  solve(t(X) %*% X  + diag(1/big_gamma))) #more efficient, took out n 
    mu = variance %*% t(X) %*% Y_Star
    beta = MASS::mvrnorm(1, mu, 1/G*variance) #added sigma multiplication
    } else {
    u = rnorm(K, 0, sqrt(1/G * big_gamma))
    delta = rnorm(n)
    vv = (X*sqrt(G)) %*% u + delta
    GX = t(X * outer(rep.int(1L, nrow(X)), 1/sqrt(G) * big_gamma))
    ww = solve((X*sqrt(G)) %*% GX + diag(n), Y_Star*sqrt(G) - vv)
    beta = u + GX %*% ww
    }
    # Resample lambda
    # for (k in  1:K) 
    loop_fun = function(lambda, beta) {
        vj = 1/rgamma(1, 1, 1 + 1/lambda)
        lambda = 1/rgamma(1, 1, 1/vj + beta^2/(2*tau*1/G))
        return(lambda)
    }
    lambda = mapply(loop_fun, lambda, beta)
    
    #resample tau
    xi = 1/rgamma(1, 1, 1 + 1/tau)
    tau = 1/rgamma(1, (K+1)/2, 1/xi + 1/(2*1/G)*sum(beta^2/lambda))
    
    # Resample sigma
    sigma = rgamma(1, b1 + (v/2), rate =  b2 + (v*G/2))
    G = rgamma(1, (v+n+K)/2, (sum((Y_Star - (X %*% beta ))^2) + sum(beta^2/big_gamma) + v*sigma)/2)

    # #Find v
    # delta = sum((Y_Star - X %*% beta) ^ 2) / sigma
    # v_fn = function(v) {
    #   v = exp(v)
    #   w = (v + n) / (v + delta)
    #   return(1 - digamma(v/2) + log(v/2) + log(w) - w + digamma((v+n)/2) - log((v+n) / 2))
    # }
    # v = exp(uniroot(v_fn, lower=log(1e-6), upper=log(20), extendInt="yes")$root)
    
    # Set new value of Big_
    big_gamma = tau * lambda
    
    # store results of this iteration
    beta.store[i,]= beta
    tau.store[i] = tau
    lambda.store[i,] = lambda
    sigma.store[i]= sigma
    G.store[i] = G
    v.store[i] = v
  }
  
  return(list('beta'=beta.store, 'lambda'=lambda.store, 'tau'=tau.store, 'sigma'=sigma.store, 'G'=G.store, 'v'=v.store))
}
