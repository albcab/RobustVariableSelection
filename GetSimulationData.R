require(MASS)
Simulate_Data <- function(n, K, rho, variance, seed=FALSE, noise=TRUE )
{
  if(missing(n)) {
    n   = 200
  }
  if(missing(K)) {
    K   = 200
  }
  if(missing(rho)) {
    rho = 0.9
  }
  if(missing(variance)) {
    variance = 1
  }
  if(seed != FALSE){
    set.seed(seed)
  }
  
  CovarX=matrix(rep(0,K*K),nrow=K)
  
  for (i in 1:K)
  {
    for (j in 1:K)
    {
      CovarX[i,j]=rho^(abs(i-j))
    }
  }
  
  X=mvrnorm(n,rep(0,K),CovarX)
  
  for (i in 1:K)
  {
    sdX=sd(X[,i])
    mea=mean(X[,i])
    for (j in 1:n)
    {
      X[j,i]= (X[j,i]-mea)/sdX
    }
  }
  
  beta1=rep(0,7)
  beta1[4]=4^1.25
  for (i in 1:3)
  {
    beta1[i]=i^1.25
    beta1[8-i]=i^1.25
  }
  
  beta0=rep(0,18)
  while(length(beta0) <K)
  {
    beta0=c(beta0,beta1,rep(0,18))
  }
  beta0=beta0[1:K]
  
  if(noise == TRUE){
    # Add a small amount of noise to betas if desired
    # Why do this??
    beta0 = beta0 + rnorm(K,0,0.5)
  }
  
  length(beta0)
  
  y = X%*%beta0
  
  if (variance == 1)
  {
    y = y + rnorm(n, 0, 1)
  }
  else if(variance == 3)
  {
    y[1:round(n/3)]                        = y[1:round(n/3)]                   + rnorm(round(n/3), 0, 0.5)
    y[round((n/3)+1): round((2*n)/3)]      = y[round((n/3)+1): round((2*n)/3)] + rnorm( length(y[round((n/3)+1): round((2*n)/3)]) , 0, 1)
    y[round( ((2*n)+1) /3):n]              = y[round( ((2*n)+1) /3):n]         + rnorm( length(y[round( ((2*n)+1) /3):n]   ) , 0, 1.5)
  }
  else if(variance == 3.1)
  {
    y[1:round(10*n/50)]                        = y[1:round(10*n/50)]                   + rnorm(round(10*n/50), 0, sqrt(0.5))
    y[round((10*n/50)+1):round((20*n)/50)]     = y[round((10*n/50)+1):round((20*n)/50)]+ rnorm(round(10*n/50), 0, 1)
    y[round((20*n/50)+1):round((30*n)/50)]     = y[round((20*n/50)+1):round((30*n)/50)]+ rnorm(round(10*n/50), 0, sqrt(1.5))
    y[round((30*n/50)+1):round((40*n)/50)]     = y[round((30*n/50)+1):round((40*n)/50)]+ rnorm(round(10*n/50), 0, sqrt(2))
    y[round((40*n/50)+1):n]                    = y[round((40*n/50)+1):n]               + rnorm(round(10*n/50), 0, sqrt(2.5))
    if (round(n/50) >= 1) #add 1,2,4 outliers for 50,100,200
        for (i in 1:round(n/50))
            y[round((6-i)*10*n/50)] = sqrt(10/((6-i)*0.5))*y[round((6-i)*10*n/50)]
  }
  
  varianceY=sd(y)
  
  X = as.data.frame(X)
  
  return(list(X=X,y=y,varianceY=varianceY,beta0=beta0))
}
