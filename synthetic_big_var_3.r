# Script to run over range of n and p
# Make sure you create the folder 'Results' in your working directory

source("BVVA.R")
source("hBVVA.R")
source("SVSS.R")
source("hSVSS.R")

source('GetSimulationData.R')

nlist <- c(10,20,50,100,200); 
plist <- c(20,50,100,200);

rho <- 0
variance <- 3.1
Niter <- 10000
Nsam <- 5000
K = 10

for (n in nlist){
  for (p in plist){
    sim = vector(mode = 'list', length = K)
    mSVSS = sim; h_SVSS = sim; mBVVA = sim; hBVVA = sim
    set.seed(999)
    for (k in 1:K)
    {
    data      = Simulate_Data(n,p,rho,variance, seed = F, noise=F)
    sim[[k]] = list()
    sim[[k]]$X         = data[1][[1]]
    sim[[k]]$y         = data[2][[1]]
    sim[[k]]$var   = data[3][[1]]
    sim[[k]]$betas  = data[4][[1]]
    }
    
    for (k in 1:K)
    {
    print(paste("Currently at n", n, "p", p, "k", k))
    # Ishwaran Method
	set.seed(999)
    start_time <- Sys.time()
    mSVSS[[k]] = SVSS(X=sim[[k]]$X, y=sim[[k]]$y, iter = Niter, a1=2.01, a2=1, b1=2.01, b2=1, verbose=T)
    end_time <- Sys.time()
    mSVSS[[k]]$time = end_time - start_time
	
    # horseshoe Method
	set.seed(999)
    start_time <- Sys.time()
    h_SVSS[[k]] = hSVSS(X=sim[[k]]$X, y=sim[[k]]$y, iter = Niter, b1=2.01, b2=1, verbose=T)
    end_time <- Sys.time()
    h_SVSS[[k]]$time = end_time - start_time
    
    # Our method
	set.seed(999)
    start_time <- Sys.time()
    mBVVA[[k]] = BVVA2(X=sim[[k]]$X, y=sim[[k]]$y, iter = Niter, a1=2.01, a2=1, b1=2.01, b2=1, d1=1, d2=1, verbose=T)
    end_time <- Sys.time()
    mBVVA[[k]]$time = end_time - start_time
	
    # Our method
	set.seed(999)
	# set.seed(666)
    start_time <- Sys.time()
    hBVVA[[k]] = hBVVA2(X=sim[[k]]$X, y=sim[[k]]$y, iter = Niter, b1=2.01, b2=1, d1=1, d2=1, verbose=T)
    end_time <- Sys.time()
    hBVVA[[k]]$time = end_time - start_time
    }
    
    save(sim,file=paste0('Results/sim_',n,'_',p,'_',variance))
    save(mSVSS,file=paste0('Results/SVSS_samples_',n,'_',p,'_',variance))
	save(h_SVSS,file=paste0('Results/hSVSS_samples_',n,'_',p,'_',variance))
    save(mBVVA,file=paste0('Results/BVVA_samples_',n,'_',p,'_',variance))
	save(hBVVA,file=paste0('Results/hBVVA_samples_',n,'_',p,'_',variance))
  }
}