# Script to run over range of n and p
# Make sure you create the folder 'Results' in your working directory

source("BVVA_tstud.r")
source("hBVVA_tstud.r")
source("SVSS_tstud.r")
source("hSVSS_tstud.r")

source('GetSimulationData_tstud.R')

nlist <- c(10,20,50,100,200); 
# nlist = as.numeric(readline("n: "))
plist <- c(50,100, 200, 500, 1000, 2000);
# plist = as.numeric(readline("p: "))

rho <- 0
variance <- 3.1
Niter <- 10000
Nsam <- 5000
K = 100

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
    # mSVSS[[k]] = SVSS(X=sim[[k]]$X, y=sim[[k]]$y, iter = Niter, a1=2.01, a2=1, b1=2.01, b2=1, verbose=T)
    mSVSS_ = SVSS(X=sim[[k]]$X, y=sim[[k]]$y, iter = Niter, a1=2.01, a2=1, b1=2.01, b2=1, verbose=T)
    end_time <- Sys.time()
    mSVSS[[k]]$beta = apply(mSVSS_$beta[-(1:Nsam), ], 2, mean)
    mSVSS[[k]]$tau = apply(mSVSS_$tau[-(1:Nsam), ], 2, mean)
    mSVSS[[k]]$w2 = apply(mSVSS_$w2[-(1:Nsam), ], 2, mean)
    mSVSS[[k]]$sigma = mSVSS_$sigma[-(1:Nsam)]
    mSVSS[[k]]$time = end_time - start_time
    
    # horseshoe Method
    set.seed(999)
    start_time <- Sys.time()
    # h_SVSS[[k]] = hSVSS(X=sim[[k]]$X, y=sim[[k]]$y, iter = Niter, b1=2.01, b2=1, verbose=T)
    h_SVSS_ = hSVSS(X=sim[[k]]$X, y=sim[[k]]$y, iter = Niter, b1=2.01, b2=1, verbose=T)
    end_time <- Sys.time()
    h_SVSS[[k]]$beta = apply(h_SVSS_$beta[-(1:Nsam), ], 2, mean)
    h_SVSS[[k]]$lambda = apply(h_SVSS_$lambda[-(1:Nsam), ], 2, mean)
    h_SVSS[[k]]$tau = h_SVSS_$tau[-(1:Nsam)]
    h_SVSS[[k]]$sigma = h_SVSS_$sigma[-(1:Nsam)]
    h_SVSS[[k]]$time = end_time - start_time
    
    # Our method
    set.seed(999)
    start_time <- Sys.time()
    # mBVVA[[k]] = BVVA2(X=sim[[k]]$X, y=sim[[k]]$y, iter = Niter, a1=2.01, a2=1, b1=2.01, b2=1, d1=1, d2=1, verbose=T)
    mBVVA_ = BVVA2(X=sim[[k]]$X, y=sim[[k]]$y, iter = Niter, a1=2.01, a2=1, b1=2.01, b2=1, d1=1, d2=1, verbose=T)
    end_time <- Sys.time()
    mBVVA[[k]]$alpha = mBVVA_$alpha[-(1:Nsam)]
    mBVVA[[k]]$beta = apply(mBVVA_$beta[-(1:Nsam), ], 2, mean)
    mBVVA[[k]]$tau = apply(mBVVA_$tau[-(1:Nsam), ], 2, mean)
    mBVVA[[k]]$sigma = apply(mBVVA_$sigma[-(1:Nsam), ], 2, mean)
    mBVVA[[k]]$c = mBVVA_$c[-(1:Nsam), ]
    mBVVA[[k]]$K = mBVVA_$K[-(1:Nsam)]
    mBVVA[[k]]$w2 = apply(mBVVA_$w2[-(1:Nsam), ], 2, mean)
    mBVVA[[k]]$time = end_time - start_time
    
    # Our method
    set.seed(999)
    start_time <- Sys.time()
    # hBVVA[[k]] = hBVVA2(X=sim[[k]]$X, y=sim[[k]]$y, iter = Niter, b1=2.01, b2=1, d1=1, d2=1, verbose=T)
    hBVVA_ = hBVVA2(X=sim[[k]]$X, y=sim[[k]]$y, iter = Niter, b1=2.01, b2=1, d1=1, d2=1, verbose=T)
    end_time <- Sys.time()
    hBVVA[[k]]$alpha = hBVVA_$alpha[-(1:Nsam)]
    hBVVA[[k]]$beta = apply(hBVVA_$beta[-(1:Nsam), ], 2, mean)
    hBVVA[[k]]$lambda = apply(hBVVA_$lambda[-(1:Nsam), ], 2, mean)
    hBVVA[[k]]$tau = hBVVA_$tau[-(1:Nsam)]
    hBVVA[[k]]$sigma = apply(hBVVA_$sigma[-(1:Nsam), ], 2, mean)
    hBVVA[[k]]$c = hBVVA_$c[-(1:Nsam), ]
    hBVVA[[k]]$K = hBVVA_$K[-(1:Nsam)]
    hBVVA[[k]]$time = end_time - start_time

    rm(mSVSS_, h_SVSS_, mBVVA_, hBVVA_)
    gc()
    }
    
    saveRDS(sim,file=paste0('Results/sim_',n,'_',p,'_',variance,'.rds'), compress='xz')
    saveRDS(mSVSS,file=paste0('Results/SVSS_samples_',n,'_',p,'_',variance,'.rds'), compress='xz')
    saveRDS(h_SVSS,file=paste0('Results/hSVSS_samples_',n,'_',p,'_',variance,'.rds'), compress='xz')
    saveRDS(mBVVA,file=paste0('Results/BVVA_samples_',n,'_',p,'_',variance,'.rds'), compress='xz')
    saveRDS(hBVVA,file=paste0('Results/hBVVA_samples_',n,'_',p,'_',variance,'.rds'), compress='xz')

    rm(sim, mSVSS, h_SVSS, mBVVA, hBVVA)
    gc()
  }
}