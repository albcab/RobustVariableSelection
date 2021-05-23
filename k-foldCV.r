#k-fold Cross Validation of predictive accracy

source("BVVA.R")
source("hBVVA.R")
source("SVSS.R")
source("hSVSS.R")

#NP density calc.
dens = function(bts, d, y, X, a, b) 
{
    t = length(bts)
    uniques = unique(bts[(d+2):t])
    n = c(counts(bts[(d+2):t], uniques), bts[d+1])
    n = n/sum(n)
    dens = c(dnorm(y, sum(X*bts[1:d]), sqrt(uniques)), b^a/sqrt(2*pi)*gamma(a+1/2)/gamma(a)*((y-sum(X*bts[1:d]))^2/2+b)^(-a-1/2))
    return(sum(n*dens))
}

log_pred_HSVS = c(); log_pred_hHSVS = c(); log_pred_SVS = c(); log_pred_hSVS = c()
CI_HSVS = c(); CI_hHSVS = c(); CI_SVS = c(); CI_hSVS = c()

K = 10 #k-fold

###PREPARE DATASETS
file_names = c('pollution', 'diabetes', 'baseball', 'temp')
train_size = c(55, 400, 300, 1800)
y_var = c(16, 11, 1, 82)
#pollution
pollutionp = read.table("pollution.data.txt", header = T, strip.white = T)
#diabetes
diabetesp = read.delim("diabetes.tab.txt")
diabetesp$SEX = as.numeric(diabetesp$SEX == 1)
#baseball
baseballp = read.csv("baseball.txt")
#temp
temp = read.csv('train.csv')
set.seed(123)
temp = temp[sample(nrow(temp), 2000, replace = F), ]
#all data
data = list(pollutionp, diabetesp, baseballp, temp)

for (D in 1:length(data))
{
set.seed(1000)
n_row = nrow(data[[D]])
train = lapply(1:K, function(k) sample(n_row, train_size[D], replace = F))
HSVS = vector(mode = 'list', length = K)
hHSVS = HSVS; SVS = HSVS; hSVS = HSVS

ll_HSVS = rep(0, K); ll_hHSVS = rep(0, K); ll_SVS = rep(0, K); ll_hSVS = rep(0, K)
gc()
for (k in 1:K)
{
    print(paste(file_names[D], "at iteration", k, "of CV"))
    
    data_train = data[[D]][train[[k]], ]
    data_test = data[[D]][-train[[k]], ]

    set.seed(100)
    HSVS[[k]] = BVVA2(10000, data_train[, -y_var[D]], log(data_train[, y_var[D]]), a1=2.01, a2=1, b1=2.01, b2=1, d1=0, d2=0, verbose=T, scaling = F)

    set.seed(100)
    hHSVS[[k]] = hBVVA2(10000, data_train[, -y_var[D]], log(data_train[, y_var[D]]), b1=2.01, b2=1, d1=0, d2=0, verbose=T, scaling = F)

    set.seed(100)
    SVS[[k]] = SVSS(10000, data_train[, -y_var[D]], log(data_train[, y_var[D]]), a1=2.01, a2=1, b1=2.01, b2=1, verbose=T, scaling = F)

    set.seed(100)
    hSVS[[k]] = hSVSS(10000, data_train[, -y_var[D]], log(data_train[, y_var[D]]), b1=2.01, b2=1, verbose=T, scaling = F)

    where = -c(1:5000)
    y = scale(log(data_test[, y_var[D]]), scale = F)
    X = scale(data_test[, -y_var[D]])
    d = ncol(X)
    
    for (i in 1:length(y))
    {
        cat('\r', i)
        ll_HSVS[k] = ll_HSVS[k] + log(mean(apply(cbind(HSVS[[k]]$beta, HSVS[[k]]$alpha, HSVS[[k]]$sigma)[where, ], 1, function(b) dens(b, d, y[i], X[i, ], 2.01, 1))))
        ll_hHSVS[k] = ll_hHSVS[k] + log(mean(apply(cbind(hHSVS[[k]]$beta, hHSVS[[k]]$alpha, hHSVS[[k]]$sigma)[where, ], 1, function(b) dens(b, d, y[i], X[i, ], 2.01, 1))))
        ll_SVS[k] = ll_SVS[k] + log(mean(apply(cbind(SVS[[k]]$beta, SVS[[k]]$sigma)[where, ], 1, function(bs) dnorm(y[i], sum(X[i, ]*bs[1:d]), sqrt(bs[d+1])))))
        ll_hSVS[k] = ll_hSVS[k] + log(mean(apply(cbind(hSVS[[k]]$beta, hSVS[[k]]$sigma)[where, ], 1, function(bs) dnorm(y[i], sum(X[i, ]*bs[1:d]), sqrt(bs[d+1])))))
    }
}

save(HSVS, hHSVS, SVS, hSVS, data_train, data_test, ll_HSVS, ll_hHSVS, ll_SVS, ll_hSVS, file = paste0("Results/", file_names[D], "_", ".RData"))

log_pred_HSVS = c(log_pred_HSVS, mean(ll_HSVS, na.rm = T))
log_pred_hHSVS = c(log_pred_hHSVS, mean(ll_hHSVS, na.rm = T))
log_pred_SVS = c(log_pred_SVS, mean(ll_SVS, na.rm = T))
log_pred_hSVS = c(log_pred_hSVS, mean(ll_hSVS, na.rm = T))

CI_HSVS = c(CI_HSVS, paste0('[', paste(round(quantile(ll_HSVS, probs = c(0.05, 0.95), na.rm = T), digits = 2), collapse = ','), ']'))
CI_hHSVS = c(CI_hHSVS, paste0('[', paste(round(quantile(ll_hHSVS, probs = c(0.05, 0.95), na.rm = T), digits = 2), collapse = ','), ']'))
CI_SVS = c(CI_SVS, paste0('[', paste(round(quantile(ll_SVS, probs = c(0.05, 0.95), na.rm = T), digits = 2), collapse = ','), ']'))
CI_hSVS = c(CI_hSVS, paste0('[', paste(round(quantile(ll_hSVS, probs = c(0.05, 0.95), na.rm = T), digits = 2), collapse = ','), ']'))
}
