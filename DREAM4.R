
source("BVVA.R")
source("hBVVA.R")
source("SVSS.R")
source("hSVSS.R")
library(glmnet)

files = list.files('DREAM4')
HSVS = list()
hHSVS = list()
SVS = list()
hSVS = list()
lasso = list()
for (i in 1:5)
{
	print(i)
	dream = read.delim(paste('DREAM4/', files[i], sep = ''))
	
	HSVS[[i]] = list()
	hHSVS[[i]] = list()
	SVS[[i]] = list()
	hSVS[[i]] = list()
	lasso[[i]] = list()
	for (j in 1:100)
	{
		cat('-', j)
		
		set.seed(1010)
		HSVS[[i]][[j]] = BVVA2(10000, dream[, -j], dream[, j], a1=2.01, a2=1, b1=2.01, b2=1, d1=1, d2=1/2)

		set.seed(1010)
		hHSVS[[i]][[j]] = hBVVA2(10000, dream[, -j], dream[, j], b1=2.01, b2=1, d1=1, d2=1/2)

		set.seed(1010)
		SVS[[i]][[j]] = SVSS(10000, dream[, -j], dream[, j], a1=2.01, a2=1, b1=2.01, b2=1)

		set.seed(1010)
		hSVS[[i]][[j]] = hSVSS(10000, dream[, -j], dream[, j], b1=2.01, b2=1)

		lambda_seq <- 10^seq(2, -2, by = -.1)
		set.seed(1010)
		cv_output <- cv.glmnet(as.matrix(dream[, -j]), dream[, j],
							   family = 'gaussian', 
							   alpha = 1, lambda = lambda_seq, 
							   nfolds = 5)
		best_lam <- cv_output$lambda.min
		lasso[[i]][[j]] <- glmnet(as.matrix(dream[, -j]), dream[, j], family = 'gaussian', alpha = 1, lambda = best_lam)
	}
}

save(HSVS, hHSVS, SVS, hSVS, file = 'DREAM.RData')

log_loss = function(pred, actual, ep = 1e-15)
{
	pred = pmin(pmax(pred, ep), 1-ep)
	return(-sum(actual*log(pred) + (1-actual)*log(1-pred))/length(actual))
}

###TABLES and FIGURES
library(PRROC)
library(tibble)
library(xtable)
library(latex2exp)
files = list.files('DREAM4gold')
delta = 0.1
res = tibble(' ' = c(rep('log-loss', 5), rep('ROC', 5), rep('PR', 5)),
			 'DPSS' = rep(0.1, 15), 'DPHS' = rep(0.1, 15), 'SS' = rep(0.1, 15), 'HS' = rep(0.1, 15))
for (i in 1:length(files))
{
	gold = read.delim(paste('DREAM4gold/', files[i], sep = ''), header = F)
	# pdf(paste0("DREAM_",i,"_plots.pdf"), width = 12)
	pos = matrix(ncol = 4, nrow = 0)
	neg = matrix(ncol = 4, nrow = 0)
	all = matrix(ncol = 4, nrow = 0)
	class = rep(NA, 99*100)
	
	for (j in 1:100)
	{
		edges = gold[which(gold[[1]] == j & gold[[3]] == 1), 2]
		for (e in seq(along.with = edges))
			if (edges[e] > j)
				edges[e] = edges[e] - 1
		noedges = gold[which(gold[[1]] == j & gold[[3]] == 0), 2]
		for (e in seq(along.with = noedges))
			if (noedges[e] > j)
				noedges[e] = noedges[e] - 1
				
		class[(j-1)*99+edges] = 1
		class[(j-1)*99+noedges] = 0
		# class = rep(0, 99)
		# class[edges] = 1
		
		uno = apply(HSVS[[i]][[j]][-(1:5000), ], 2, mean)
		cuatro = apply(hHSVS[[i]][[j]][-(1:5000), ], 2, function(x) mean(abs(x) > delta))
		cinco = apply(SVS[[i]][[j]][-(1:5000), ], 2, mean)
		ocho = apply(hSVS[[i]][[j]][-(1:5000), ], 2, function(x) mean(abs(x) > delta))
		
		all = rbind(all, matrix(c(uno, cuatro, cinco, ocho), ncol = 4))
										  
		# uno = apply(HSVS[[i]][[j]]$beta[-(1:5000), ], 2, mean)
		# cuatro = apply(hHSVS[[i]][[j]]$beta[-(1:5000), ], 2, mean)
		# cinco = apply(SVS[[i]][[j]]$beta[-(1:5000), ], 2, mean)
		# ocho = apply(hSVS[[i]][[j]]$beta[-(1:5000), ], 2, mean)
								  
		# a = ggplot(data = tibble(ss = ocho, hss = cinco, dss = uno, hdss = cuatro, c = factor(class))) +
			# geom_point(aes(x = 1:99, y = ss, shape ="SS", color = c)) +
			# geom_point(aes(x = 1:99, y = hss, shape ="HS", color = c)) +
			# geom_point(aes(x = 1:99, y = dss, shape ="DPSS", color = c)) +
			# geom_point(aes(x = 1:99, y = hdss, shape ="DPHS", color = c)) +
			# scale_shape_manual(values=2:5) +
			# scale_color_manual(values = c('black', 'red')) +
			# theme(legend.position = "bottom") +
			# labs(title = TeX(paste0("Values for $\\beta$ parameters coded according to edges in the network.")), x = TeX("$\\beta$"), y = "", color = "Edge", shape = "Model")
		# print(a)
		# # readline('next')
	}
	# dev.off()
	
	if (any(is.na(class)))
		stop('fucked up')
	
	for (k in 1:ncol(all))
	{		
		res[c(i, 5+i, 10+i), k+1] = c(log_loss(all[, k], class), 
			roc.curve(scores.class0 = all[class == 1, k], scores.class1 = all[class == 0, k])$auc,
			pr.curve(scores.class0 = all[class == 1, k],, scores.class1 = all[class == 0, k])$auc.davis.goadrich)
	}
	# readline('next')
}
print(xtable(res, caption = "Classification errors in the five DREAM Gene Network detection datasets.", digits = 4), include.rownames=FALSE)

# > for(i in 1:5)
# + for(j in 1:100)
# + {print(quantile(apply(hHSVS[[i]][[j]]$sigma, 2, mean), probs = c(0.05, 0.1, 0.2, 0.8, 0.9, 0.95))); readline()}

#    ` `        DPSS   DPHS     SS     HS
#    <chr>     <dbl>  <dbl>  <dbl>  <dbl>
#  1 log-loss 0.0928 0.0873 0.114  0.0920
#  2 log-loss 0.138  0.120  0.150  0.125
#  3 log-loss 0.119  0.0914 0.161  0.120
#  4 log-loss 0.115  0.0964 0.159  0.112
#  5 log-loss 0.108  0.0939 0.131  0.104
#  6 ROC      0.616  0.590  0.605  0.604
#  7 ROC      0.599  0.580  0.601  0.609
#  8 ROC      0.644  0.675  0.681  0.692
#  9 ROC      0.645  0.693  0.659  0.672
# 10 ROC      0.597  0.621  0.633  0.660
# 11 PR       0.0874 0.0960 0.0933 0.101
# 12 PR       0.0577 0.0565 0.0708 0.0770
# 13 PR       0.123  0.124  0.126  0.139
# 14 PR       0.0942 0.114  0.0870 0.108
# 15 PR       0.0627 0.0856 0.0873 0.0939

# \begin{table}[ht]
# \centering
# \begin{tabular}{lrrrr}
#   \hline
#   & DPSS & DPHS & SS & HS \\
#   \hline
# log-loss & 0.0928 & 0.0873 & 0.1143 & 0.0920 \\
#   log-loss & 0.1382 & 0.1202 & 0.1498 & 0.1248 \\
#   log-loss & 0.1185 & 0.0914 & 0.1610 & 0.1196 \\
#   log-loss & 0.1152 & 0.0964 & 0.1588 & 0.1124 \\
#   log-loss & 0.1079 & 0.0939 & 0.1314 & 0.1044 \\
#   ROC & 0.6164 & 0.5901 & 0.6047 & 0.6038 \\
#   ROC & 0.5986 & 0.5798 & 0.6010 & 0.6089 \\
#   ROC & 0.6440 & 0.6746 & 0.6807 & 0.6920 \\
#   ROC & 0.6454 & 0.6934 & 0.6588 & 0.6724 \\
#   ROC & 0.5972 & 0.6208 & 0.6333 & 0.6605 \\
#   PR & 0.0874 & 0.0960 & 0.0933 & 0.1006 \\
#   PR & 0.0577 & 0.0565 & 0.0708 & 0.0770 \\
#   PR & 0.1229 & 0.1239 & 0.1258 & 0.1392 \\
#   PR & 0.0942 & 0.1139 & 0.0870 & 0.1082 \\
#   PR & 0.0627 & 0.0856 & 0.0873 & 0.0939 \\
#    \hline
# \end{tabular}
# \caption{Classification errors in the five DREAM Gene Network detection datasets.}
# \end{table}