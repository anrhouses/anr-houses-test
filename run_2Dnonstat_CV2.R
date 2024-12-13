library(lhs)
library(dgpsi)
library(ggplot2)
library(plotly)
library(tidyverse)
library(ranger)
library(arf)
library(scp)
library(verification)
library(scoringRules)
library(MultivCalibration)
library("WeightedForecastVerification")
library(deepgp)

source("utils.R")

###########################################################################
# construct a 2D non-stationary function that takes a matrix as the input
f <- function(x) {
sin(1/((0.7*x[,1,drop=F]+0.3)*(0.7*x[,2,drop=F]+0.3)))
}

Ntr = 10000
X <- expand.grid(seq(0,1,length.out=101),seq(0,1,length.out=101))#matrix(runif(2*Ntr),ncol=2)#maximinLHS(Ntr,2)
names(X) = c("X1","X2")
Y <- f(X)

df.plt = data.frame(X,Y=Y)
names(df.plt)[3] = "Y"
ggplot(df.plt,aes(X1,X2,color=Y))+geom_point()

p <- plot_ly(
  df.plt, x = ~X1, y = ~X2, z = ~Y, 
  color = ~Y
  ) %>%
  add_markers()

###########################################################################
# DATA TRAINING
Ntr = 100
set.seed(12345)

X.tr0 <- maximinLHS(Ntr,2)
X.tr0 <- cluster_sample_mult(Ntr,20,0.15,5)

Y.tr0 <- f(X.tr0)
df.tr0 = data.frame(X.tr0,Y=Y.tr0)

df.plt = data.frame(X=X.tr0,Y=Y.tr0)
ggplot(df.plt,aes(X.1,X.2,color=Y))+geom_point()

NCV = 10
nmc = 250

# random
ncv = list()
ncv0 = matrix(sample(1:Ntr,Ntr),ncol=NCV)
for (ii in 1:NCV) ncv[[ii]] = ncv0[,ii] 

# kmeans

cl = kmeans(X.tr0,NCV)
for (ii in 1:NCV) ncv[[ii]] = which(cl$cluster ==ii)


scores = resamples = list()
for (icv in 1:1){

	X.tr <- X.tr0[-ncv[[icv]],]
	Y.tr <- Y.tr0[-ncv[[icv]]]
	df.tr = data.frame(X1=X.tr[,1],X2=X.tr[,2],Y=Y.tr)

	X.te <- X.tr0[ncv[[icv]],]
	Y.te <- Y.tr0[ncv[[icv]]]
	df.te = data.frame(X1=X.te[,1],X2=X.te[,2],Y=Y.te)

###########################################################################
# GP
	'
	mod <- gp(X.tr, Y.tr)
	RSAMPLE = predict(object=mod,x=X.te,method="sampling",sample_size=nmc)$results
	scores[["gp"]] = rbind(scores[["gp"]],scores_fct(RSAMPLE,ytrue=Y.te,0.95))
	'
	
	fit <- fit_one_layer(x=X.tr, y=Y.tr, nmcmc = 1000)
	#fit <- trim(fit, 500, 2)
	pred = predict(fit, X.te, cores = 1, lite=FALSE)
      Sigma_smooth <- pred$Sigma - diag(mean(pred$g * pred$tau2), nrow(X.te))
      RSAMPLE <- t(mvtnorm::rmvnorm(nmc, pred$mean, Sigma_smooth, checkSymmetry = F))
	resamples[["gp"]] = rbind(resamples[["gp"]],RSAMPLE)
	

###########################################################################
# TWO LAYER GP
	'
	mod <- dgp(X.tr, Y.tr, B = nmc/50)
	RSAMPLE = predict(object=mod,x=X.te,method="sampling",sample_size=50)$results$output1
	scores[["dgp"]] = rbind(scores[["dgp"]],scores_fct(RSAMPLE,ytrue=Y.te,0.95))
	'
	
	fit <- fit_two_layer(x=X.tr, y=Y.tr, nmcmc = 1000)
	#fit <- trim(fit, 500, 2)
	pred = predict(fit, X.te, cores = 1, lite=FALSE)
      Sigma_smooth <- pred$Sigma - diag(mean(pred$g * pred$tau2), nrow(X.te))
      RSAMPLE <- t(mvtnorm::rmvnorm(nmc, pred$mean, Sigma_smooth, checkSymmetry = F))
	resamples[["dgp"]] = rbind(resamples[["dgp"]],RSAMPLE)
	

###########################################################################
# QRF
	mod <- ranger(Y~.,df.tr, keep.inbag = TRUE, quantreg=T, num.trees=1000)
	mod.rg <- ranger(Y~.,df.tr, num.trees=1000)
	qq = seq(0,1,length.out=nmc)
	RSAMPLE = unlist(predict(mod, df.te, type = "quantiles", quantiles = qq)$predictions)
	resamples[["qrf"]] = rbind(resamples[["qrf"]],RSAMPLE)

###########################################################################
# QRF - spatial
	df.tr2 = add_dist(df.tr)
	### ADD distances
	mod <- mod.rg <- ranger(Y~.,df.tr2, keep.inbag = TRUE, quantreg=T, num.trees=1000)
	qq = seq(0,1,length.out=nmc)
	df.te2 = add_dist(df.te)

	RSAMPLE = unlist(predict(mod, df.te2, type = "quantiles", quantiles = qq)$predictions)
	resamples[["qrf_spa"]] = rbind(resamples[["qrf_spa"]],RSAMPLE)

###########################################################################
# CONFORMAL
	'
	funRG = function(s0,s,Y,model=mod.rg){
		yhat.rg = predict(object=model,data=data.frame(X1=s0[1],X2=s0[2]) )$predictions
		return(yhat.rg)
	}

	s0 = X.tr[1,]
	system.time(mod.conf <- scp(s0=X.te[1:10,],s=X.tr,Y=Y.tr,pred_fun=funRG, global=TRUE,alpha=0.1))
	mod.conf = scp(s0=X.te[1:2,],s=X.tr,Y=Y.tr,pred_fun=funRG, global=FALSE,m=50,alpha=0.05)
	'

	#pit = pit_sample(Y.te, RSAMPLE)


###########################################################################
# GENERATIVE
	arf_mod <- adversarial_rf(df.tr, delta = 0, num_trees = 250L)

	# Estimate parameters
	params_mod <- forde(arf_mod, df.tr)

	# Evidence
	RSAMPLE = matrix(0,nrow(X.te),nmc)
	for (ii in 1:nrow(X.te)){
		evi <- data.frame(X1=X.te[ii,1],X2=X.te[ii,2])
		RSAMPLE[ii,] <- forge(params_mod, n_synth = nmc, evidence = evi)[,3]
	}
	resamples[["arf"]] = rbind(resamples[["arf"]],RSAMPLE)

}
###########################################################################
# POST - TRAIT
Ytrue = Y.tr0[unlist(ncv[[1]])]
for (ii in 1:length(resamples)){
	scores[[ii]] = scores_fct2(resamples[[ii]],Ytrue,.95)
}

###########################################################################
# PLOT
df.plt = data.frame(
		method =  NULL,
		mae =  NULL,
		maxe = NULL,
		crps = NULL,
		crps_reli = NULL,
		crps_pot = NULL,
		width = NULL,
		coverage = NULL,
		accuracy_cov = NULL,
		energy = NULL,
		variog = NULL
)

for (ii in 1:length(scores)){
	df.plt0 = data.frame(
		method = rep(names(resamples)[ii],1),
		mae =  scores[[ii]]["mae"],
		maxe = scores[[ii]]["maxe"],
		crps = scores[[ii]]["crps"],
		crps_reli = scores[[ii]]["crps_reli"],
		crps_pot = scores[[ii]]["crps_pot"],
		width = scores[[ii]]["width"],
		coverage = scores[[ii]]["coverage"],
		accuracy_cov = scores[[ii]]["accuracy_cov"],
		energy = scores[[ii]]["Energy"],
		variog = scores[[ii]]["Variog"]
	)
	df.plt = rbind(df.plt,df.plt0)
}

ggplot(df.plt,aes(as.factor(method),crps_reli))+geom_bar(stat="identity")
x11();ggplot(df00,aes(as.factor(method),crps_reli))+geom_bar(stat="identity")


