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
library(akima)

rm(list=ls())

source("utils.R")

###########################################################################
# construct a 2D non-stationary function that takes a matrix as the input
x100 = seq(10,870,by=10)
x200 = seq(10,610,by=10)

f <- function(x,x0=x100,y0=x200) {
	#x <- matrix(x, nrow = 1)
	res = bilinear(x=x0, y=y0, z=volcano, x0=x[,1], y0=x[,2])$z
	res = (res - mean(volcano))/sd(volcano)
}

x10 = seq(10,870,by=10)
x20 = seq(10,610,by=10)

X.te <- X.te_uns <- matrix(0,length(x10)*length(x20),2)
c = 0
for (i in 1:length(x10)){
	for (j in 1:length(x20)){
		c = c + 1
		X.te[c,1] = x10[i]
		X.te[c,2] = x20[j]
	}
}

X.te_uns[,1] = (X.te[,1] - 10) /(870-10)
X.te_uns[,2] = (X.te[,2] -10)/ (610-10)

Y.te <- f(X.te,x0=x100,y0=x200)
df.te = data.frame(X1=X.te[,1],X2=X.te[,2],Y=Y.te)

df.plt = data.frame(X.te,Y=Y.te)
names(df.plt)[3] = "Y"
ggplot(df.plt,aes(X1,X2,color=Y))+geom_point()+scale_colour_viridis_c(direction=1)

p <- plot_ly(
  df.plt, x = ~X1, y = ~X2, z = ~Y, 
  color = ~Y
  ) %>%
  add_markers()
print(p)

###########################################################################
# DATA TRAINING
set.seed(12345)

Ntr_list = c(50,100)
Nit = 25
Ncv = 5
nmc = 1000

for (SAMPL in c("unif-unif","clust-unif","clust-clust")){

for (Ntr in Ntr_list){

for (it in 1:Nit){

if (SAMPL == "unif-unif") X.tr_uns <- maximinLHS(Ntr,2)
if (SAMPL == "clust-unif" | SAMPL == "clust-clust") X.tr_uns <- cluster_sample_mult(Ntr,20,0.15,5)

X.tr0_uns = X.tr0 = X.tr_uns
X.tr0[,1] = X.tr0[,1]*(870-10)+10
X.tr0[,2] = X.tr0[,2]*(610-10)+10

Y.tr0 <- f(X.tr0,x0=x100,y0=x200)
df.tr0 = data.frame(X1=X.tr0_uns[,1],X2=X.tr0_uns[,2],Y=Y.tr0)

### CV
set.seed(12345)
ncv = list()
# random
if (SAMPL == "unif-unif" | SAMPL == "clust-unif"){
	ncv0 = matrix(sample(1:Ntr,Ntr),ncol=Ncv)
	for (ii in 1:Ncv){ncv[[ii]] = ncv0[,ii]}
}
# kmeans
if (SAMPL == "clust-clust"){
	cl = kmeans(X.tr0_uns,Ncv)
	for (ii in 1:Ncv){ncv[[ii]] = which(cl$cluster ==ii)}
}

for (icv in 1:Ncv){

	X.tr_uns <- X.tr0_uns[-ncv[[icv]],]
	Y.tr <- Y.tr0[-ncv[[icv]]]
	df.tr = data.frame(X1=X.tr_uns[,1],X2=X.tr_uns[,2],Y=Y.tr)

	X.te <- X.tr0[ncv[[icv]],]
	X.te_uns <- X.tr0_uns[ncv[[icv]],]
	Y.te <- Y.tr0[ncv[[icv]]]
	df.te = data.frame(X1=X.te_uns[,1],X2=X.te_uns[,2],Y=Y.te)

# random
scores = rsample = list()
###########################################################################
# GP
	mod <- dgp(X.tr_uns, Y.tr, depth=2, nugget_est = T, name = "matern2.5", struc = NULL, B = 50)
	RSAMPLE = predict(object=mod,x=X.te_uns,method="sampling",sample_size=nmc/50)$results$output1
	#rsample[["dgp"]] = (RSAMPLE[,sample(1:ncol(RSAMPLE),nmc)])
	rsample[["dgp"]] = (RSAMPLE)
	scores[["dgp"]] = rbind(scores_fct2(RSAMPLE,ytrue=Y.te,0.95))

	km = DiceKriging::km(~.,X.tr_uns,Y.tr,nugget.estim=T)
	RSAMPLE = DiceKriging::simulate(km, nsim=nmc, seed=NULL, newdata=X.te_uns,
						cond=TRUE, nugget.sim=0, checkNames=F)
	rsample[["km"]] = t(RSAMPLE)
	scores[["km"]] = scores_fct2(rsample[["km"]],ytrue=Y.te,0.95)

###########################################################################
# QRF
	mod <- mod.rg <- ranger(Y~.,df.tr, keep.inbag = TRUE, quantreg=T, num.trees=1000)
	qq = seq(0,1,length.out=nmc)
	RSAMPLE = unlist(predict(mod, df.te, type = "quantiles", quantiles = qq)$predictions)
	rsample[["qrf"]] = RSAMPLE
	scores[["qrf"]] = rbind(scores[["qrf"]],scores_fct2(RSAMPLE,ytrue=Y.te,0.95))

###########################################################################
# QRF - spatial
	df.tr2 = add_dist(df.tr,mini=c(10,10),maxi=c(870,610))
	### ADD distances
	mod <- mod.rg <- ranger(Y~.,df.tr2, keep.inbag = TRUE, quantreg=T, num.trees=1000)
	qq = seq(0,1,length.out=nmc)
	df.te2 = add_dist(df.te,mini=c(10,10),maxi=c(870,610))

	RSAMPLE = unlist(predict(mod, df.te2, type = "quantiles", quantiles = qq)$predictions)
	rsample[["qrf_spa"]] = RSAMPLE
	scores[["qrf_spa"]] = scores_fct2(RSAMPLE,ytrue=Y.te,0.95)

###########################################################################
# GENERATIVE
	arf_mod <- adversarial_rf(df.tr, delta = 0, num_trees = 250L)

	# Estimate parameters
	params_mod <- forde(arf_mod, df.tr)

	# Evidence
	RSAMPLE = matrix(0,nrow(X.te_uns),nmc)
	for (ii in 1:nrow(X.te_uns)){
		evi <- data.frame(X1=X.te_uns[ii,1],X2=X.te_uns[ii,2])
		RSAMPLE[ii,] <- forge(params_mod, n_synth = nmc, evidence = evi)[,3]
	}
	rsample[["arf"]] = RSAMPLE
	scores[["arf"]] = scores_fct2(RSAMPLE,ytrue=Y.te,0.95)

###########################################################################
# GENERATIVE - SPA
	arf_mod <- adversarial_rf(df.tr2, delta = 0, num_trees = 250L)

	# Estimate parameters
	params_mod <- forde(arf_mod, df.tr2)

	# Evidence
	RSAMPLE = matrix(0,nrow(X.te),nmc)
	for (ii in 1:nrow(X.te_uns)){
		evi <- df.te2[ii,-3]#data.frame(X1=X.te[ii,1],X2=X.te[ii,2])
		RSAMPLE[ii,] <- forge(params_mod, n_synth = nmc, evidence = evi)[,3]
	}
	rsample[["arf_spa"]] = RSAMPLE
	scores[["arf_spa"]] = scores_fct2(RSAMPLE,ytrue=Y.te,0.95)

###########################################################################
# PLOT

save(scores,df.te,Y.te,rsample,file=paste0("./cv/scores_Ntr",Ntr,"_it",it,"_cv",icv,"_sampl",SAMPL,".RData"))

}## CV

}## IT

}## Ntr

}## sample


###############################################################
#### POST TRAIT
###############################################################
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(median = mean(x[[col]], na.rm=TRUE)
      #q1 = quantile(x[[col]], na.rm=TRUE,0.25),
      #q3 = quantile(x[[col]], na.rm=TRUE,0.75)
	)
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  return(data_sum)
}


df.plt0 = NULL

Ntr0 = Sample0 = ss0 = NULL

for (it in 1:25){

print(it)

for (SAMPL in c("unif-unif","clust-unif","clust-clust")){

for (Ntr in c(50,100)){

df.plt00 = df.plt.init

for (icv in 1:5){

load(file=paste0("./cv/scores_Ntr",Ntr,"_it",it,"_cv",icv,"_sampl",SAMPL,".RData"))

for (ii in 1:length(scores)){
	ss = scores[[ii]]
	if (is.null(nrow(ss))==FALSE){
		ss = data.frame(apply(scores[[ii]],2,unlist))
		ss$method = rep(names(scores)[ii],nrow(ss))
	}else{
		ss = data.frame(ss)
		ss$method = rep(names(scores)[ii],nrow(ss))
	}
	ss0 = rbind(ss0,ss)
	Ntr0 = c(Ntr0,rep(Ntr,nrow(ss)))
	Sample0 = c(Sample0,rep(SAMPL,nrow(ss)))
}


}##icv

}## Ntr

}## sample

	df.plt = data.frame(ss0)
	df.plt$Ntr = Ntr0
	df.plt$Sampling = Sample0
	df.plt00 = NULL
	for (ii in 1:15) df.plt00 = cbind(df.plt00,data_summary(df.plt,colnames(df.plt)[ii],c("method","Ntr","Sampling"))$median)
	df.plt00 = data.frame(df.plt00)
	names(df.plt00) = colnames(df.plt)[1:15]
	df.plt00$method = data_summary(df.plt,colnames(df.plt)[1],c("method","Ntr","Sampling"))$method
	df.plt00$Ntr = data_summary(df.plt,colnames(df.plt)[1],c("method","Ntr","Sampling"))$Ntr
	df.plt00$Sampling = data_summary(df.plt,colnames(df.plt)[1],c("method","Ntr","Sampling"))$Sampling
	df.plt0 = rbind(df.plt0,df.plt00)
	

}##it


#################
## PLOT
tt0 = theme(
            legend.key.size = unit(0.75, 'cm'), #change legend key size
            #legend.key.height = unit(0.5, 'cm'), #change legend key height
            #legend.key.width = unit(0.5, 'cm'), #change legend key width
            legend.title = element_blank(), #change legend title font size
            #legend.text = element_text(size=8),
            text = element_text(size = 14),
            legend.position="top"
) #change legend text font size

x11();ggplot(df.plt0,aes(method,crps,color=as.factor(Ntr)))+geom_boxplot()+facet_wrap(~Sampling)


f = which(df.plt0$Sampling == "clust-clust")
plt = list()
Metric = c("mae","int_score","q1_score","q3_score","dispersion","crps","Energy","Variog")
Nom_Metric = c("mae [m]","interval score","q1 score","q3 score","dispersion","crps","energy","variogram score")
for (i in 1:length(Metric)){
	df.pltA = df.plt0[f,c("method","Ntr")]
	df.pltA$score =  df.plt0[f,Metric[i]]
	plt[[i]] = ggplot(df.pltA,aes(method,score,color=as.factor(Ntr)))+geom_boxplot()+theme_bw()+ylab(Nom_Metric[i])+tt0+xlab("")
}
ggpubr::ggarrange(plotlist=plt,ncol=4, nrow=2, common.legend = TRUE, legend="top")


