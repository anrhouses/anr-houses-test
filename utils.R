require(scoringutils)

Q2<-function(y,yhat){
	return(1-sum((y-yhat)^2)/sum((y-mean(y))^2))
}

RMSE<-function(y,yhat){
	return(sqrt(mean((y-yhat)^2)))
}

MAE<-function(y,yhat){
	return((mean(abs(y-yhat))))
}

MAXE <- function(y,yhat){
	return(max(abs(y-yhat)))
}

CA <- function(y,lo,up){
	id <- NULL
	for (i in 1:length(y)){
		id[i] <- ifelse(y[i] >= lo[i] & y[i] <= up[i], 1, 0)
	}
	sum(id)/length(y)
}

CRPS <- function(Y.te,RSAMPLE){
	crps = crpsDecomposition(obs=Y.te, eps=RSAMPLE)
}

scores_fct = function(rsample,ytrue,conf){

	yhat = apply(rsample,1,mean)
	ytrue = as.vector(unlist(ytrue))

	## error metrics
	mae = MAE(ytrue,yhat)
	maxe = MAXE(ytrue,yhat)
	rmse = RMSE(ytrue,yhat)

	## prediction interval
	qq = seq(0,1,by=0.005)[-c(1,201)]
	qpred = t(apply(rsample,1,quantile,qq))
	lo = qpred[,which.min(abs(qq-(1-conf)/2))]
	up = qpred[,which.min(abs(qq-(1+conf)/2))]
	w.PI = mean((up-lo)/2)
	cov.PI = mean(sapply(ytrue,CA,lo,up))

	## Accuracy plot
	cov.pi = NULL
	qqlim = seq(0.05,0.95,by=0.05)
	for (i in 1:length(qqlim)){
		lo = qpred[,which.min(abs(qq-(1-qqlim[i])/2))]
		up = qpred[,which.min(abs(qq-(1+qqlim[i])/2))]
		cov.pi0 = sapply(ytrue,CA,lo,up)
		cov.pi[i] = mean(cov.pi0)
	}
	Mcov.PI = MAE(qqlim,cov.pi)

	## Crps
	crps0 <- crpsDecomposition(obs=ytrue, eps=rsample)
	crps <- crps0$CRPS
	crps_pot <- crps0$CRPSpot
	crps_reli <- crps0$Reli

	## Multivariate
	Energy <- es_sample(y=ytrue, dat=rsample)
	Variog <- vs_sample(y=ytrue, dat=rsample)
	#mmds <- mmds_sample(y=ytrue, dat=rsample)
	
	## New Multivariate --> correct packg Matrix for EPIT
	'
	ranki = get_prerank(ytrue, rsample, return_rank = F, prerank =  "average_rank")
	pit_hist(ranki,bins=101)
	evals <- epit::e_rank_histogram(r = r, h = lag, m = m,
                                  options = list(n0 = n0),
                                  strategy = strategy)$evalues_h
	evals <- epit:::evalue_combine_h(lapply(evals, function(x) x$e))
	'

	return(list=c(
			mae=mae,rmse=rmse,maxe=maxe,
			crps=crps,crps_pot=crps_pot,crps_reli=crps_reli,
			width=w.PI,coverage=cov.PI,accuracy_cov = Mcov.PI,
			Energy=Energy,Variog=Variog
			)
		)

}

#####################################################################
scores_fct2 = function(rsample,ytrue,conf){

	yhat = apply(rsample,1,mean)
	ytrue = as.vector(unlist(ytrue))


	## error metrics
	mae = MAE(ytrue,yhat)
	maxe = MAXE(ytrue,yhat)
	rmse = RMSE(ytrue,yhat)
	q2 = Q2(ytrue,yhat)

	## prediction interval
	qq = seq(0,1,by=0.005)[-c(1,201)]
	qpred = t(apply(rsample,1,quantile,qq))
	lo = qpred[,which.min(abs(qq-(1-conf)/2))]
	up = qpred[,which.min(abs(qq-(1+conf)/2))]

	int_score <- interval_score(
		true_values = ytrue,
		lower = lo,
		upper = up,
		interval_range = conf*100
		)
	w.PI = (up-lo)/2
	cov.PI = sapply(ytrue,CA,lo,up)

	q1_score <- quantile_score(ytrue,lo, (1-conf)/2, weigh = TRUE)
	q3_score <- quantile_score(ytrue,lo, (1+conf)/2, weigh = TRUE)

	## Accuracy plot
	cov.pi = NULL
	qqlim = seq(0.05,0.95,by=0.05)
	for (i in 1:length(qqlim)){
		lo = qpred[,which.min(abs(qq-(1-qqlim[i])/2))]
		up = qpred[,which.min(abs(qq-(1+qqlim[i])/2))]
		cov.pi0 = sapply(ytrue,CA,lo,up)
		cov.pi[i] = mean(cov.pi0)
	}
	Mcov.PI = MAE(qqlim,cov.pi)

	## Crps
	crps <- crps_sample(ytrue,rsample)
	dss <- dss_sample(ytrue,rsample)
	logs <- logs_sample(ytrue,rsample)
	dispersion <- mad_sample(rsample) ##dispersion

	## Multivariate
	Energy <- es_sample(y=ytrue, dat=rsample)
	Variog <- vs_sample(y=ytrue, dat=rsample)
	#mmds <- mmds_sample(y=ytrue, dat=rsample)

	SS = data.frame(
			mae = rep(mae,length(crps)),
			rmse = rep(rmse,length(crps)),
			q2 = rep(q2,length(crps)),
			maxe = rep(maxe,length(crps)),
			crps,
			logs,
			dss,
			dispersion,
			w.PI,
			cov.PI,
			int_score,
			q1_score,
			q3_score,
			Energy = rep(Energy,length(crps)),
			Variog = rep(Variog,length(crps))
			)

	return(list=c(
			SS
			)
		)

}

cluster_sample <- function(Ntr,Nce,tol){
	center = runif(1)
	X.ce <- maximinLHS(Nce,2)*(2*tol)+center - tol
	Nou = Ntr - Nce
	X.ou <- maximinLHS(Nou,2)
	XX = rbind(X.ce,X.ou)
	return(XX)
}

cluster_sample_mult <- function(Ntr,Nce,tol,mult){
	center = maximinLHS(mult,2)
	X.ce <- NULL
	for (ii in 1:mult){
		X.ce0 = maximinLHS(Nce,2)
		for (j in 1:2) X.ce0[,j] = X.ce0[,j]*(2*tol)+center[ii,j] - tol
		X.ce <- rbind(X.ce,X.ce0)
	}
	if (nrow(X.ce)<Ntr){
		Nou = Ntr - Nce
		X.ou <- maximinLHS(Nou,2)
		XX = rbind(X.ce,X.ou)
	}else{
		XX = X.ce
	}
	return(XX)
}

add_dist <- function(df.tr,mini,maxi){
	df.tr2 = df.tr
	corner = matrix(0,5,2)
	corner[1, ] = c(mini[1],mini[2])
	corner[2, ] = c(mini[1],maxi[2])
	corner[3, ] = c(maxi[1],mini[2])
	corner[4, ] = c(maxi[1],maxi[2])
	corner[5, ] = c((mini[1]+mini[2])/2,(maxi[1]+maxi[2])/2)

	mydist<-function(u,v=rep(0,2)){
		sqrt((u[1]-v[1])^2+(u[2]-v[2])^2)
	}
	DD = matrix(0,nrow(df.tr),5)
	for (ii in 1:5) DD[,ii] = apply(df.tr[,c("X1","X2")],1,mydist,v=corner[ii,])
	df.tr2$D1 = DD[,1]
	df.tr2$D2 = DD[,2]
	df.tr2$D3 = DD[,3]
	df.tr2$D4 = DD[,4]
	df.tr2$D5 = DD[,5]
	return(df.tr2)
}


mydist<-function(u,v=rep(0,2)){
		sqrt((u[1]-v[1])^2+(u[2]-v[2])^2)
}

add_dist2 <- function(df.te,df.tr,mydist=NULL){
	nte = nrow(df.te)
	ntr = nrow(df.tr)
	DD = matrix(0,nte,ntr)
	for (ii in 1:ntr){
		DD[,ii] = unlist(apply(df.te[,c("X1","X2")],1,mydist,v=df.tr[ii,c("X1","X2")]))
	}
	df.te2 = data.frame(df.te,DD)
	names(df.te2) = c(names(df.te),paste0("D",1:ntr))
	return(df.te2)
}

