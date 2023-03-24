Cov <- function(x,xmin,xmax){
  ifelse(x >= xmin & x <= xmax,1,0)
}

extract <- function(pred){
  alpha = na.omit(as.numeric(unlist(strsplit(colnames(predB[[i]]$predictions),"quantile= "))))
  return(alpha)
}

CA <- function(true,pred,alpha){
  ntest = length(true)
  nalpha = length(alpha)
  ll = extract(pred)
  
  cov = matrix(0,ntest,nalpha)
  for (i in 1:ntest){
    for (k in 1:nalpha){
      ou0 = which.min(abs(ll - (1-alpha[k])/2))
      ou1 = which.min(abs(ll - (1+alpha[k])/2))                    
      pred0 = pred[i,ou0]
      pred1 = pred[i,ou1]
      cov[i,k] = Cov(true[i],pred0,pred1)
    }
  }
  return(cov)
}

CA.KM <- function(true,pred,alpha){
  ntest = length(true)
  nalpha = length(alpha)
  mu = pred$var1.pred
  s = sqrt(pred$var1.var)
  
  cov = matrix(0,ntest,nalpha)
  for (i in 1:ntest){
    for (k in 1:nalpha){
      pred0 = mu[i] - qnorm((alpha[k]+1)/2)*s[i]
      pred1 = mu[i] + qnorm((alpha[k]+1)/2)*s[i]
      cov[i,k] = Cov(true[i],pred0,pred1)
    }
  }
  return(cov)
}

CA.sampl <- function(true,real,alpha){
  ntest = length(true)
  nalpha = length(alpha)
  cov = matrix(0,ntest,nalpha)
  for (i in 1:ntest){
    for (k in 1:nalpha){
      pred0 = quantile(real[i,],abs(alpha[k]-1)/2)
      pred1 = quantile(real[i,],(alpha[k]+1)/2)
      cov[i,k] = Cov(true[i],pred0,pred1)
    }
  }
  return(cov)
}

Q2 <- function(true,pred){
  1 - sum((true - pred)^2) / sum((true - mean(true))^2) 
}

RMSE <- function(true,pred){
  sqrt(mean((true - pred)^2))
}

MAE <- function(true,pred){
  mean(abs(true - pred))
}

MMDfct <- function(true,real){
  real0 = t(real)
  true = matrix(true,nrow=1,ncol=ncol(real0))
  dmat <- as.matrix(dist(rbind(real0, true)))  # Euclidean distance matrix
  kmat <- exp(-(dmat^2) / mean(dmat)^2)                      # build a gaussian kernel matrix
  label  <- c(rep(1,nrow(real0)), rep(2,1))
  ulabel = unique(label)
  id1 = which(label==ulabel[1]); m = length(id1)
  id2 = which(label==ulabel[2]); n = length(id2)
  XX=kmat[id1,id1]
  YY=kmat[id2,id2]
  XY=kmat[id1,id2]
  m = nrow(XX)
  n = length(YY)
  mmd = (sum(XX)/(m^2)) + (sum(YY)/(n^2)) - ((2/(m*n))*sum(XY))
  return(mmd)
}