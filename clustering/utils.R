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

Q2 <- function(true,pred){
  1 - sum((true - pred)^2) / sum((true - mean(true))^2) 
}

RMSE <- function(true,pred){
  sqrt(mean((true - pred)^2))
}

MAE <- function(true,pred){
  mean(abs(true - pred))
}