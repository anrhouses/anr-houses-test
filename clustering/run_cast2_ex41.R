library(CAST)
library(sf)
library(mapview)
library(caret)
library(lubridate)
library(ggplot2)
library(raster)
library(ranger)
library(scoringRules)
library(sp)
library(virtualspecies)
library(viridis)
library(latticeExtra)
library(gridExtra)
library(spatialsample)
library(gstat)
library(automap)

rm(list = ls())

source("utils.R")
source("sim_utils.R")

set.seed(15)

qq = seq(0,1,by=0.01)

#######################################################
## data
#######################################################
data("meuse")
# Read data
trainDat <- as.data.frame(meuse[,c("x","y","dist","ffreq","soil","zinc")])
names(trainDat)[ncol(trainDat)] = "response"
response = trainDat[,ncol(trainDat)]

train_points <- st_as_sf(meuse[,c("x","y")], coords = c("x", "y"), crs = 28992, remove = F)

data("meuse.area")
meuse.area <- SpatialPolygons(list(Polygons(list(Polygon(meuse.area)), "area")))
spoly <- st_as_sf(meuse.area)
st_crs(spoly) <- 28992

data("meuse.grid")
wclim <- st_as_sf(meuse.grid, coords = c("x", "y"), crs = 28992, remove = F)

'
base_cntrl <- trainControl(method="none")
base_mod <- train(response~., data=trainDat, method="rf",ntree=100)
train_points$res <- trainDat$response - predict(base_mod)
empvar_res <- variogram(as.formula("res~1"), cutoff=600000, train_points)
fitvar_res <- quiet(suppressWarnings(
  fit.variogram(empvar_res, vgm(model="Sph", nugget=0), fit.sills=c(F,T), fit.method = 1)))
lrange_res <- fitvar_res$range[fitvar_res$model=="Sph"]
folds_bLOO_res <- bLOO(train_points, lrange_res, 0.5)
'

#######################################################
## valid
#######################################################
NCV = 10
B = 250

resu = list()

VV = c("classic","LOL")
MOD = c("RF","KM")

for (mod in MOD){

for (vv in 1:length(VV)){

  if (VV[vv] == "LOL"){
    folds <- spatial_clustering_cv(train_points, v = NCV)
    #autoplot(folds)
    ## leave one location
    indices <- NULL
    indices$index <- list()
    indices$indexOut <- list()
    for (i in 1:NCV){
      indices$index[[i]] <- id0 <- folds$splits[[i]]$in_id
      indices$indexOut[[i]] <- (1:nrow(trainDat))[-id0]
    }
  }else{
    ## classical
    indices0 <- caret::createFolds(trainDat$response,k = NCV,returnTrain = TRUE)
    indices <- NULL
    indices$index <- indices0
    indices$indexOut <- list()
    for (i in 1:NCV){
      indices$indexOut[[i]] <- (1:nrow(trainDat))[-indices0[[i]]]
    }
  }

predB = randB = predB0 = list()
for (i in 1:NCV){
  
  trainDatB = trainDat[indices$index[[i]],]
  testDatB =trainDat[indices$indexOut[[i]],]
  ntest = nrow(testDatB)

  if (mod == "RF"){  
    modB = ranger(response ~., data = trainDatB, mtry = 3, quantreg = TRUE)
    ## quantile
    predB[[i]] = predict(modB, data = testDatB,
                      type = "quantiles",
                      quantiles = qq)
    ## random realisations
    rB = list()
    for (k in 1:ntest){
      U = runif(B)
      rB[[k]] = predict(modB, data = testDatB[k,],
                       type = "quantiles",
                       quantiles = U)$predictions
    }
    randB[[i]] = rB
    ## best estimate
    predB0[[i]] = predict(modB, data = testDatB)
  }else{
    coordinates(trainDatB) = ~x+y
    coordinates(testDatB) = ~x+y
    variogram = autofitVariogram(response ~ soil + ffreq + dist + x + y, trainDatB)
    #plot(variogram)
    rB = list()
    for (k in 1:ntest){
      rB[[k]] = krige(response ~ soil + ffreq + dist, trainDatB, newdata = testDatB[k,], 
                          model = variogram$var_model,
                          nsim = B)@data
    }
    randB[[i]] = rB
    ## best estimate
    tt = krige(response ~ soil + ffreq + dist, trainDatB, newdata = testDatB, 
               model = variogram$var_model)
    predB[[i]] = tt@data
    predB0[[i]] = tt@data$var1.pred
    #spplot( zn.condsim)
  }
  
    
} ## cv

#######################################################
## Indicators
#######################################################
alpha = seq(0,1,by=0.1)
Coverage = ES = VS = q2 = rmse = mae = CRPS = NULL
for (ii in 1:NCV){

  ## coverage
  if (mod == "RF"){ 
    pred = predB[[ii]]$predictions
    true = trainDat[indices$indexOut[[ii]],"response"]
    CC = CA(true,pred,alpha)
    Coverage[ii] = mean(abs(apply(CC,2,mean) - alpha))
  
    ## scores
    ntest = length(randB[[ii]])
    B = length(randB[[1]][[1]])
    real = matrix(0,ntest,B)
    for (k in 1:ntest){
      real[k,1:B] = as.vector(unlist(randB[[ii]][[k]]))
    }
    ES[ii] = es_sample(y = true, dat = real)
    VS[ii] = vs_sample(y = true, dat = real)
    CRPS[ii] = crps_sample(y = true, dat = real)
  
    ## Q2, rmse, mae
    q2[ii] = Q2(true, predB0[[ii]]$predictions)
    rmse[ii] = RMSE(true, predB0[[ii]]$predictions)
    mae[ii] = MAE(true, predB0[[ii]]$predictions)
  }else{
    pred = predB[[ii]]
    true = trainDat[indices$indexOut[[ii]],"response"]
    CC = CA.KM(true,pred,alpha)
    Coverage[ii] = mean(abs(apply(CC,2,mean) - alpha))
    
    ## scores
    ntest = length(randB[[ii]])
    B = length(randB[[1]][[1]])
    real = matrix(0,ntest,B)
    for (k in 1:ntest){
      real[k,1:B] = as.vector(unlist(randB[[ii]][[k]]))
    }
    ES[ii] = es_sample(y = true, dat = real)
    VS[ii] = vs_sample(y = true, dat = real)
    CRPS[ii] = crps_sample(y = true, dat = real)
    
    ## Q2, rmse, mae
    q2[ii] = Q2(true, predB0[[ii]])
    rmse[ii] = RMSE(true, predB0[[ii]])
    mae[ii] = MAE(true, predB0[[ii]])
  }
}

resu[[paste(VV[vv],mod)]] = data.frame(ES=ES,VS,Coverage,q2,rmse,mae,CRPS)

}## valid type
  
} ## mod

x11()
par(mfrow=c(1,2),mar=c(4,4,2,2))
df = data.frame(resu[[1]]$q2,resu[[2]]$q2,resu[[3]]$q2,resu[[4]]$q2)
names(df) = names(resu)
boxplot(df,main="q2")
df = data.frame(resu[[1]]$mae,resu[[2]]$mae,resu[[3]]$mae,resu[[4]]$mae)
names(df) = names(resu)
boxplot(df,main="mae")


x11()
par(mfrow=c(1,2),mar=c(4,4,2,2))
df = data.frame(resu[[1]]$Coverage,resu[[2]]$Coverage,resu[[3]]$Coverage,resu[[4]]$Coverage)
names(df) = names(resu)
boxplot(df,main="C")
df = data.frame(resu[[1]]$CRPS,resu[[2]]$CRPS,resu[[3]]$CRPS,resu[[4]]$CRPS)
names(df) = names(resu)
boxplot(df,main="CRPS")

x11()
par(mfrow=c(1,2),mar=c(4,4,2,2))
df = data.frame(resu[[1]]$ES,resu[[2]]$ES,resu[[3]]$ES,resu[[4]]$ES)
names(df) = names(resu)
boxplot(df,main="ES")
df = data.frame(resu[[1]]$VS,resu[[2]]$VS,resu[[3]]$VS,resu[[4]]$VS)
names(df) = names(resu)
boxplot(df,main="VS")


