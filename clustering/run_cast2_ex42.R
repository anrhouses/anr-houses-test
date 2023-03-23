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

#set.seed(123)

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

#######################################################
## valid
#######################################################
NCV = 10
B = 100
R = 50

resu = list()

VV = c("classic","LOL")
MOD = c("RF","KM")

for (mod in MOD){

for (vv in 1:length(VV)){

for (rr in 1:R){# repetitions

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
      rB[[k]] = krige(response ~ soil + ffreq + dist + x + y, trainDatB, newdata = testDatB[k,], 
                          model = variogram$var_model,
                          nsim = B)@data
    }
    randB[[i]] = rB
    ## best estimate
    tt = krige(response ~ soil + ffreq + dist + x + y, trainDatB, newdata = testDatB, 
               model = variogram$var_model)
    predB[[i]] = tt@data
    predB0[[i]] = tt@data$var1.pred
  }

} ## cv

#######################################################
## Indicators
#######################################################
alpha = seq(0,1,by=0.1)
Coverage = ES = VS = q2 = rmse = mae = CRPS = NULL

real = CC = true = pred = NULL

for (ii in 1:NCV){

  ## coverage
  if (mod == "RF"){ 
    
    pred0 = predB[[ii]]$predictions
    true0 = trainDat[indices$indexOut[[ii]],"response"]
    CC0 = CA(true0,pred0,alpha)
    
    pred = c(pred,pred0)
    true = c(true,true0)
    CC = rbind(CC,CC0)

    ntest = length(randB[[ii]])
    B = length(randB[[1]][[1]])
    real0 = matrix(0,ntest,B)
    for (k in 1:ntest){
      real0[k,1:B] = as.vector(unlist(randB[[ii]][[k]]))
    }
    real = rbind(real,real0)
    
  }else{

    pred0 = predB0[[ii]]
    true0 = trainDat[indices$indexOut[[ii]],"response"]
    CC0 = CA.KM(true0,predB[[ii]],alpha)

    pred = c(pred,pred0)
    true = c(true,true0)
    CC = rbind(CC,CC0)
    
    ntest = length(randB[[ii]])
    B = length(randB[[1]][[1]])
    real0 = matrix(0,ntest,B)
    for (k in 1:ntest){
      real0[k,1:B] = as.vector(unlist(randB[[ii]][[k]]))
    }
    real = rbind(real,real0)
    
  }
}

ES = es_sample(y = true, dat = real)
VS = vs_sample(y = true, dat = real)
MMD = MMDfct(true, real)

CRPS = mean(crps_sample(y = true, dat = real))
LogS = mean(logs_sample(y = true, dat = real))
DSS = mean(dss_sample(y = true, dat = real))
## Q2, rmse, mae
q2 = Q2(true, pred)
rmse = RMSE(true, pred)
mae = MAE(true, pred)
Coverage = mean(abs(apply(CC,2,mean) - alpha))
#plot(alpha,apply(CC,2,mean))
#abline(0,1)

resu0 = data.frame(ES,VS,MMD,Coverage,q2,rmse,mae,CRPS,LogS,DSS)

resu[[paste(VV[vv],mod)]] = rbind(resu[[paste(VV[vv],mod)]], resu0)

} ## repeat

}## valid type
  
} ## mod

save(resu,file="Meuse_10CV_Rep50.RData")

x11()
par(mfrow=c(1,2),mar=c(4,4,2,2))
df = data.frame(resu[[1]]$rmse,resu[[3]]$rmse,resu[[2]]$rmse,resu[[4]]$rmse)
names(df) = names(resu)[c(1,3,2,4)]
boxplot(df,main="rmse")
df = data.frame(resu[[1]]$mae,resu[[3]]$mae,resu[[2]]$mae,resu[[4]]$mae)
names(df) = names(resu)[c(1,3,2,4)]
boxplot(df,main="mae")


x11()
par(mfrow=c(2,2),mar=c(4,4,2,2))
df = data.frame(resu[[1]]$Coverage,resu[[3]]$Coverage,resu[[2]]$Coverage,resu[[4]]$Coverage)
names(df) = names(resu)[c(1,3,2,4)]
boxplot(df,main="C")
df = data.frame(resu[[1]]$CRPS,resu[[3]]$CRPS,resu[[2]]$CRPS,resu[[4]]$CRPS)
names(df) = names(resu)[c(1,3,2,4)]
boxplot(df,main="CRPS")
df = data.frame(resu[[1]]$LogS,resu[[3]]$LogS,resu[[2]]$LogS,resu[[4]]$LogS)
names(df) = names(resu)[c(1,3,2,4)]
boxplot(df,main="logS")
df = data.frame(resu[[1]]$DSS,resu[[3]]$DSS,resu[[2]]$DSS,resu[[4]]$DSS)
names(df) = names(resu)[c(1,3,2,4)]
boxplot(df,main="DSS")

x11()
par(mfrow=c(1,2),mar=c(4,4,2,2))
df = data.frame(resu[[1]]$ES,resu[[3]]$ES,resu[[2]]$ES,resu[[4]]$ES)
names(df) = names(resu)[c(1,3,2,4)]
boxplot(df,main="ES")
df = data.frame(resu[[1]]$VS,resu[[3]]$VS,resu[[2]]$VS,resu[[4]]$VS)
names(df) = names(resu)[c(1,3,2,4)]
boxplot(df,main="VS")
#df = data.frame(resu[[1]]$MMD,resu[[2]]$MMD,resu[[3]]$MMD,resu[[4]]$MMD)
#names(df) = names(resu)
#boxplot(df,main="MMD")

