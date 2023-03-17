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
set.seed(10)
NCV = 10
B = 100

resu = list()

VV = c("classic","LOL")

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
}

#######################################################
## Indicators
#######################################################
alpha = seq(0,1,by=0.1)
Coverage = ES = VS = q2 = rmse = mae = NULL
for (ii in 1:NCV){

  ## coverage
  pred = predB[[ii]]$predictions
  true = trainDat[indices$indexOut[[ii]],"response"]
  CC = CA(true,pred,alpha)
  Coverage[ii] = mean(abs(apply(CC,2,mean) - alpha))

  ## scores
  ntest = length(randB[[ii]])
  B = length(randB[[1]][[1]])
  real = matrix(0,ntest,B)
  for (k in 1:ntest){
    real[k,] = randB[[ii]][[k]]
  }
  ES[ii] = es_sample(y = true, dat = real)
  VS[ii] = vs_sample(y = true, dat = real)

  ## Q2, rmse, mae
  q2[ii] = Q2(true, predB0[[ii]]$predictions)
  rmse[ii] = RMSE(true, predB0[[ii]]$predictions)
  mae[ii] = MAE(true, predB0[[ii]]$predictions)

}

resu[[VV[vv]]] = data.frame(ES=ES,VS,Coverage,q2,rmse,mae)

}## valid type


par(mfrow=c(2,2),mar=c(4,4,2,2))
boxplot(cbind(resu[[1]]$q2,resu[[2]]$q2),main="q2")
boxplot(cbind(resu[[1]]$ES,resu[[2]]$ES),main="ES")
boxplot(cbind(resu[[1]]$VS,resu[[2]]$VS),main="VS")
boxplot(cbind(resu[[1]]$Coverage,resu[[2]]$Coverage),main="C")

