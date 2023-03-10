library(CAST)
library(sf)
library(mapview)
library(caret)
library(lubridate)
library(ggplot2)
library(raster)
library(ranger)
library(scoringRules)

rm(list = ls())

source("utils.R")

#######################################################
## data
#######################################################
data <- get(load(system.file("extdata","Cookfarm.RData",package="CAST")))

#######################################################
## train
#######################################################
qq = seq(0,1,by=0.01)
trainDat <- data[data$altitude==-0.3&
                   year(data$Date)==2012&
                   week(data$Date)%in%c(10:12),]

predictors <- c("DEM","TWI","Precip_cum","cday",
                "MaxT_wrcc","Precip_wrcc","BLD",
                "Northing","Easting","NDRE.M")

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
  ## leave one location
    indices <- CreateSpacetimeFolds(trainDat,spacevar = "SOURCEID",k = NCV)
  }else{
    ## classical
    indices0 <- caret::createFolds(trainDat$VW,k = NCV,returnTrain = TRUE)
    indices <- NULL
    indices$index <- indices0
    indices$indexOut <- list()
    for (i in 1:NCV){
      indices$indexOut[[i]] <- (1:nrow(trainDat))[-indices0[[i]]]
    }
  }

predB = randB = predB0 = list()
for (i in 1:NCV){
  trainDatB = trainDat[indices$index[[i]],c(predictors,"VW")]
  testDatB =trainDat[indices$indexOut[[i]],predictors]
  ntest = nrow(testDatB)
  modB = ranger(VW ~., data = trainDatB, mtry = 2, quantreg = TRUE)

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
  true = trainDat[indices$indexOut[[ii]],"VW"]
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