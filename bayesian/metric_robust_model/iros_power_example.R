graphics.off()
rm(list=ls())
library(tidyr)
library(dplyr)

source("metric_robust_functions.R")


fileNameRoot = "power_testing"
dir.create(fileNameRoot)
fileNameRoot = paste0(fileNameRoot, "/", fileNameRoot)

ropeEffSz = c(-0.1, 0.1) # ROPE around null value of effect size

# Loading the Data In for the power_analysis
load("paired_rmse/paired_rmse-Mcmc.Rdata")
rmse <- codaSamples
mcmcMat_rmse = as.matrix(rmse,chains=TRUE)


load("paired_controleffort/paired_controleffort-Mcmc.Rdata")
ce <- codaSamples
mcmcMat_ce = as.matrix(ce,chains=TRUE)

load("paired_simtime/paired_simtime-Mcmc.Rdata")
st <- codaSamples
mcmcMat_st = as.matrix(st,chains=TRUE)

gName = "Controller"
yName = "RMSE"
groupNames = "FBL-QP - QP"
powerEstimation(mcmcChain = mcmcMat_rmse, muNULL = 0, effROPE = ropeEffSz, effHDImaxWid = 0.1, N = 500 ,nRep = 1000, yName = yName, gName = gName, groupNames = groupNames, saveName = fileNameRoot)