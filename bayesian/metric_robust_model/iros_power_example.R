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

load("paired_simtime/paired_simtime_summary_stats.RData")
st_median = 0.1*abs(summary_stats$diff_median[summary_stats$Controller=="QP"])

load("paired_rmse/paired_rmse_summary_stats.RData")
rmse_median = 0.1*abs(summary_stats$diff_median[summary_stats$Controller=="QP"])

load("paired_controleffort/paired_controleffort_summary_stats.RData")
ce_median = 0.1*abs(summary_stats$diff_median[summary_stats$Controller=="QP"])

effHDIMaxWid = 0.2;

gName = "Controller"
yName = "RMSE"
chain = mcmcMat_rmse
groupNames = "FBL-QP - QP"

MuHDImaxWid = 0.8*(2*rmse_median)

powerEstimation(mcmcChain = chain, muNULL = 0, effROPE = ropeEffSz, effHDImaxWid = effHDIMaxWid, MuHDImaxWid= MuHDImaxWid, N = 100 ,nRep = 1000, yName = yName, gName = gName, groupNames = groupNames, saveName = fileNameRoot)
powerEstimation(mcmcChain = chain, muNULL = 0, effROPE = ropeEffSz, effHDImaxWid = effHDIMaxWid, MuHDImaxWid= MuHDImaxWid, N = 200 ,nRep = 1000, yName = yName, gName = gName, groupNames = groupNames, saveName = fileNameRoot)
powerEstimation(mcmcChain = chain, muNULL = 0, effROPE = ropeEffSz, effHDImaxWid = effHDIMaxWid,MuHDImaxWid= MuHDImaxWid, N = 400 ,nRep = 1000, yName = yName, gName = gName, groupNames = groupNames, saveName = fileNameRoot)
powerEstimation(mcmcChain = chain, muNULL = 0, effROPE = ropeEffSz, effHDImaxWid = effHDIMaxWid,MuHDImaxWid= MuHDImaxWid, N = 800 ,nRep = 1000, yName = yName, gName = gName, groupNames = groupNames, saveName = fileNameRoot)
powerEstimation(mcmcChain = chain, muNULL = 0, effROPE = ropeEffSz, effHDImaxWid = effHDIMaxWid, MuHDImaxWid= MuHDImaxWid,N = 1000 ,nRep = 1000, yName = yName, gName = gName, groupNames = groupNames, saveName = fileNameRoot)

gName = "Controller"
yName = "SimTime"
chain = mcmcMat_st
groupNames = "FBL-QP - QP"
MuHDImaxWid = 0.8*(2*st_median)

powerEstimation(mcmcChain = chain, muNULL = 0, effROPE = ropeEffSz, effHDImaxWid = effHDIMaxWid,MuHDImaxWid= MuHDImaxWid, N = 100 ,nRep = 1000, yName = yName, gName = gName, groupNames = groupNames, saveName = fileNameRoot)
powerEstimation(mcmcChain = chain, muNULL = 0, effROPE = ropeEffSz, effHDImaxWid = effHDIMaxWid,MuHDImaxWid= MuHDImaxWid, N = 200 ,nRep = 1000, yName = yName, gName = gName, groupNames = groupNames, saveName = fileNameRoot)
powerEstimation(mcmcChain = chain, muNULL = 0, effROPE = ropeEffSz, effHDImaxWid = effHDIMaxWid,MuHDImaxWid= MuHDImaxWid, N = 400 ,nRep = 1000, yName = yName, gName = gName, groupNames = groupNames, saveName = fileNameRoot)
powerEstimation(mcmcChain = chain, muNULL = 0, effROPE = ropeEffSz, effHDImaxWid = effHDIMaxWid,MuHDImaxWid= MuHDImaxWid, N = 800 ,nRep = 1000, yName = yName, gName = gName, groupNames = groupNames, saveName = fileNameRoot)
powerEstimation(mcmcChain = chain, muNULL = 0, effROPE = ropeEffSz, effHDImaxWid = effHDIMaxWid,MuHDImaxWid= MuHDImaxWid, N = 1000 ,nRep = 1000, yName = yName, gName = gName, groupNames = groupNames, saveName = fileNameRoot)

gName = "Controller"
yName = "ControlEffort"
chain = mcmcMat_ce
groupNames = "FBL-QP - QP"
MuHDImaxWid = 0.8*(2*ce_median)

powerEstimation(mcmcChain = chain, muNULL = 0, effROPE = ropeEffSz, effHDImaxWid = effHDIMaxWid, MuHDImaxWid= MuHDImaxWid,N = 100 ,nRep = 1000, yName = yName, gName = gName, groupNames = groupNames, saveName = fileNameRoot)
powerEstimation(mcmcChain = chain, muNULL = 0, effROPE = ropeEffSz, effHDImaxWid = effHDIMaxWid, MuHDImaxWid= MuHDImaxWid,N = 200 ,nRep = 1000, yName = yName, gName = gName, groupNames = groupNames, saveName = fileNameRoot)
powerEstimation(mcmcChain = chain, muNULL = 0, effROPE = ropeEffSz, effHDImaxWid = effHDIMaxWid, MuHDImaxWid= MuHDImaxWid,N = 400 ,nRep = 1000, yName = yName, gName = gName, groupNames = groupNames, saveName = fileNameRoot)
powerEstimation(mcmcChain = chain, muNULL = 0, effROPE = ropeEffSz, effHDImaxWid = effHDIMaxWid, MuHDImaxWid= MuHDImaxWid,N = 800 ,nRep = 1000, yName = yName, gName = gName, groupNames = groupNames, saveName = fileNameRoot)
powerEstimation(mcmcChain = chain, muNULL = 0, effROPE = ropeEffSz, effHDImaxWid = effHDIMaxWid, MuHDImaxWid= MuHDImaxWid,N = 1000 ,nRep = 1000, yName = yName, gName = gName, groupNames = groupNames, saveName = fileNameRoot)