graphics.off()
rm(list=ls())

source("metric_robust_functions.R")
source("lognormal_functions.R")





datFrm = read.csv(file="iros_data_6.csv")

datFrm = subset(datFrm, 
                Solved == 1 & 
                  IsStable == 1 & 
                  ConstraintsOK == 1 & JointsOK == 1 & SimTime <= 100)

gName = "Controller"
yName = "SimTime"

# Defining directory for results and prefix of the files names:
fileNameRoot = "iros_test_6_lnorm_ce_test"
dir.create(fileNameRoot)
fileNameRoot = paste0(fileNameRoot, "/", fileNameRoot)
graphFileType = "pdf"

analysis = "independent"
# ======= Independent/two-groups analysis =====================================
if (analysis == "independent") {
  # Setting the comparison values and ROPE intervals:
  compValMu = NULL # comparison value to show in plots of means
  compValMuDiff = NULL # comparison value to show in plot of means difference
  compValSigma = NULL # comparison value to show in plots of scales
  compValSigmaDiff = NULL # comparison value to show in plot of scales difference
  compValNu = NULL # comparison value to show in plots of normality
  compValEff = NULL # comparison value to show in plot of effect size
  ropeEffSz = NULL # ROPE around null value of effect size
  ropeMu = NULL
  ropeSigma = NULL
  
  diagnostics = TRUE # generate diagnostics plots or not
  computeEffsz = TRUE # calculate and plot effect size or not
  groupNames = c("FBL-QP","QP") # names of the groups
  subscript = "" # subscript to identify measured variable
  subsEffsz = "" # subscript to identify the effect size expression
  
  # Getting number of groups from data:
  g = as.numeric(as.factor(datFrm[,gName]))
  nG = max(g) # number of groups
  
  # Setting parameters for prior distributions:
  y = (as.numeric(datFrm[,yName]))
  muPrior = c( mean(y[ g==1 ]), mean(y[ g==2 ]) ) # mean of means priors (normal distribution)
  muSdPrior = 100*c( sd(y[ g==1 ]), sd(y[ g==2 ]) ) # s.d. of means priors (normal distribution)
  sigmaPriorLow = c( sd(y[ g==1 ]), sd(y[ g==2 ]) )/1000 # low limits of scales priors (uniform distribution)
  sigmaPriorHigh = 1000*c( sd(y[ g==1 ]), sd(y[ g==2 ]) ) # high limits of scales priors (uniform distribution)
  
  print(c( mean(y[ g==1 ]), mean(y[ g==2 ]) ))
  print(c( sd(y[ g==1 ]), sd(y[ g==2 ]) ))
  
  
  y = log(as.numeric(datFrm[,yName]))
  
  ropeMu = c(-0.05*mean(log(y[ g==2 ])), 0.05*mean(log(y[ g==2 ])))  
  
  # Generating MCMC chain:
  codaSamples = genMCMC_lnorm(datFrm, shapePriorMean=muPrior, shapePriorSD=muSdPrior, 
                                scalePriorLow=sigmaPriorLow, 
                                scalePriorHigh=sigmaPriorHigh,
                                saveName=fileNameRoot, yName=yName, gName = gName)
  
  # Summary and diagnostics of the MCMC chain:
  smryMCMC_lnorm(codaSamples, nG=nG, saveName=fileNameRoot, 
                   diagnostics=diagnostics, graphFileType=graphFileType)
  
  # Plotting the posterior distributions:
  plotMCMC_lnorm(codaSamples, datFrm, compValMu=compValMu, 
                   compValMuDiff=compValMuDiff, compValSigma=compValSigma,
                   compValSigmaDiff=compValSigmaDiff, compValNu=compValNu, 
                   compValEff=compValEff, ropeEffSz=ropeEffSz,ropeSigma = ropeSigma, ropeMu=ropeMu,
                   graphFileType=graphFileType, saveName=fileNameRoot, 
                   groupNames=groupNames, yName=yName, gName = gName)
  
  
  
}