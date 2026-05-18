graphics.off()
rm(list=ls())

source("metric_robust_functions.R")
source("weibull_functions.R")



# Generating a data set using a uniform distribution just to test the script:

datFrm = read.csv(file="iros_data_6.csv")

datFrm = subset(datFrm, 
                Solved == 1 & 
                  IsStable == 1 & 
                  ConstraintsOK == 1 & JointsOK == 1 & SimTime <= 100)

gName = "Controller"
yName = "SimTime"

# Defining directory for results and prefix of the files names:
fileNameRoot = "iros_test_6_weib_rmse_test_s"
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
  y = as.numeric(datFrm[,yName])
  muPrior = c( mean(y[ g==1 ]), mean(y[ g==2 ]) ) # mean of means priors (normal distribution)
  muSdPrior = c( sd(y[ g==1 ]), sd(y[ g==2 ]) ) # s.d. of means priors (normal distribution)
  sigmaPriorLow = c( 0, 0 ) # low limits of scales priors (uniform distribution)
  sigmaPriorHigh = 10*c( sd(y[ g==1 ]), sd(y[ g==2 ]) ) # high limits of scales priors (uniform distribution)
  
  print(c( mean(y[ g==1 ]), mean(y[ g==2 ]) ))
  print(c( sd(y[ g==1 ]), sd(y[ g==2 ]) ))
  
  ropeMu = c(-0.05*mean(y[ g==2 ]), 0.05*mean(y[ g==2 ]))  
  # Generating MCMC chain:
  codaSamples = genMCMC_weibull(datFrm, shapePriorMean=muPrior, shapePriorSD=muSdPrior, 
                        scalePriorLow=sigmaPriorLow, 
                        scalePriorHigh=sigmaPriorHigh,
                        saveName=fileNameRoot, yName=yName, gName = gName)
  
  # Summary and diagnostics of the MCMC chain:
  smryMCMC_weibull(codaSamples, nG=nG, saveName=fileNameRoot, 
           diagnostics=diagnostics, graphFileType=graphFileType)
  
  # Plotting the posterior distributions:
  plotMCMC_weibull(codaSamples, datFrm, compValMu=compValMu, 
           compValMuDiff=compValMuDiff, compValSigma=compValSigma,
           compValSigmaDiff=compValSigmaDiff, compValNu=compValNu, 
           compValEff=compValEff, ropeEffSz=ropeEffSz,ropeSigma = ropeSigma, ropeMu=ropeMu,
           graphFileType=graphFileType, saveName=fileNameRoot, 
           groupNames=groupNames, yName=yName, gName = gName)
  
  mcmcMat = as.matrix(codaSamples)
  powerEstimation_weibull(mcmcChain = mcmcMat, muNULL = meanPrior, effROPE = ropeEffSz, effHDImaxWid = 0.1, N = 30 ,nRep = 1, yName = yName, groupNames=groupNames, muROPE = ropeMu)
  
  

}