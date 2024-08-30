graphics.off()
rm(list=ls())

source("metric_robust_functions.R")

# Loading data:
load("data_metric_robust_example.Rdata")
analysis_options = c("paired", "independent") # choose one option
analysis = analysis_options[2]

# Prefix for results files names and file type:
fileNameRoot = "results"
dir.create(fileNameRoot)
fileNameRoot = paste0(fileNameRoot, "/", fileNameRoot)
graphFileType = "pdf"

# Some parameters for the plots:
groupNames = c("1","2") # names of the groups for independent analysis
subscript = "" # subscript to identify variable (ex: "time", "x", etc)
subsEffsz = "" # subscript to identify the effect size expression (ex: "one-group", "two-groups")

# Paired analysis (one group):
if (analysis == "paired") {
  datFrm = datFrmDiff
  
  # Setting the null values and the ROPE intervals for each parameter:
  nullHypValMu = 0.0 # null hypothesis value to calculate effect size
  compValMu = 0.0 # real value to show in posterior plot
  compValSigma = NULL # real value to show in posterior plot
  compValNu = NULL # comparative value to show in posterior plot
  compValEff = NULL # real value to show in posterior plot
  ropeEffSz = c(-0.1, 0.1)
  # Extra parameters:
  compValMuDiff = NULL
  compValSigmaDiff = NULL
}

# Independent analysis (two groups):
if (analysis == "independent") {
  compValMu = NULL # real values to show in posterior plots
  compValMuDiff = 0.0 # real value to show in posterior plots
  compValSigma = NULL # real values to show in posterior plots
  compValSigmaDiff = NULL # real value to show in posterior plots
  compValNu = NULL # comparative values to show in posterior plots
  compValEff = NULL # real value to show in posterior plots
  ropeEffSz = c(-0.1, 0.1)
  # Extra parameters:
  nullHypValMu = NULL 
}

# Getting number of groups from data:
g = as.numeric(as.factor(datFrm[,gName]))
nG = max(g) # number of groups

# Setting priors:
y = as.numeric(datFrm[,yName])
if (nG == 1) {
  muPrior = mean(y)
  muSdPrior = 100*sd(y)
  sigmaPriorLow = sd(y)/1000
  sigmaPriorHigh = 1000*sd(y)
}
if (nG == 2) {
  muPrior = c( mean(y[ g==1 ]), mean(y[ g==2 ]) )
  muSdPrior = 100*c( sd(y[ g==1 ]), sd(y[ g==2 ]) )
  sigmaPriorLow = c( sd(y[ g==1 ]), sd(y[ g==2 ]) )/1000
  sigmaPriorHigh = 1000*c( sd(y[ g==1 ]), sd(y[ g==2 ]) )
}

# Generating MCMC chain:
codaSamples = genMCMC(datFrm, muPrior=muPrior, muSdPrior=muSdPrior, 
                      sigmaPriorLow=sigmaPriorLow, sigmaPriorHigh=sigmaPriorHigh,
                      saveName=fileNameRoot, yName="resp")

# Summary and diagnostics of the MCMC chain:
smryMCMC(codaSamples, nG=nG, nullHypValMu=nullHypValMu, saveName=fileNameRoot, 
         diagnostics=FALSE, graphFileType=graphFileType)

# Plotting the posterior distributions:
plotMCMC(codaSamples, datFrm, nullHypValMu=nullHypValMu, compValMu=compValMu, 
         compValMuDiff=compValMuDiff, compValSigma=compValSigma,
         compValSigmaDiff=compValSigmaDiff, compValNu=compValNu, 
         compValEff=compValEff, ropeEffSz=ropeEffSz, graphFileType=graphFileType,
         saveName=fileNameRoot, groupNames=groupNames, yName="resp", 
         subscript=subscript, subsEffsz=subsEffsz)

graphics.off()