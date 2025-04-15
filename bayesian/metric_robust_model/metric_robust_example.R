graphics.off()
rm(list=ls())

source("metric_robust_functions.R")

# Generating a data set using a uniform distribution just to test the script:
datFrm = data.frame (resp=integer(), cond=integer(), subID=integer())
nG = 2 # number of groups
nS = 30 # number of subjects
minValue = 40
maxValue = 120
for (group in 1:nG) {
  for (subject in 1:nS) {
    y = runif(1, min=minValue, max=maxValue) # sampling from uniform distribution
    datFrm[nrow(datFrm) + 1,] <- c(y, group, subject)
  }
}
gName = "cond"
yName = "resp"
sName = "subID"

# Choosing type of analysis ("paired" or "independent"):
analysis = "independent"

# Defining directory for results and prefix of the files names:
fileNameRoot = "test"
dir.create(fileNameRoot)
fileNameRoot = paste0(fileNameRoot, "/", fileNameRoot)
graphFileType = "pdf"

diagnostics = TRUE
computeEffsz = TRUE

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
  compValMu = NULL # comparison value to show in plots of means
  compValMuDiff = 0.0 # comparison value to show in plot of means difference
  compValSigma = NULL # comparison value to show in plots of scales
  compValSigmaDiff = NULL # comparison value to show in plot of scales difference
  compValNu = NULL # comparison value to show in plots of normality
  compValEff = NULL # comparison value to show in plot of effect size
  ropeEffSz = c(-0.1, 0.1) # ROPE around the null value of effect size
  nullValEff = 0 # null value to calculate effect size
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
                      saveName=fileNameRoot, yName=yName)

# Summary and diagnostics of the MCMC chain:
smryMCMC(codaSamples, nG=nG, nullValEff=nullValEff, saveName=fileNameRoot, 
         diagnostics=diagnostics, graphFileType=graphFileType, 
         computeEffsz=computeEffsz)

# Plotting the posterior distributions:
plotMCMC(codaSamples, datFrm, nullValEff=nullValEff, compValMu=compValMu, 
         compValMuDiff=compValMuDiff, compValSigma=compValSigma,
         compValSigmaDiff=compValSigmaDiff, compValNu=compValNu, 
         compValEff=compValEff, ropeEffSz=ropeEffSz, graphFileType=graphFileType,
         saveName=fileNameRoot, groupNames=groupNames, yName=yName, 
         subscript=subscript, subsEffsz=subsEffsz, plotEffsz=computeEffsz)

graphics.off()