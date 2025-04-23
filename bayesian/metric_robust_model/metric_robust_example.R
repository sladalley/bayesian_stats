graphics.off()
rm(list=ls())

source("metric_robust_functions.R")

# Generating a data set using a uniform distribution just to test the script:
datFrm = data.frame (resp=integer(), cond=integer(), subID=integer())
groups = 2 # number of groups
subjects = 30 # number of subjects
minValue = 40 # minimum value for the uniform distribution
maxValue = 120 # maximum value for the uniform distribution
set.seed(7)
for (g in 1:groups) {
  for (s in 1:subjects) {
    y = runif(1, min=minValue, max=maxValue) # sampling from uniform distribution
    datFrm[nrow(datFrm) + 1,] <- c(y, g, s)
  }
}
gName = "cond"
yName = "resp"
sName = "subID"

# Defining directory for results and prefix of the files names:
fileNameRoot = "test"
dir.create(fileNameRoot)
fileNameRoot = paste0(fileNameRoot, "/", fileNameRoot)
graphFileType = "pdf"

# Choosing type of analysis ("paired" or "independent"):
analysis = "independent"


# ======= Paired/one-group analysis ===========================================
if (analysis == "paired") {
  # Generating data frame with differences between groups:
  diff = data.frame(resp=integer(), cond=integer(), subID=integer())
  for (subject in 1:subjects) {
    y1 = datFrm[datFrm$subID == subject & datFrm$cond == 1, yName]
    y2 = datFrm[datFrm$subID == subject & datFrm$cond == 2, yName]
    diff[nrow(diff) + 1,] <- c(y2 - y1, 1, subject)
  }
  datFrm = diff
  
  # Setting the comparison values and ROPE intervals:
  compValMu = 0.0 # comparison value to show in plots of means
  compValSigma = NULL # comparison value to show in plots of scales
  compValNu = NULL # comparison value to show in plots of normality
  compValEff = NULL # comparison value to show in plot of effect size
  ropeEffSz = c(-0.1, 0.1) # ROPE around null value of effect size
  nullValEff = 0 # null value to calculate effect size
  
  diagnostics = TRUE # generate diagnostics plots or not
  computeEffsz = TRUE # calculate and plot effect size or not
  subscript = "diff" # subscript to identify measured variable
  subsEffsz = "1-group" # subscript to identify the effect size expression
 
  # Getting number of groups from data:
  g = as.numeric(as.factor(datFrm[,gName]))
  nG = max(g) # number of groups
  
  # Setting parameters for prior distributions:
  y = as.numeric(datFrm[,yName])
  muPrior = mean(y) # mean of mean prior (normal distribution)
  muSdPrior = 100*sd(y) # s.d. of mean prior (normal distribution)
  sigmaPriorLow = sd(y)/1000 # low limit of scale prior (uniform distribution)
  sigmaPriorHigh = 1000*sd(y) # high limit of scale prior (uniform distribution)
  
  # Generating MCMC chain:
  codaSamples = genMCMC(datFrm, muPrior=muPrior, muSdPrior=muSdPrior, 
                        sigmaPriorLow=sigmaPriorLow, 
                        sigmaPriorHigh=sigmaPriorHigh,
                        saveName=fileNameRoot, yName=yName)
  
  # Summary and diagnostics of the MCMC chain:
  smryMCMC(codaSamples, nG=nG, nullValEff=nullValEff, saveName=fileNameRoot, 
           diagnostics=diagnostics, graphFileType=graphFileType, 
           computeEffsz=computeEffsz)
  
  # Plotting the posterior distributions:
  plotMCMC(codaSamples, datFrm, nullValEff=nullValEff, compValMu=compValMu,
           compValSigma=compValSigma,
           compValNu=compValNu, 
           compValEff=compValEff, ropeEffSz=ropeEffSz, 
           graphFileType=graphFileType, saveName=fileNameRoot, 
           yName=yName, 
           subscript=subscript, subsEffsz=subsEffsz, plotEffsz=computeEffsz)
}

# ======= Independent/two-groups analysis =====================================
if (analysis == "independent") {
  # Setting the comparison values and ROPE intervals:
  compValMu = NULL # comparison value to show in plots of means
  compValMuDiff = 0.0 # comparison value to show in plot of means difference
  compValSigma = NULL # comparison value to show in plots of scales
  compValSigmaDiff = NULL # comparison value to show in plot of scales difference
  compValNu = NULL # comparison value to show in plots of normality
  compValEff = NULL # comparison value to show in plot of effect size
  ropeEffSz = c(-0.1, 0.1) # ROPE around null value of effect size
  
  diagnostics = TRUE # generate diagnostics plots or not
  computeEffsz = TRUE # calculate and plot effect size or not
  groupNames = c("1","2") # names of the groups
  subscript = "" # subscript to identify measured variable
  subsEffsz = "" # subscript to identify the effect size expression
  
  # Getting number of groups from data:
  g = as.numeric(as.factor(datFrm[,gName]))
  nG = max(g) # number of groups
  
  # Setting parameters for prior distributions:
  y = as.numeric(datFrm[,yName])
  muPrior = c( mean(y[ g==1 ]), mean(y[ g==2 ]) ) # mean of means priors (normal distribution)
  muSdPrior = 100*c( sd(y[ g==1 ]), sd(y[ g==2 ]) ) # s.d. of means priors (normal distribution)
  sigmaPriorLow = c( sd(y[ g==1 ]), sd(y[ g==2 ]) )/1000 # low limits of scales priors (uniform distribution)
  sigmaPriorHigh = 1000*c( sd(y[ g==1 ]), sd(y[ g==2 ]) ) # high limits of scales priors (uniform distribution)
  
  # Generating MCMC chain:
  codaSamples = genMCMC(datFrm, muPrior=muPrior, muSdPrior=muSdPrior, 
                        sigmaPriorLow=sigmaPriorLow, 
                        sigmaPriorHigh=sigmaPriorHigh,
                        saveName=fileNameRoot, yName=yName)
  
  # Summary and diagnostics of the MCMC chain:
  smryMCMC(codaSamples, nG=nG, saveName=fileNameRoot, 
           diagnostics=diagnostics, graphFileType=graphFileType, 
           computeEffsz=computeEffsz)
  
  # Plotting the posterior distributions:
  plotMCMC(codaSamples, datFrm, compValMu=compValMu, 
           compValMuDiff=compValMuDiff, compValSigma=compValSigma,
           compValSigmaDiff=compValSigmaDiff, compValNu=compValNu, 
           compValEff=compValEff, ropeEffSz=ropeEffSz, 
           graphFileType=graphFileType, saveName=fileNameRoot, 
           groupNames=groupNames, yName=yName, 
           subscript=subscript, subsEffsz=subsEffsz, plotEffsz=computeEffsz)
}

# =============================================================================

graphics.off()