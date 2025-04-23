graphics.off()
rm(list=ls())

source("ordinal_robust_functions.R")

# Generating a data set using a uniform distribution just to test the script:
datFrm = data.frame (resp=integer(), item=integer(), 
                     cond=integer(), subID=integer())
groups = 2 # number of groups
subjects = 30 # number of subjects
items = 4 # number of items
minLevel = 1 # minimum level of response scale
maxLevel = 5 # maximum level of response scale
set.seed(7)
for (g in 1:groups) {
  for (s in 1:subjects) {
    for (i in 1:items) {
      y = round(runif(1, min=minLevel, max=maxLevel)) # sampling from uniform distribution
      datFrm[nrow(datFrm) + 1,] <- c(y, i, g, s)
    }
  }
}
gName = "cond"
yName = "resp"
sName = "subID"
qName = "item"

# Defining directory for results and prefix of the files names:
fileNameRoot = "test"
dir.create(fileNameRoot)
fileNameRoot = paste0(fileNameRoot, "/", fileNameRoot)
graphFileType = "pdf"

# Choosing type of analysis ("paired" or "independent"):
analysis = "independent"

# Choosing to treat the data with empty levels or not
# If TRUE, a single answer will be added to each empty level found
treatEmpty = TRUE

# ======= Paired/one-group analysis =====================================
if (analysis == "paired") {
  # The difference will be calculated for each item and the final value will
  # be added to a constant (maxLevel) to make all the possible values
  # positive. Example: if minLevel = 1 and maxLevel = 5, the minimum value
  # of the difference is -4 (1-5) and the maximum is 4 (5-1). So we will add 
  # this to 5 so the new possible values are from 1 to 9.
  diff = data.frame (resp=integer(), item=integer(), 
                     cond=integer(), subID=integer())
  for (subject in 1:subjects) {
    for (item in 1:items) {
      y1 = datFrm[datFrm$subID == subject & datFrm$item == item & datFrm$cond == 1, yName]
      y2 = datFrm[datFrm$subID == subject & datFrm$item == item & datFrm$cond == 2, yName]
      diff[nrow(diff) + 1,] <- c(y2 - y1 + maxLevel, item, 1, subject)
    }
  }
  datFrm = diff
  
  minLevel = 1 # minimum level level
  maxLevel = 9 # maximum level level
  compValMu = 5.0 # comparison value to show in plots of means
  compValSigma = NULL # comparison value to show in plots of scales
  compValNu = NULL # comparison value to show in plots of normality
  compValEff = NULL # comparison value to show in plot of effect size
  ropeEffSz = c(-0.1, 0.1) # ROPE around null value of effect size
  nullValEff = 5.0 # # null value to calculate effect size
  
  diagnostics = TRUE # generate diagnostics plots or not
  computeEffsz = TRUE # calculate and plot effect size or not
  subscript = "diff" # subscript to identify measured variable
  subsEffsz = "1-group" # subscript to identify the effect size expression
  
  # Getting number of items and groups from data:
  q = as.numeric(as.factor(datFrm[,qName]))
  g = as.numeric(as.factor(datFrm[,gName]))
  nQ = max(q) # number of items
  nG = max(g) # number of groups
  
  if (treatEmpty) {
    # Plotting data histograms before adding extra data:
    plotDataHistograms(datFrm, Nlevels=maxLevel, graphFileType=graphFileType,
                       saveName=fileNameRoot)
    
    # Adding data to remove empty levels:
    extraList <- onlyEmptyExtraData(datFrm, Nlevels=maxLevel)
    datFrm = extraList$data
    extraInfo = extraList$extra
  } else { 
    extraInfo = NULL
  }
  
  # Setting parameters for prior distributions:
  muPrior = rep(0.0,nG)+((minLevel+maxLevel)/2) # mean of means priors (normal distribution)
  muSdPrior = rep(0.0,nG)+maxLevel # s.d. of means priors (normal distribution)
  sigmaPriorLow = rep(0.0,nG)+0.01 # low limits of scales priors (uniform distribution)
  sigmaPriorHigh = rep(0.0,nG)+(10*maxLevel) # high limits of scales priors (uniform distribution)
  threshPrior = matrix( rep( seq( minLevel+0.5 , maxLevel-0.5 , 1 ) , nQ ),
                        nrow=nQ , byrow=TRUE ) # mean of thresholds priors (normal distribution)
  threshSdPrior = 2.0 # s.d. of thresholds priors (normal distribution)
  
  # Generating MCMC chain:
  ordCodaSamples = genMCMC(datFrm, maxLevel=maxLevel, muPrior=muPrior, 
                           muSdPrior=muSdPrior, sigmaPriorLow=sigmaPriorLow, 
                           sigmaPriorHigh=sigmaPriorHigh, 
                           threshPrior=threshPrior, threshSdPrior=threshSdPrior, 
                           saveName=fileNameRoot)
  
  # Summary and diagnostics of the MCMC chain:
  smryMCMC(ordCodaSamples, nG=nG, nullValEff=nullValEff, 
           saveName=fileNameRoot, diagnostics=diagnostics, 
           graphFileType=graphFileType)
  
  # Plotting the posterior distributions:
  plotMCMC(ordCodaSamples, datFrm, nullValEff=nullValEff, 
           compValMu=compValMu,
           compValSigma=compValSigma,
           compValNu=compValNu, compValEff=compValEff, ropeEffSz=ropeEffSz, 
           minLevel=minLevel, maxLevel=maxLevel, graphFileType=graphFileType,
           saveName=fileNameRoot,
           subscript = subscript, subsEffsz = subsEffsz, 
           extraInfo = extraInfo)
}

# ======= Independent/two-groups analysis =====================================
if (analysis == "independent") {
  # Setting the comparison values and ROPE intervals:
  minLevel = 1 # minimum response level
  maxLevel = 5 # maximum response level
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
  
  # Getting number of items and groups from data:
  q = as.numeric(as.factor(datFrm[,qName]))
  g = as.numeric(as.factor(datFrm[,gName]))
  nQ = max(q) # number of items
  nG = max(g) # number of groups
  
  if (treatEmpty) {
    # Plotting data histograms before adding extra data:
    plotDataHistograms(datFrm, Nlevels=maxLevel, graphFileType=graphFileType,
                       saveName=fileNameRoot, groupNames=groupNames)
    
    # Adding data to remove empty levels:
    extraList1 = onlyEmptyExtraData(datFrm[g==1,], Nlevels=maxLevel)
    extraList2 = onlyEmptyExtraData(datFrm[g==2,], Nlevels=maxLevel)
    
    datFrm = rbind(extraList1$data, extraList2$data)
    extraInfo = list(extraList1$extra, extraList2$extra)
  } else { 
    extraInfo = NULL
  }
  
  # Setting parameters for prior distributions:
  muPrior = rep(0.0,nG)+((minLevel+maxLevel)/2) # mean of means priors (normal distribution)
  muSdPrior = rep(0.0,nG)+maxLevel # s.d. of means priors (normal distribution)
  sigmaPriorLow = rep(0.0,nG)+0.01 # low limits of scales priors (uniform distribution)
  sigmaPriorHigh = rep(0.0,nG)+(10*maxLevel) # high limits of scales priors (uniform distribution)
  threshPrior = matrix( rep( seq( minLevel+0.5 , maxLevel-0.5 , 1 ) , nQ ),
                        nrow=nQ , byrow=TRUE ) # mean of thresholds priors (normal distribution)
  threshSdPrior = 2.0 # s.d. of thresholds priors (normal distribution)
  
  # Generating MCMC chain:
  ordCodaSamples = genMCMC(datFrm, maxLevel=maxLevel, muPrior=muPrior, 
                           muSdPrior=muSdPrior, sigmaPriorLow=sigmaPriorLow, 
                           sigmaPriorHigh=sigmaPriorHigh, 
                           threshPrior=threshPrior, threshSdPrior=threshSdPrior, 
                           saveName=fileNameRoot)
  
  # Summary and diagnostics of the MCMC chain:
  smryMCMC(ordCodaSamples, nG=nG, 
           saveName=fileNameRoot, diagnostics=diagnostics, 
           graphFileType=graphFileType)
  
  # Plotting the posterior distributions:
  plotMCMC(ordCodaSamples, datFrm, 
           compValMu=compValMu, compValMuDiff=compValMuDiff, 
           compValSigma=compValSigma, compValSigmaDiff=compValSigmaDiff, 
           compValNu=compValNu, compValEff=compValEff, ropeEffSz=ropeEffSz, 
           minLevel=minLevel, maxLevel=maxLevel, graphFileType=graphFileType,
           saveName=fileNameRoot, groupNames=groupNames,
           extraInfo = extraInfo)
  
}

graphics.off()