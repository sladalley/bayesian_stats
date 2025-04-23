# Functions for the Bayesian analysis of ordinal data with robust 
# ordered-probit model.
# Based on Chapter 23 of Kruschke (2015) and Lidell&Kruschke (2018) and
# the scripts accompanying them.
#
# * Kruschke, J. K. (2015). Doing Bayesian Data Analysis, Second Edition: 
#   A Tutorial with R, JAGS, and Stan. Academic Press / Elsevier.
# * Liddell, T. M., Kruschke, J. K. (2018). Analyzing ordinal data with 
#   metric models: What could possibly go wrong? Journal of Experimental 
#   Social Psychology, Volume 79, Pages 328-348, ISSN 0022-1031,
#   https://doi.org/10.1016/j.jesp.2018.08.009.

source("DBDA2E-utilities.R")
library('latex2exp')

# ====== Extra data function ==================================================
onlyEmptyExtraData = function(data, yName="resp", qName="item", sName="subID", 
                              gName="cond", Nlevels) {
  # Add one extra answer for each empty level of each item.
  # List of parameters:
  #   - data: data set of a single group.
  #   - yName: name of the column with the answers (ordinal values).
  #   - qName: name of the column with the item indexes.
  #   - sName: name of the column with the subject indexes.
  #   - gName: name of the column with the group index.
  #   - Nlevels: number of possible ordinal levels.
  
  q = as.numeric(as.factor(data[,qName])) # items
  s = as.numeric(data[,sName]) # subjects
  g = as.numeric(data[,gName]) # group
  nQ = max(q) # number of items
  lastS = max(s) # last subject's id
  nS = length(s)/nQ # total number of subjects
  gIdx = g[1] # group index
  
  # One extra answer for each empty level of each item.
  extra = matrix(rep(0, Nlevels*nQ), nrow=nQ, byrow=TRUE)
  for (item in seq(1, nQ)) {
    for (level in seq(1,Nlevels)) {
      if (length( which(data[q==item,yName]==level) ) == 0)
        extra[item, level] = 1
    }
  }
  extra_return = extra
  
  # Add new answers in the data frame.
  for (item in seq(1, nQ)) {
    for (subj in seq(lastS+1, lastS+sum(extra[item,]))) {
      for (level in seq(1, Nlevels)) {
        if (extra[item, level] > 0) {
          data = rbind(data, c(0, 0, 0, 0))
          c(level, item, subj, gIdx)
          data[nrow(data), yName] = level
          data[nrow(data), qName] = item
          data[nrow(data), gName] = gIdx
          data[nrow(data), sName] = subj
          extra[item, level] = extra[item, level] - 1
          break
        }
      }
    }
  }
  
  return(list("data"=data, "extra"=extra_return))
}

# ====== Data histograms =======================================================
plotDataHistograms = function(data, yName="resp", qName="item", sName="subID", 
                              gName="cond", Nlevels, graphFileType="png", 
                              saveName=NULL, groupNames=c(1,2)) {
  # Plot histograms of ordinal responses.
  # List of parameters:
  #   - data: data set of a single group.
  #   - yName: name of the column with the answers (ordinal values).
  #   - qName: name of the column with the item indexes.
  #   - sName: name of the column with the subject indexes.
  #   - gName: name of the column with the group index.
  #   - Nlevels: number of possible ordinal levels.
  #   - graphFileType: type of the image output files.
  #   - saveName: prefix of the output files.
  #   - groupNames: names of each group to be displayed in the plots.
  y = as.numeric(datFrm[,yName])
  Ntotal = length(y)
  
  # Question and Group data vectors:
  q = as.numeric(as.factor(datFrm[,qName]))
  g = as.numeric(as.factor(datFrm[,gName]))
  qLevels = levels(as.factor(datFrm[,qName]))
  gLevels = levels(as.factor(datFrm[,gName]))
  nQ = max(q)
  nG = max(g)
  
  dataScaleDensMax = 0.7
  openGraph(height=min(2.5*nQ+0.75,14),width=3.5*nG)
  par( mar=c(3.5,3.5,2.5,0.5) , mgp=c(2.0,0.7,0) , oma=c(0,0,3.5,0) ) # , xpd=NA )
  layout( matrix(1:(nG*nQ),nrow=nQ,ncol=nG,byrow=FALSE) )
  for ( gIdx in 1:nG ) {
    for ( qIdx in 1:nQ ) {
      probInfo = NULL
      # Data histogram:
      thisY = y[ g==gIdx & q==qIdx ]
      xLim = c( minLevel-0.5 , maxLevel+0.5 )
      xBreaks = seq( xLim[1] , xLim[2] , 1 )  
      if (nG == 1) {
        if (qIdx < nQ) {
          # Only y label
          histInfo = hist( thisY , prob=TRUE , ylim=c(0,dataScaleDensMax) , 
                           xlim=xLim , breaks=xBreaks ,
                           xlab="" , ylab="p(response)" , cex.axis=1.2 , 
                           cex.lab=1.5 ,
                           col="#f4e3d7ff" , border="white" , xaxt="n" ,
                           main=paste0("Item ",qIdx), cex.main=1.5)
        } else {
          # Both labels
          histInfo = hist( thisY , prob=TRUE , ylim=c(0,dataScaleDensMax) , 
                           xlim=xLim , breaks=xBreaks ,
                           xlab="Ordinal responses" , ylab="p(response)" , 
                           cex.axis=1.2, cex.lab=1.5 , # yaxt="n" ,
                           col="#f4e3d7ff" , border="white" , xaxt="n" ,
                           main=paste0("Item ",qIdx), cex.main=1.5)
        }
      } else {
        if (gIdx > 1 & qIdx == nQ) {
          # Only x label
          histInfo = hist( thisY , prob=TRUE , ylim=c(0,dataScaleDensMax) , 
                           xlim=xLim , breaks=xBreaks ,
                           xlab="Ordinal responses" , ylab="" , cex.axis=1.2,
                           cex.lab=1.5 , yaxt="n" , 
                           col="#f4e3d7ff" , border="white" , xaxt="n" ,
                           main=paste0("Item ",qIdx,", Group ",groupNames[gIdx]),
                           cex.main=1.5)
        } else if (gIdx == 1 & qIdx < nQ) {
          # Only y label
          histInfo = hist( thisY , prob=TRUE , ylim=c(0,dataScaleDensMax) , 
                           xlim=xLim , breaks=xBreaks ,
                           xlab="" , ylab="p(response)" , cex.axis=1.2,
                           cex.lab=1.5 , # yaxt="n" ,
                           col="#f4e3d7ff" , border="white" , xaxt="n" ,
                           main=paste0("Item ",qIdx,", Group ",groupNames[gIdx]),
                           cex.main=1.5)
        } else if (gIdx > 1 & qIdx < nQ) {
          # No labels
          histInfo = hist( thisY , prob=TRUE , ylim=c(0,dataScaleDensMax) , 
                           xlim=xLim , breaks=xBreaks ,
                           xlab="" , ylab="" , cex.axis=1.2, 
                           cex.lab=1.5 , yaxt="n" ,
                           col="#f4e3d7ff" , border="white" , xaxt="n" ,
                           main=paste0("Item ",qIdx,", Group ",groupNames[gIdx]),
                           cex.main=1.5)
        } else {
          # Both labels
          histInfo = hist( thisY , prob=TRUE , ylim=c(0,dataScaleDensMax) , 
                           xlim=xLim , breaks=xBreaks ,
                           xlab="Ordinal responses" , ylab="p(response)" , 
                           cex.axis=1.2, cex.lab=1.5 , #yaxt="n" ,
                           col="#f4e3d7ff" , border="white" , xaxt="n" ,
                           main=paste0("Item ",qIdx,", Group ",groupNames[gIdx]),
                           cex.main=1.5)
        }
      }
      axis(side=1, seq(minLevel, maxLevel, 1), seq(minLevel, maxLevel, 1))
    }
  }
  saveGraph( paste0(saveName,"-DataHistogram") , type=graphFileType )
}

# ====== MCMC chain ===========================================================
genMCMC = function (datFrm, maxLevel, yName="resp", qName="item", gName="cond", 
                    muPrior, muSdPrior, sigmaPriorLow, sigmaPriorHigh, 
                    threshPrior, threshSdPrior, saveName=NULL) {
  # Generate the MCMC chain using JAGS.
  # List of parameters:
  #   - datFrm: data set.
  #   - maxLevel: the maximum ordinal level.
  #   - yName: name of the column with the answers (ordinal values).
  #   - qName: name of the column with the item indexes.
  #   - gName: name of the column with the group indexes.
  #   - muPrior: mean of the prior for the mean of each group (normal dist.).
  #   - muSdPrior: std. dev. of the prior for the mean of each group (normal dist.).
  #   - sigmaPriorLow: low limit of the prior for the scale of each group (uniform dist.).
  #   - sigmaPriorHigh: high limit of the prior for the scale of each group (uniform dist.).
  #   - threshPrior: mean of the priors of the thresholds (normal dist.).
  #   - threshSdPrior: std. dev. of the priors of the thresholds (normal dist.).
  #   - saveName: prefix of the output files.
  # ------------------------------------------------------------------------- #
  # ASSEMBLE THE DATA FOR JAGS.
  
  # Rename and reclass y values for convenience:
  y = as.numeric(datFrm[,yName])
  # Do some checking that data make sense:
  if ( any( y!=round(y) ) ) { stop("All y values must be integers (whole numbers).") }
  if ( any( y < 1 ) ) { stop("All y values must be 1 or larger.") }
  
  Ntotal = length(y)
  # Question and Group data vectors:
  q = as.numeric(as.factor(datFrm[,qName]))
  g = as.numeric(as.factor(datFrm[,gName]))
  qLevels = levels(as.factor(datFrm[,qName]))
  gLevels = levels(as.factor(datFrm[,gName]))
  nQ = max(q)
  nG = max(g)
  
  # Define number of Y levels for each question:
  nYlevels = rep( maxLevel, nQ )
  
  # Create threshold matrix. For the first question, threshold 1 and nYlevels-1
  # are fixed; other interior thresholds are estimated. For other questions, 
  # all thresholds are estimated.
  thresh = matrix( NA , nrow=nQ , ncol=max(nYlevels)-1 ) # default to NA
  # Fix low thresh of item 1 at 1.5:
  thresh[1,1] = 1.5
  # Fix upper thresh of item 1 at K-0.5:
  thresh[1,nYlevels[1]-1] = maxLevel-0.5
  
  # Specify the data in a list, for ORDERED-PROBIT model, for later shipment to
  # JAGS:
  ordDataList = list(
    y = y , # answers
    q = q , # items
    g = g , # groups
    thresh = thresh , # thresholds
    nYlevels = nYlevels , # number of levels
    nQ = nQ , # number of items
    nG = nG , # number of groups
    Ntotal = Ntotal , # total size of data set
    muPrior = muPrior , # mean of the normal prior of the mean
    muSdPrior = muSdPrior , # std. dev. of the normal prior of the mean
    sigmaPriorLow = sigmaPriorLow , # low limit of the uniform prior of sigma
    sigmaPriorHigh = sigmaPriorHigh , # high limit of the uniform prior of sigma
    threshPrior = threshPrior , # mean of the normal priors of the thresholds
    threshSdPrior = threshSdPrior # std. dev. of the normal priors of the thresholds
  )
  
  # --------------------------------------------------------------------------- #
  # THE ORDERED-PROBIT MODEL FOR MULTIPLE ITEMS.
  # N.B.: THIS MODEL ASSUMES THAT ALL ITEMS ARE SCORED WITH
  # POSITIVE CORRELATIONS WITH EACH OTHER.
  
  ordModelString = "
  model {
    for ( i in 1:Ntotal ) {
      y[i] ~ dcat( pr[i,1:nYlevels[q[i]]] )
      pr[i,1] <- pt( thresh[q[i],1] ,
                     mu[g[i]] , 1/sigma[g[i]]^2 , nu[g[i]] )
      for ( k in 2:(nYlevels[q[i]]-1) ) {
        pr[i,k] <- max( 0 ,  pt( thresh[q[i], k]  ,
                                 mu[g[i]] , 1/sigma[g[i]]^2 , nu[g[i]] )
                           - pt( thresh[q[i],k-1] ,
                                 mu[g[i]] , 1/sigma[g[i]]^2 , nu[g[i]] ) )
      }
      pr[i,nYlevels[q[i]]] <- (
         1 - pt( thresh[q[i], nYlevels[q[i]]-1] ,
                 mu[g[i]] , 1/sigma[g[i]]^2, nu[g[i]] ) )
    }
    for ( gIdx in 1:nG ) { 
      mu[gIdx] ~ dnorm( muPrior[gIdx] , 1/(muSdPrior[gIdx])^2 )
      sigma[gIdx] ~ dunif( sigmaPriorLow[gIdx] , sigmaPriorHigh[gIdx] )
      nu[gIdx] ~ dexp(1/30.0)
    }
    # Prior on thresh[q,k]. Stochastic for all except thresh[1,1] and thresh[1,last].
    for ( qIdx in 1 ) {
      for ( kIdx in 2:(nYlevels[qIdx]-2) ) { # 1 and nYlevels-1 are fixed
        thresh[qIdx,kIdx] ~ dnorm( threshPrior[qIdx, kIdx] , 1/(threshSdPrior)^2 )
      }
    }
    for ( qIdx in 2:nQ ) {
      for ( kIdx in 1:(nYlevels[qIdx]-1) ) { 
        thresh[qIdx,kIdx] ~ dnorm( threshPrior[qIdx, kIdx] , 1/(threshSdPrior)^2 )
      }
    }
}" # close quote for ordModelString
  # Write out ordModelString to a text file
  writeLines( ordModelString , con=paste0(saveName, "-Model.txt") )
  
  # --------------------------------------------------------------------------- #
  # RUN THE CHAINS FOR ORDERED-PROBIT MODEL.
  
  parameters = c( "mu" , "sigma" , "thresh" , "nu")
  numSavedSteps = 20000 # 20000 
  thinSteps = 5 # 5 
  adaptSteps = 500               # Number of steps to "tune" the samplers
  burnInSteps = 1000
  #saveName=fileNameRoot 
  runjagsMethod=runjagsMethodDefault # from DBDA2E-utilities
  nChains=nChainsDefault # from DBDA2E-utilities
  
  ordRunJagsOut <- run.jags( method="parallel" , # runjagsMethod ,
                             model=paste0(saveName, "-Model.txt") , 
                             monitor=parameters , 
                             data=ordDataList ,  
                             #inits=initsList , 
                             n.chains=nChains ,
                             adapt=adaptSteps ,
                             burnin=burnInSteps , 
                             sample=ceiling(numSavedSteps/nChains) ,
                             thin=thinSteps ,
                             summarise=FALSE ,
                             plots=FALSE )
  ordCodaSamples = as.mcmc.list( ordRunJagsOut )
  # resulting codaSamples object has these indices: 
  #   codaSamples[[ chainIdx ]][ stepIdx , paramIdx ]
  if ( !is.null(saveName) ) {
    save( ordCodaSamples , file=paste(saveName,"-Ord-Mcmc.Rdata",sep="") )
  }
  
  return(ordCodaSamples)
}

# ====== Summary of data and diagnostics ======================================
smryMCMC = function (ordCodaSamples, nG, nullValEff=0, saveName=NULL, 
                     diagnostics=TRUE, graphFileType="png", computeEffsz=TRUE) {
  # Compute summary statistics of the chain and generate diagnostics plots.
  # List of parameters:
  #   - ordCodaSamples: codaSamples object with the MCMC chain.
  #   - nG: number of groups.
  #   - nullValEff: null value to calculate effect size when single group.
  #   - saveName: prefix of the output files.
  #   - diagnostics: to define if it should generate diagnostics plots.
  #   - graphFileType: type of the image output files.
  #   - computeEffsz: if effect size should be calculated (TRUE) or not (FALSE)
  # --------------------------------------------------------------------------- #
  # SUMMARIZING THE MCMC POSTERIOR.
  
  summaryInfo = NULL
  mcmcMat = as.matrix(ordCodaSamples, chains=TRUE)
  # One group:
  if (nG == 1) {
    # Mean:
    summaryInfo = rbind( summaryInfo , 
                         "mu" = summarizePost( mcmcMat[,"mu"] , 
                                               compVal=NULL , ROPE=NULL ) )
    # Scale parameter:
    summaryInfo = rbind( summaryInfo , 
                         "sigma" = summarizePost( mcmcMat[,"sigma"] , 
                                                  compVal=NULL , ROPE=NULL ) )
    # Normality parameter:
    summaryInfo = rbind( summaryInfo , 
                         "nu" = summarizePost( mcmcMat[,"nu"] , 
                                               compVal=NULL , ROPE=NULL ) )
    summaryInfo = rbind( summaryInfo , 
                         "log10(nu)" = summarizePost( log10(mcmcMat[,"nu"]) , 
                                                      compVal=NULL , ROPE=NULL ) )
    # Effect size, (mu - mu0) / sigma:
    if (computeEffsz) {
      summaryInfo = rbind( summaryInfo , 
                           "effSz" = summarizePost( 
                             ( mcmcMat[,"mu"] - nullValEff ) / mcmcMat[,"sigma"] ,
                             compVal=NULL , ROPE=NULL ) )
    }
  }
  
  # Two groups:
  if (nG == 2) {
    # Mean of group 1: 
    summaryInfo = rbind( summaryInfo , 
                         "mu[1]" = summarizePost( mcmcMat[,"mu[1]"] , 
                                                  compVal=NULL , ROPE=NULL ) )
    # Mean of group 2:
    summaryInfo = rbind( summaryInfo , 
                         "mu[2]" = summarizePost( mcmcMat[,"mu[2]"] , 
                                                  compVal=NULL , ROPE=NULL ) )
    # Difference between means, mu2 - mu1:
    summaryInfo = rbind( summaryInfo , 
                         "muDiff" = summarizePost( 
                           mcmcMat[,"mu[2]"] - mcmcMat[,"mu[1]"] , 
                           compVal=NULL , ROPE=NULL ) )
    # Scale parameter of group 1:
    summaryInfo = rbind( summaryInfo , 
                         "sigma[1]" = summarizePost( mcmcMat[,"sigma[1]"] , 
                                                     compVal=NULL , ROPE=NULL ) )
    # Scale parameter of group 2:
    summaryInfo = rbind( summaryInfo , 
                         "sigma[2]" = summarizePost( mcmcMat[,"sigma[2]"] , 
                                                     compVal=NULL , ROPE=NULL ) )
    # Difference between scales, sigma2 - sigma1:
    summaryInfo = rbind( summaryInfo , 
                         "sigmaDiff" = summarizePost( 
                           mcmcMat[,"sigma[2]"] - mcmcMat[,"sigma[1]"] , 
                           compVal=NULL , ROPE=NULL ) )
    # Normality parameter of group 1:
    summaryInfo = rbind( summaryInfo , 
                         "nu[1]" = summarizePost( mcmcMat[,"nu[1]"] , 
                                                  compVal=NULL , ROPE=NULL ) )
    summaryInfo = rbind( summaryInfo , 
                         "log10(nu[1])" = summarizePost( log10(mcmcMat[,"nu[1]"]) , 
                                                         compVal=NULL , ROPE=NULL ) )
    # Normality parameter of group 2:
    summaryInfo = rbind( summaryInfo , 
                         "nu[2]" = summarizePost( mcmcMat[,"nu[2]"] , 
                                                  compVal=NULL , ROPE=NULL ) )
    summaryInfo = rbind( summaryInfo , 
                         "log10(nu[2])" = summarizePost( log10(mcmcMat[,"nu[2]"]) , 
                                                         compVal=NULL , ROPE=NULL ) )
    # Effect size, (mu2 - mu1) / sqrt((sigma1^2 + sigma2^2)/2):
    if (computeEffsz) {
      summaryInfo = rbind( summaryInfo , 
                           "effSz" = summarizePost( 
                             ( mcmcMat[,"mu[2]"]-mcmcMat[,"mu[1]"] ) 
                             / sqrt((mcmcMat[,"sigma[1]"]^2+mcmcMat[,"sigma[2]"]^2)/2) ,
                             compVal=NULL , ROPE=NULL ) )
    }
  }
  # Thresholds (the same set for all groups)
  for ( colName in grep( "thresh" , colnames(mcmcMat) , value=TRUE ) ) {
    summaryInfo = rbind( summaryInfo , 
                         summarizePost( mcmcMat[,colName] , 
                                        compVal=NULL , ROPE=NULL ) )
    rownames(summaryInfo)[nrow(summaryInfo)] = colName
  }
  
  if ( !is.null(saveName) ) {
    write.csv( summaryInfo , file=paste0(saveName,"-SummaryInfo.csv") )
  }
  
  # --------------------------------------------------------------------------- #
  # DIAGNOSTICS
  
  if ( diagnostics == TRUE ) {
    diagnosticsPath = paste0(sub("/.*", "", saveName), "/diagnostics")
    dir.create(diagnosticsPath)
    ordParameterNames = varnames(ordCodaSamples) 
    for ( parName in ordParameterNames ) {
      diagMCMC( codaObject=ordCodaSamples , parName=parName ,  
                saveName=paste0(diagnosticsPath, "/", sub(".*/", "", saveName), "-") , 
                saveType=graphFileType )
    }
  }
}

# ====== Posterior distributions ==============================================
plotMCMC = function (ordCodaSamples, datFrm, yName="resp", qName="item", 
                     gName="cond", nullValEff=0, compValMu=NULL, 
                     compValMuDiff=NULL, compValSigma=NULL, compValSigmaDiff=NULL, 
                     compValNu=NULL, compValEff=NULL, ropeEffSz=NULL, 
                     minLevel, maxLevel, graphFileType="png", saveName=NULL, 
                     groupNames=c(1,2), subscript="", subsEffsz="", 
                     extraInfo=NULL, plotEffsz=TRUE) {
  # Display posterior information.
  # List of parameters:
  #   - ordCodaSamples: codaSamples object with the MCMC chain.
  #   - datFrm: data set.
  #   - yName: name of the column with the answers (ordinal values).
  #   - qName: name of the column with the item indexes.
  #   - gName: name of the column with the group indexes. 
  #   - nullValEff: null value to calculate effect size when single group.
  #   - compValMu: comparison value in the posterior mu plot(s).
  #   - compValMuDiff: comparison value in the posterior of mu difference plot (two groups).
  #   - compValSigma: comparison value in the posterior sigma plot(s).
  #   - compValSigmaDiff: comparison value in the posterior of sigma difference plot (two groups).
  #   - compValNu: comparison value in the posterior nu plot(s).
  #   - compValEff: comparison value in the posterior effect size plot.
  #   - ropeEffSz: ROPE in the posterior effect size plot.
  #   - minLevel: the minimum ordinal level.
  #   - maxLevel: the maximum ordinal level.
  #   - graphFileType: type of the image output files.
  #   - saveName: prefix of the output files.
  #   - groupNames: names of each group to be displayed in the plots.
  #   - subscript: subscript to indicate the variable.
  #   - subsEffsz: subscript to indicate the expression for the effect size.
  #   - extraInfo: information about extra data added.
  #   - plotEffSz: if effect size should be plotted (TRUE) or not (FALSE).
  ordMcmcMat = as.matrix( ordCodaSamples )
  ordChainLength = nrow(ordMcmcMat)
  
  y = as.numeric(datFrm[,yName])
  Ntotal = length(y)
  # Question and Group data vectors:
  q = as.numeric(as.factor(datFrm[,qName]))
  g = as.numeric(as.factor(datFrm[,gName]))
  qLevels = levels(as.factor(datFrm[,qName]))
  gLevels = levels(as.factor(datFrm[,gName]))
  nQ = max(q)
  nG = max(g)
  
  cex_plotPost = 1.0
  
  # Posterior mu and sigma of one or two groups:
  if (nG == 1) {
    # One group:
    openGraph(height=2.5*nG,width=7)
    par( mar=c(3.5,3.5,3,1) , mgp=c(2.25,0.7,0) )
    layout( matrix(1:(2*nG),nrow=nG,byrow=TRUE) )
    muLim = range( ordMcmcMat[, grep("^mu",colnames(ordMcmcMat)) ] )
    sigmaLim = range( ordMcmcMat[, grep("^sigma",colnames(ordMcmcMat)) ] )
    if (subscript != "") {
      plotPost( ordMcmcMat[,"mu"] , xlab=TeX(sprintf("$\\mu_{%s}$", subscript)) ,
                xlim=muLim , compVal=compValMu , cex=cex_plotPost , col="#0088aa8a")
      plotPost( ordMcmcMat[,"sigma"] , xlab=TeX(sprintf("$\\tau_{%s}$", subscript)) , 
                xlim=sigmaLim , compVal=compValSigma , cex=cex_plotPost , col="#0088aa8a")
    } else {
      plotPost( ordMcmcMat[,"mu"] , xlab=TeX(sprintf("$\\mu$")) ,
                xlim=muLim , compVal=compValMu , cex=cex_plotPost , col="#0088aa8a")
      plotPost( ordMcmcMat[,"sigma"] , xlab=TeX(sprintf("$\\tau$")) , 
                xlim=sigmaLim , compVal=compValSigma , cex=cex_plotPost , col="#0088aa8a")
    }
  }
  if (nG == 2) {
    # Two groups:
    openGraph(height=2.5*(nG+1),width=7)
    par( mar=c(3.5,3.5,3,1) , mgp=c(2.25,0.7,0) )
    layout( matrix(1:(2*(nG+1)),nrow=nG+1,byrow=TRUE) )
    muLim = range( ordMcmcMat[, grep("^mu\\[",colnames(ordMcmcMat)) ] )
    sigmaLim = range( ordMcmcMat[, grep("^sigma\\[",colnames(ordMcmcMat)) ] )
    if (subscript != "") {
      plotPost( ordMcmcMat[,"mu[1]"] , 
                xlab=TeX(sprintf("$\\mu_{%s,%s}$", subscript, groupNames[1])) , 
                xlim=muLim , compVal=compValMu[1] , 
                col="#0088aa8a" , cex=cex_plotPost )
      plotPost( ordMcmcMat[,"sigma[1]"] , 
                xlab=TeX(sprintf("$\\tau_{%s,%s}$", subscript, groupNames[1])) , 
                xlim=sigmaLim , compVal=compValSigma[1] , 
                col="#0088aa8a" , cex=cex_plotPost )
      plotPost( ordMcmcMat[,"mu[2]"] , 
                xlab=TeX(sprintf("$\\mu_{%s,%s}$", subscript, groupNames[2])) , 
                xlim=muLim , compVal=compValMu[2] , 
                col="#0088aa8a" , cex=cex_plotPost )
      plotPost( ordMcmcMat[,"sigma[2]"] , 
                xlab=TeX(sprintf("$\\tau_{%s,%s}$", subscript, groupNames[2])) , 
                xlim=sigmaLim , compVal=compValSigma[2] , 
                col="#0088aa8a" , cex=cex_plotPost )
      plotPost( ordMcmcMat[,"mu[2]"]-ordMcmcMat[,"mu[1]"] , compVal=compValMuDiff ,
                xlab=TeX(sprintf("$\\mu_{%s,%s} - \\mu_{%s,%s}$", subscript, groupNames[2], subscript, groupNames[1])) , 
                cex=cex_plotPost , col="#0088aa8a" )
      plotPost( ordMcmcMat[,"sigma[2]"]-ordMcmcMat[,"sigma[1]"] , compVal=compValSigmaDiff , 
                xlab=TeX(sprintf("$\\tau_{%s,%s} - \\tau_{%s,%s}$", subscript, groupNames[2], subscript, groupNames[1])) , 
                cex=cex_plotPost , col="#0088aa8a"  )
    } else {
      plotPost( ordMcmcMat[,"mu[1]"] , 
                xlab=TeX(sprintf("$\\mu_{%s}$", groupNames[1])) , 
                xlim=muLim , compVal=compValMu[1] , 
                col="#0088aa8a" , cex=cex_plotPost )
      plotPost( ordMcmcMat[,"sigma[1]"] , 
                xlab=TeX(sprintf("$\\tau_{%s}$", groupNames[1])) , 
                xlim=sigmaLim , compVal=compValSigma[1] , 
                col="#0088aa8a" , cex=cex_plotPost )
      plotPost( ordMcmcMat[,"mu[2]"] , 
                xlab=TeX(sprintf("$\\mu_{%s}$", groupNames[2])) , 
                xlim=muLim , compVal=compValMu[2] , 
                col="#0088aa8a" , cex=cex_plotPost )
      plotPost( ordMcmcMat[,"sigma[2]"] , 
                xlab=TeX(sprintf("$\\tau_{%s}$", groupNames[2])) , 
                xlim=sigmaLim , compVal=compValSigma[2] , 
                col="#0088aa8a" , cex=cex_plotPost )
      plotPost( ordMcmcMat[,"mu[2]"]-ordMcmcMat[,"mu[1]"] , compVal=compValMuDiff ,
                xlab=TeX(sprintf("$\\mu_{%s} - \\mu_{%s}$", groupNames[2], groupNames[1])) , 
                cex=cex_plotPost , col="#0088aa8a" )
      plotPost( ordMcmcMat[,"sigma[2]"]-ordMcmcMat[,"sigma[1]"] , compVal=compValSigmaDiff , 
                xlab=TeX(sprintf("$\\tau_{%s} - \\tau_{%s}$", groupNames[2], groupNames[1])) , 
                cex=cex_plotPost , col="#0088aa8a"  )
    }
  } 
  saveGraph( paste0(saveName,"-MuSigma") , type=graphFileType )
  
  # Posterior nu of one or two groups:
  if (nG == 1) {
    # One group:
    openGraph(width=3.5,height=4)
    if (subscript != "") {
      plotPost( log10(ordMcmcMat[,"nu"]) , 
                xlab=TeX(sprintf("$\\log_{10}\\nu_{%s}$", subscript)) ,
                compVal=compValNu , ROPE=NULL , cex=cex_plotPost , col="#0088aa8a")
    } else {
      plotPost( log10(ordMcmcMat[,"nu"]) , 
                xlab=TeX(sprintf("$\\log_{10}\\nu$")) ,
                compVal=compValNu , ROPE=NULL , cex=cex_plotPost , col="#0088aa8a")
    }
  }
  if (nG == 2) {
    # Two groups:
    openGraph(width=3.5*nG,height=4)
    par( mar=c(3.5,3.5,3,1) , mgp=c(2.25,0.7,0) )
    layout( matrix(1:2,nrow=1,byrow=TRUE) )
    if (subscript != "") {
      plotPost( log10(ordMcmcMat[,"nu[1]"]) , 
                xlab=TeX(sprintf("$\\log_{10}\\nu_{%s,%s}$", subscript, groupNames[1])) , 
                compVal=compValNu[1] , ROPE=NULL , cex=cex_plotPost , col="#0088aa8a")
      plotPost( log10(ordMcmcMat[,"nu[2]"]) , 
                xlab=TeX(sprintf("$\\log_{10}\\nu_{%s,%s}$", subscript, groupNames[2])) , 
                compVal=compValNu[2] , ROPE=NULL , cex=cex_plotPost , col="#0088aa8a")
    } else {
      plotPost( log10(ordMcmcMat[,"nu[1]"]) , 
                xlab=TeX(sprintf("$\\log_{10}\\nu_{%s}$", groupNames[1])) , 
                compVal=compValNu[1] , ROPE=NULL , cex=cex_plotPost , col="#0088aa8a")
      plotPost( log10(ordMcmcMat[,"nu[2]"]) , 
                xlab=TeX(sprintf("$\\log_{10}\\nu_{%s}$", groupNames[2])) , 
                compVal=compValNu[2] , ROPE=NULL , cex=cex_plotPost , col="#0088aa8a")
    }
  }
  saveGraph( paste0(saveName,"-nu") , type=graphFileType )
  
  # Posterior effect size:
  if (plotEffsz) {
    if ( nG == 1 ) {
      # One group:
      postEffSz = ( ordMcmcMat[,"mu"] - nullValEff ) / ordMcmcMat[,"sigma"]
      openGraph(width=3.5,height=4)
      if (subsEffsz != "" && subscript != "") {
        plotPost( postEffSz , 
                  xlab=TeX(sprintf("$d_{%s,%s}$", subsEffsz, subscript)) , 
                  compVal=compValEff , ROPE=ropeEffSz ,
                  col="#0088aa8a" , cex=cex_plotPost )
      }
      if (subsEffsz == "" && subscript != "") {
        plotPost( postEffSz , 
                  xlab=TeX(sprintf("$d_{%s}$", subscript)) , 
                  compVal=compValEff , ROPE=ropeEffSz ,
                  col="#0088aa8a" , cex=cex_plotPost )
      }
      if (subsEffsz != "" && subscript == "") {
        plotPost( postEffSz , 
                  xlab=TeX(sprintf("$d_{%s}$", subsEffsz)) , 
                  compVal=compValEff , ROPE=ropeEffSz ,
                  col="#0088aa8a" , cex=cex_plotPost )
      }
      if (subsEffsz == "" && subscript == "") {
        plotPost( postEffSz , 
                  xlab=TeX(sprintf("$d$")) , 
                  compVal=compValEff , ROPE=ropeEffSz ,
                  col="#0088aa8a" , cex=cex_plotPost )
      }
    }
    if ( nG == 2 ) {
      # Two groups:
      postEffSz = ( ( ordMcmcMat[,"mu[2]"] - ordMcmcMat[,"mu[1]"] ) 
                    / sqrt( ( ordMcmcMat[,"sigma[1]"]^2 + ordMcmcMat[,"sigma[2]"]^2 ) / 2 ) )
      openGraph(width=7,height=4)
      if (subsEffsz != "" && subscript != "") {
        plotPost( postEffSz , compVal=compValEff , ROPE=ropeEffSz , 
                  col="#0088aa8a" , cex=cex_plotPost , cex.lab=1.2,
                  xlab=TeX(sprintf("$d_{%s,%s}$", subsEffsz, subscript)) )
      }
      if (subsEffsz == "" && subscript != "") {
        plotPost( postEffSz , compVal=compValEff , ROPE=ropeEffSz , 
                  col="#0088aa8a" , cex=cex_plotPost , cex.lab=1.2,
                  xlab=TeX(sprintf("$d_{%s}$", subscript)) )
      }
      if (subsEffsz != "" && subscript == "") {
        plotPost( postEffSz , compVal=compValEff , ROPE=ropeEffSz , 
                  col="#0088aa8a" , cex=cex_plotPost , cex.lab=1.2,
                  xlab=TeX(sprintf("$d_{%s}$", subsEffsz)) )
      }
      if (subsEffsz == "" && subscript == "") {
        plotPost( postEffSz , compVal=compValEff , ROPE=ropeEffSz , 
                  col="#0088aa8a" , cex=cex_plotPost , cex.lab=1.2,
                  xlab=TeX(sprintf("$d$")) )
      }
    }
    saveGraph( paste0(saveName,"-EffSz") , type=graphFileType )
  }
  
  # Thresholds:
  openGraph(height=min(2.5*nQ,14),width=7)
  par( mar=c(3.5,3.5,2,1) , mgp=c(2.25,0.7,0) )
  layout( matrix(1:nQ,nrow=nQ) )
  for ( qIdx in 1:nQ ) {
    threshCols = paste0("thresh[",qIdx,",",1:(maxLevel-1),"]")
    threshMean = rowMeans( ordMcmcMat[,threshCols] )
    xLim = range(ordMcmcMat[,threshCols])
    nPtToPlot = 1000
    plotIdx = floor(seq(1,nrow(ordMcmcMat),length=nPtToPlot))
    plot( ordMcmcMat[plotIdx,threshCols[1]] , threshMean[plotIdx] , col="#0088aa8a" ,
          xlim=xLim , xlab="Thresholds" , ylab="Mean of thresholds" , 
          main=paste0("Item ",qIdx) )
    abline(v=mean(ordMcmcMat[plotIdx,threshCols[1]]),lty="dashed",col="#0088aa")
    axis(side=1, mean(ordMcmcMat[plotIdx,threshCols[1]]), 
         col="#0088aa", col.axis="#0088aa", 
         labels = round(mean(ordMcmcMat[plotIdx,threshCols[1]]), 1))
    for ( i in 2:length(threshCols) ) {
        points( ordMcmcMat[plotIdx,threshCols[i]] , threshMean[plotIdx] , col="#0088aa8a" )
      abline(v=mean(ordMcmcMat[plotIdx,threshCols[i]]),lty="dashed",col="#0088aa")
      axis(side=1, mean(ordMcmcMat[plotIdx,threshCols[i]]), 
           col="#0088aa", col.axis="#0088aa", 
           labels = round(mean(ordMcmcMat[plotIdx,threshCols[i]]), 1))
    }
  }  
  saveGraph( paste0(saveName,"-Thresh") , type=graphFileType )
  
  # Posterior predictive (histograms of ordinal responses with posterior 
  # predicted probabilities superimposed): 
  dataScaleDensMax = 0.7
  openGraph(height=min(2.5*nQ+0.75,14),width=3.5*nG)
  par( mar=c(3.5,3.5,2.5,0.5) , mgp=c(2.0,0.7,0) , oma=c(0,0,3.5,0) ) # , xpd=NA )
  layout( matrix(1:(nG*nQ),nrow=nQ,ncol=nG,byrow=FALSE) )
  for ( gIdx in 1:nG ) {
    for ( qIdx in 1:nQ ) {
      probInfo = NULL
      # Data histogram:
      thisY = y[ g==gIdx & q==qIdx ]
      xLim = c( minLevel-0.5 , maxLevel+0.5 )
      xBreaks = seq( xLim[1] , xLim[2] , 1 )  
      if (nG == 1) {
        if (qIdx < nQ) {
          # Only y label
          histInfo = hist( thisY , prob=TRUE , ylim=c(0,dataScaleDensMax) , 
                           xlim=xLim , breaks=xBreaks ,
                           xlab="" , ylab="p(response)" , cex.axis=1.2 , 
                           cex.lab=1.5 ,
                           col="#f4e3d7ff" , border="white" , xaxt="n" ,
                           main=paste0("Item ",qIdx), cex.main=1.5)
        } else {
          # Both labels
          histInfo = hist( thisY , prob=TRUE , ylim=c(0,dataScaleDensMax) , 
                           xlim=xLim , breaks=xBreaks ,
                           xlab="Ordinal responses" , ylab="p(response)" , 
                           cex.axis=1.2, cex.lab=1.5 , # yaxt="n" ,
                           col="#f4e3d7ff" , border="white" , xaxt="n" ,
                           main=paste0("Item ",qIdx), cex.main=1.5)
        }
        extra = extraInfo
      } else {
        if (gIdx > 1 & qIdx == nQ) {
          # Only x label
          histInfo = hist( thisY , prob=TRUE , ylim=c(0,dataScaleDensMax) , 
                           xlim=xLim , breaks=xBreaks ,
                           xlab="Ordinal responses" , ylab="" , cex.axis=1.2,
                           cex.lab=1.5 , yaxt="n" , 
                           col="#f4e3d7ff" , border="white" , xaxt="n" ,
                           main=paste0("Item ",qIdx,", Group ",groupNames[gIdx]),
                           cex.main=1.5)
        } else if (gIdx == 1 & qIdx < nQ) {
          # Only y label
          histInfo = hist( thisY , prob=TRUE , ylim=c(0,dataScaleDensMax) , 
                           xlim=xLim , breaks=xBreaks ,
                           xlab="" , ylab="p(response)" , cex.axis=1.2,
                           cex.lab=1.5 , # yaxt="n" ,
                           col="#f4e3d7ff" , border="white" , xaxt="n" ,
                           main=paste0("Item ",qIdx,", Group ",groupNames[gIdx]),
                           cex.main=1.5)
        } else if (gIdx > 1 & qIdx < nQ) {
          # No labels
          histInfo = hist( thisY , prob=TRUE , ylim=c(0,dataScaleDensMax) , 
                           xlim=xLim , breaks=xBreaks ,
                           xlab="" , ylab="" , cex.axis=1.2, 
                           cex.lab=1.5 , yaxt="n" ,
                           col="#f4e3d7ff" , border="white" , xaxt="n" ,
                           main=paste0("Item ",qIdx,", Group ",groupNames[gIdx]),
                           cex.main=1.5)
        } else {
          # Both labels
          histInfo = hist( thisY , prob=TRUE , ylim=c(0,dataScaleDensMax) , 
                           xlim=xLim , breaks=xBreaks ,
                           xlab="Ordinal responses" , ylab="p(response)" , 
                           cex.axis=1.2, cex.lab=1.5 , #yaxt="n" ,
                           col="#f4e3d7ff" , border="white" , xaxt="n" ,
                           main=paste0("Item ",qIdx,", Group ",groupNames[gIdx]),
                           cex.main=1.5)
        }
        #histInfo = hist( thisY , prob=TRUE , ylim=c(0,dataScaleDensMax) , 
        #                 xlim=xLim , breaks=xBreaks ,
        #                 xlab="Ordinal responses" , ylab="p(resp)" , cex.lab=1.5 , 
        #                 col="#f4e3d7ff" , border="white" , xaxt="n" ,
        #                 #yaxt="n" , 
        #                 main=paste0("Item ",qIdx,", Group ",groupNames[gIdx]) )
        extra = extraInfo[[gIdx]]
      }
      # Adding * in axis labels to indicate levels with extra data:
      if (!is.null(extraInfo)) {
        ticksLabels = NULL
        for (lv in 1:maxLevel) {
          if (extra[qIdx, lv] == 0)
            ticksLabels = cbind(ticksLabels, lv)
          else
            ticksLabels = cbind(ticksLabels, paste0(lv,"*"))
        }
      } else {
        ticksLabels = seq(minLevel, maxLevel, 1)
      }
      axis(side=1, seq(minLevel, maxLevel, 1), ticksLabels)
      
      # Posterior predicted probabilities:
      outProb=matrix(0,nrow=ordChainLength,ncol=maxLevel)
      require("LaplacesDemon")
      for ( stepIdx in 1:ordChainLength ) {
        if ( nG > 1 ) {
          threshCumProb = LaplacesDemon::pst(
            ordMcmcMat[ stepIdx , paste0("thresh[",qIdx,",",1:(maxLevel-1),"]") ] ,
            ordMcmcMat[ stepIdx , paste0("mu[",gIdx,"]") ] ,
            ordMcmcMat[ stepIdx , paste0("sigma[",gIdx,"]") ] ,
            ordMcmcMat[ stepIdx , paste0("nu[",gIdx,"]") ])
        } else { # if nG == 1
          threshCumProb = LaplacesDemon::pst(
            ordMcmcMat[ stepIdx , paste0("thresh[",qIdx,",",1:(maxLevel-1),"]") ] ,
            ordMcmcMat[ stepIdx , paste0("mu") ], 
            ordMcmcMat[ stepIdx , paste0("sigma") ], 
            ordMcmcMat[ stepIdx , paste0("nu") ])
        }
        outProb[stepIdx,] = c(threshCumProb,1) - c(0,threshCumProb)
      }
      outHdi = apply( outProb , 2 , HDIofMCMC )
      outMedian = apply( outProb , 2 , median , na.rm=TRUE )
      show(outMedian)
      points( x=1:maxLevel , y=outMedian  , pch=19 , cex=1.5 , col="#0088aab3" )
      segments( x0=1:maxLevel , y0=outHdi[1,] , 
                x1=1:maxLevel , y1=outHdi[2,] , lwd=4 , col="#0088aab3" )
      
      probInfo = cbind(probInfo, Median=outMedian)
      probInfo = cbind(probInfo, Mean=apply( outProb , 2 , mean , na.rm=TRUE ))
      probInfo = cbind(probInfo, HDIlow=outHdi[1,])
      probInfo = cbind(probInfo, HDIhigh=outHdi[2,])
      if ( !is.null(saveName) ) {
        write.csv( probInfo , file=paste0(saveName,"-ProbItem",qIdx,"-Group",gIdx,".csv") )
      }
    }
  }
  saveGraph( paste0(saveName,"-PostPred") , type=graphFileType )
  
  # Parameters pairwise, to see correlations:
  if (nG == 1) {
    openGraph(width=7*3/5,height=7*3/5)
    nPtToPlot = 1000
    plotIdx = floor(seq(1,NROW(ordMcmcMat),by=NROW(ordMcmcMat)/nPtToPlot))
    panel.cor = function(x, y, digits=2, prefix="", cex.cor, ...) {
      usr = par("usr"); on.exit(par(usr))
      par(usr = c(0, 1, 0, 1))
      r = (cor(x, y))
      txt = format(c(r, 0.123456789), digits=digits)[1]
      txt = paste(prefix, txt, sep="")
      if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
      text(0.5, 0.5, txt, cex=1.25 ) # was cex=cex.cor*r
    }
    if (subscript != "") {
      pairs( cbind( ordMcmcMat[,"mu"] , ordMcmcMat[,"sigma"] , 
                    log10(ordMcmcMat[,"nu"]) )[plotIdx,] ,
             labels=c( TeX(sprintf("$\\mu_{%s}$", subscript)) , 
                       TeX(sprintf("$\\tau_{%s}$", subscript)) , 
                       TeX(sprintf("$\\log_{10}\\nu_{%s}$", subscript)) ) , 
             lower.panel=panel.cor , col="#0088aa8a" )
    } else {
      pairs( cbind( ordMcmcMat[,"mu"] , ordMcmcMat[,"sigma"] , 
                    log10(ordMcmcMat[,"nu"]) )[plotIdx,] ,
             labels=c( TeX(sprintf("$\\mu$")) , 
                       TeX(sprintf("$\\tau$")) , 
                       TeX(sprintf("$\\log_{10}\\nu$")) ) , 
             lower.panel=panel.cor , col="#0088aa8a" )
    }
  }
  if (nG == 2) {
    openGraph(width=7,height=7)
    nPtToPlot = 1000
    plotIdx = floor(seq(1,length(ordMcmcMat[,"mu[1]"]),
                        by=length(ordMcmcMat[,"mu[1]"])/nPtToPlot))
    panel.cor = function(x, y, digits=2, prefix="", cex.cor, ...) {
      usr = par("usr"); on.exit(par(usr))
      par(usr = c(0, 1, 0, 1))
      r = (cor(x, y))
      txt = format(c(r, 0.123456789), digits=digits)[1]
      txt = paste(prefix, txt, sep="")
      if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
      text(0.5, 0.5, txt, cex=1.25 ) # was cex=cex.cor*r
    }
    if (subscript != "") {
      pairs( cbind( ordMcmcMat[,"mu[1]"] , ordMcmcMat[,"mu[2]"] , 
                    ordMcmcMat[,"sigma[1]"], ordMcmcMat[,"sigma[2]"] , 
                    log10(ordMcmcMat[,"nu[1]"]), log10(ordMcmcMat[,"nu[2]"]))[plotIdx,] ,
             labels=c( TeX(sprintf("$\\mu_{%s,%s}$", subscript, groupNames[1])) , 
                       TeX(sprintf("$\\mu_{%s,%s}$", subscript, groupNames[2])) , 
                       TeX(sprintf("$\\tau_{%s,%s}$", subscript, groupNames[1])) , 
                       TeX(sprintf("$\\tau_{%s,%s}$", subscript, groupNames[2])) , 
                       TeX(sprintf("$\\log_{10}\\nu_{%s,%s}$", subscript, groupNames[1])) , 
                       TeX(sprintf("$\\log_{10}\\nu_{%s,%s}$", subscript, groupNames[2])) ) , 
             lower.panel=panel.cor , col="#0088aa8a" )
    } else {
      pairs( cbind( ordMcmcMat[,"mu[1]"] , ordMcmcMat[,"mu[2]"] , 
                    ordMcmcMat[,"sigma[1]"], ordMcmcMat[,"sigma[2]"] , 
                    log10(ordMcmcMat[,"nu[1]"]), log10(ordMcmcMat[,"nu[2]"]))[plotIdx,] ,
             labels=c( TeX(sprintf("$\\mu_{%s}$", groupNames[1])) , 
                       TeX(sprintf("$\\mu_{%s}$", groupNames[2])) , 
                       TeX(sprintf("$\\tau_{%s}$", groupNames[1])) , 
                       TeX(sprintf("$\\tau_{%s}$", groupNames[2])) , 
                       TeX(sprintf("$\\log_{10}\\nu_{%s}$", groupNames[1])) , 
                       TeX(sprintf("$\\log_{10}\\nu_{%s}$", groupNames[2])) ) , 
             lower.panel=panel.cor , col="#0088aa8a" )
    }
  }
  saveGraph( file=paste(saveName,"-PostPairs",sep=""), type=graphFileType)
  
}