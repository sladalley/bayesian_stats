# Functions for the Bayesian analysis of metric data with robust model.
# Based on Chapter 16 of Kruschke (2015) and the scripts accompanying it.
#
# * Kruschke, J. K. (2015). Doing Bayesian Data Analysis, Second Edition: 
#   A Tutorial with R, JAGS, and Stan. Academic Press / Elsevier.

source("DBDA2E-utilities.R")
library('latex2exp')

# ====== MCMC chain ===========================================================
genMCMC = function(datFrm, yName="y", gName="cond", muPrior, muSdPrior, 
                   sigmaPriorLow, sigmaPriorHigh, saveName=NULL) { 
  # Generate the MCMC chain using JAGS.
  # List of parameters:
  #   - datFrm: data set.
  #   - yName: name of the column with the data.
  #   - gName: name of the column with the group indexes.
  #   - muPrior: mean of the prior for the mean of each group (normal dist.).
  #   - muSdPrior: std. dev. of the prior for the mean of each group (normal dist.).
  #   - sigmaPriorLow: low limit of the prior for the scale of each group (uniform dist.).
  #   - sigmaPriorHigh: high limit of the prior for the scale of each group (uniform dist.).
  #   - saveName: prefix of the output files.
  #-----------------------------------------------------------------------------

  # Creating group column if not existing:
  if (is.null(gName)) {
    datFrm = cbind(datFrm, group=rep(1,dim(datFrm)[1]))
    gName = "group"
  }
  
  # ASSEMBLE THE DATA FOR JAGS.
  y = as.numeric(datFrm[,yName])
  # Do some checking that data make sense:
  if ( any( !is.finite(y) ) ) { stop("All y values must be finite.") }
 
  Ntotal = length(y)
  g = as.numeric(as.factor(datFrm[,gName])) # group data vector
  gLevels = levels(as.factor(datFrm[,gName]))
  nG = max(g) # number of groups
  
  # Specify the data in a list, for later shipment to JAGS:
  dataList = list(
    y = y , # data
    g = g , # groups
    nG = nG , # number of groups
    Ntotal = Ntotal , # total size of data set
    muPrior = muPrior , # mean of the normal prior of the mean
    muSdPrior = muSdPrior , # std. dev. of the normal prior of the mean
    sigmaPriorLow = sigmaPriorLow , # low limit of the uniform prior of sigma
    sigmaPriorHigh = sigmaPriorHigh # high limit of the uniform prior of sigma
  )
  #-----------------------------------------------------------------------------
  # THE METRIC ROBUST MODEL.
  modelString = "
  model {
    for ( i in 1:Ntotal ) {
      y[i] ~ dt( mu[g[i]] , 1/sigma[g[i]]^2 , nu[g[i]] )
    }
    for ( gIdx in 1:nG ) { 
      mu[gIdx] ~ dnorm( muPrior[gIdx] , 1/(muSdPrior[gIdx])^2 )
      sigma[gIdx] ~ dunif( sigmaPriorLow[gIdx] , sigmaPriorHigh[gIdx] )
      nu[gIdx] ~ dexp(1/30.0)
    }
  }
  " # close quote for modelString
  # Write out modelString to a text file
  writeLines( modelString , con=paste0(saveName, "-Model.txt") )
  #-----------------------------------------------------------------------------
  # INTIALIZE THE CHAINS.
  # Initial values of MCMC chains based on data:
  #   The start of the MCMC chains is defined according to the data. 
  #   The estimated paremeters are: mean, standard deviation and normality.
  #   The initial mean and s.d. are the ones obtained from the data and
  #     the normality initializes at 5, as a moderate value.
  if (nG == 1) { 
    initsList = list( mu = mean(y) , sigma = sd(y) , nuMinusOne = 5 )
  }
  if (nG == 2) {
    initsList = list( mu = c( mean(y[g==1]) , mean(y[g==2]) ) , 
                      sigma = c( sd(y[g==1]) , sd(y[g==2]) ) , 
                      nuMinusOne = c( 5 , 5 ) )
  }
  #-----------------------------------------------------------------------------
  # RUN THE CHAINS
  parameters = c( "mu" , "sigma" , "nu" )     # The parameters to be monitored
  numSavedSteps = 20000
  thinSteps = 5
  adaptSteps = 500               # Number of steps to "tune" the samplers
  burnInSteps = 1000
  runjagsMethod=runjagsMethodDefault # from DBDA2E-utilities
  nChains=nChainsDefault # from DBDA2E-utilities
  
  runJagsOut <- run.jags( method=runjagsMethod ,
                          model=paste0(saveName, "-Model.txt") , 
                          monitor=parameters , 
                          data=dataList ,  
                          inits=initsList , 
                          n.chains=nChains ,
                          adapt=adaptSteps ,
                          burnin=burnInSteps , 
                          sample=ceiling(numSavedSteps/nChains) ,
                          thin=thinSteps ,
                          summarise=FALSE ,
                          plots=FALSE )
  codaSamples = as.mcmc.list( runJagsOut )
  # resulting codaSamples object has these indices: 
  #   codaSamples[[ chainIdx ]][ stepIdx , paramIdx ]
  if ( !is.null(saveName) ) {
    save( codaSamples , file=paste(saveName,"-Mcmc.Rdata",sep="") )
  }
  
  return(codaSamples)
}

# ====== Summary of data and diagnostics ======================================
smryMCMC = function (codaSamples, nG, nullValEff=0, saveName=NULL, 
                     diagnostics=TRUE, graphFileType="png", computeEffsz=TRUE) {
  # Compute summary statistics of the chain and generate diagnostics plot.
  # List of parameters:
  #   - codaSamples: codaSamples object with the MCMC chain.
  #   - nG: number of groups.
  #   - nullValEff: null value to calculate effect size when single group.
  #   - saveName: prefix of the output files.
  #   - diagnostics: to define if it should generate diagnostics plots.
  #   - graphFileType: type of the image output files.
  #   - computeEffsz: if effect size should be calculated (TRUE) or not (FALSE)
  # --------------------------------------------------------------------------- #
  # SUMMARIZING THE MCMC POSTERIOR.
  
  summaryInfo = NULL
  mcmcMat = as.matrix(codaSamples, chains=TRUE)
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
    # Difference between means, mu1 - mu2:
    summaryInfo = rbind( summaryInfo , 
                         "muDiff" = summarizePost( 
                           mcmcMat[,"mu[1]"] - mcmcMat[,"mu[2]"] , 
                           compVal=NULL , ROPE=NULL ) )
    # Scale parameter of group 1:
    summaryInfo = rbind( summaryInfo , 
                         "sigma[1]" = summarizePost( mcmcMat[,"sigma[1]"] , 
                                                     compVal=NULL , ROPE=NULL ) )
    # Scale parameter of group 2:
    summaryInfo = rbind( summaryInfo , 
                         "sigma[2]" = summarizePost( mcmcMat[,"sigma[2]"] , 
                                                     compVal=NULL , ROPE=NULL ) )
    # Difference between scales, sigma1 - sigma2:
    summaryInfo = rbind( summaryInfo , 
                         "sigmaDiff" = summarizePost( 
                           mcmcMat[,"sigma[1]"] - mcmcMat[,"sigma[2]"] , 
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
    # Effect size, (mu1 - mu2) / sqrt((sigma1^2 + sigma2^2)/2):
    if (computeEffsz) {
      summaryInfo = rbind( summaryInfo , 
                           "effSz" = summarizePost( 
                             ( mcmcMat[,"mu[1]"]-mcmcMat[,"mu[2]"] ) 
                             / sqrt((mcmcMat[,"sigma[1]"]^2+mcmcMat[,"sigma[2]"]^2)/2) ,
                             compVal=NULL , ROPE=NULL ) )
    }
  }
  
  if ( !is.null(saveName) ) {
    write.csv( summaryInfo , file=paste0(saveName,"-SummaryInfo.csv") )
  }
  
  # --------------------------------------------------------------------------- #
  # DIAGNOSTICS
  
  if ( diagnostics == TRUE ) {
    splitted = strsplit(saveName, "/")
    diagnosticsPath = ""
    for (i in 1:(length(splitted[[1]])-1)) {
      diagnosticsPath = paste0(diagnosticsPath, splitted[[1]][i], "/")
    }
    diagnosticsPath = paste0(diagnosticsPath, "diagnostics")
    dir.create(diagnosticsPath)
    parameterNames = varnames(codaSamples) 
    for ( parName in parameterNames ) {
      diagMCMC( codaObject=codaSamples , parName=parName ,  
                saveName=paste0(diagnosticsPath, "/", sub(".*/", "", saveName), "-") , 
                saveType=graphFileType )
    }
  }
}

# ====== Posterior distributions ==============================================
plotMCMC = function( codaSamples, datFrm, yName="y", gName="cond", 
                     nullValEff=0, compValMu=NULL, compValMuDiff=NULL, 
                     compValSigma=NULL, compValSigmaDiff=NULL, compValNu=NULL, 
                     compValEff=NULL, ropeEffSz=NULL, graphFileType="png", 
                     saveName=NULL, groupNames=c(1,2), subscript="", subsEffsz="",
                     plotEffsz=TRUE) {
  # Display posterior information.
  # List of parameters:
  #   - codaSamples: codaSamples object with the MCMC chain.
  #   - datFrm: data set.
  #   - yName: name of the column with the data.
  #   - gName: name of the column with the group indexes. 
  #   - nullValEff: null value to calculate effect size when single group.
  #   - compValMu: comparison value in the posterior mu plot(s).
  #   - compValMuDiff: comparison value in the posterior of mu difference plot (two groups).
  #   - compValSigma: comparison value in the posterior sigma plot(s).
  #   - compValSigmaDiff: comparison value in the posterior of sigma difference plot (two groups).
  #   - compValNu: comparison value in the posterior nu plot(s).
  #   - compValEff: comparison value in the posterior effect size plot.
  #   - ropeEffSz: ROPE in the posterior effect size plot.
  #   - graphFileType: type of the image output files.
  #   - saveName: prefix of the output files.
  #   - groupNames: names of each group to be displayed in the plots.
  #   - subscript: subscript to indicate the variable.
  #   - subsEffsz: subscript to indicate the expression for the effect size.
  mcmcMat = as.matrix(codaSamples,chains=TRUE)
  chainLength = NROW( mcmcMat )
  
  # Creating group column if not existing:
  if (is.null(gName)) {
    datFrm = cbind(datFrm, group=rep(1,dim(datFrm)[1]))
    gName = "group"
  }
  
  y = as.numeric(datFrm[,yName])
  Ntotal = length(y)
  g = as.numeric(as.factor(datFrm[,gName])) # group data vector
  gLevels = levels(as.factor(datFrm[,gName]))
  nG = max(g) # number of groups
  
  cex_plotPost = 1.0
  # Posterior mu and sigma of one or two groups:
  if (nG == 1) {
    # One group:
    openGraph(height=2.5*nG,width=7)
    par( mar=c(3.5,3.5,3,1) , mgp=c(2.25,0.7,0) )
    layout( matrix(1:(2*nG),nrow=nG,byrow=TRUE) )
    muLim = range( mcmcMat[, grep("^mu",colnames(mcmcMat)) ] )
    sigmaLim = range( mcmcMat[, grep("^sigma",colnames(mcmcMat)) ] )
    if (subscript != "") {
      plotPost( mcmcMat[,"mu"] , xlab=TeX(sprintf("$\\mu_{%s}$", subscript)) , xlim=muLim , 
                compVal=compValMu , cex=cex_plotPost , col="#0088aa8a")
      plotPost( mcmcMat[,"sigma"] , xlab=TeX(sprintf("$\\tau_{%s}$", subscript)) , xlim=sigmaLim , 
                compVal=compValSigma , cex=cex_plotPost , col="#0088aa8a")
    } else {
      plotPost( mcmcMat[,"mu"] , xlab=TeX(sprintf("$\\mu$")) , xlim=muLim , 
                compVal=compValMu , cex=cex_plotPost , col="#0088aa8a")
      plotPost( mcmcMat[,"sigma"] , xlab=TeX(sprintf("$\\tau$")) , xlim=sigmaLim , 
                compVal=compValSigma , cex=cex_plotPost , col="#0088aa8a")
    }
  }
  if (nG == 2) {
    # Two groups:
    openGraph(height=2.5*(nG+1),width=7)
    par( mar=c(3.5,3.5,3,1) , mgp=c(2.25,0.7,0) )
    layout( matrix(1:(2*(nG+1)),nrow=nG+1,byrow=TRUE) )
    muLim = range( mcmcMat[, grep("^mu\\[",colnames(mcmcMat)) ] )
    sigmaLim = range( mcmcMat[, grep("^sigma\\[",colnames(mcmcMat)) ] )
    if (subscript != "") {
      plotPost( mcmcMat[,"mu[1]"] , 
                xlab=TeX(sprintf("$\\mu_{%s,%s}$", subscript, groupNames[1])) , 
                xlim=muLim , compVal=compValMu[1] , 
                col="#0088aa8a" , cex=cex_plotPost )
      plotPost( mcmcMat[,"sigma[1]"] , 
                xlab=TeX(sprintf("$\\tau_{%s,%s}$", subscript, groupNames[1])) , 
                xlim=sigmaLim , compVal=compValSigma[1] , 
                col="#0088aa8a" , cex=cex_plotPost )
      plotPost( mcmcMat[,"mu[2]"] , 
                xlab=TeX(sprintf("$\\mu_{%s,%s}$", subscript, groupNames[2])) , 
                xlim=muLim , compVal=compValMu[2] , 
                col="#0088aa8a" , cex=cex_plotPost )
      plotPost( mcmcMat[,"sigma[2]"] , 
                xlab=TeX(sprintf("$\\tau_{%s,%s}$", subscript, groupNames[2])) , 
                xlim=sigmaLim , compVal=compValSigma[2] , 
                col="#0088aa8a" , cex=cex_plotPost )
      plotPost( mcmcMat[,"mu[1]"]-mcmcMat[,"mu[2]"] , compVal=compValMuDiff ,
                xlab=TeX(sprintf("$\\mu_{%s,%s} - \\mu_{%s,%s}$", subscript, groupNames[1], subscript, groupNames[2])) , 
                cex=cex_plotPost , col="#0088aa8a" )
      plotPost( mcmcMat[,"sigma[1]"]-mcmcMat[,"sigma[2]"] , compVal=compValSigmaDiff , 
                xlab=TeX(sprintf("$\\tau_{%s,%s} - \\tau_{%s,%s}$", subscript, groupNames[1], subscript, groupNames[2])) , 
                cex=cex_plotPost , col="#0088aa8a"  )
    } else {
      plotPost( mcmcMat[,"mu[1]"] , 
                xlab=TeX(sprintf("$\\mu_{%s}$", groupNames[1])) , 
                xlim=muLim , compVal=compValMu[1] , 
                col="#0088aa8a" , cex=cex_plotPost )
      plotPost( mcmcMat[,"sigma[1]"] , 
                xlab=TeX(sprintf("$\\tau_{%s}$", groupNames[1])) , 
                xlim=sigmaLim , compVal=compValSigma[1] , 
                col="#0088aa8a" , cex=cex_plotPost )
      plotPost( mcmcMat[,"mu[2]"] , 
                xlab=TeX(sprintf("$\\mu_{%s}$", groupNames[2])) , 
                xlim=muLim , compVal=compValMu[2] , 
                col="#0088aa8a" , cex=cex_plotPost )
      plotPost( mcmcMat[,"sigma[2]"] , 
                xlab=TeX(sprintf("$\\tau_{%s}$", groupNames[2])) , 
                xlim=sigmaLim , compVal=compValSigma[2] , 
                col="#0088aa8a" , cex=cex_plotPost )
      plotPost( mcmcMat[,"mu[1]"]-mcmcMat[,"mu[2]"] , compVal=compValMuDiff ,
                xlab=TeX(sprintf("$\\mu_{%s} - \\mu_{%s}$", groupNames[1], groupNames[2])) , 
                cex=cex_plotPost , col="#0088aa8a" )
      plotPost( mcmcMat[,"sigma[1]"]-mcmcMat[,"sigma[2]"] , compVal=compValSigmaDiff , 
                xlab=TeX(sprintf("$\\tau_{%s} - \\tau_{%s}$", groupNames[1], groupNames[2])) , 
                cex=cex_plotPost , col="#0088aa8a"  )
    }
  } 
  saveGraph( paste0(saveName,"-MuSigma") , type=graphFileType )
  
  # Posterior nu of one or two groups:
  if (nG == 1) {
    # One group:
    openGraph(width=3.5,height=4)
    if (subscript != "") {
      plotPost( log10(mcmcMat[,"nu"]) , 
                xlab=TeX(sprintf("$\\log_{10}\\nu_{%s}$", subscript)) , 
                compVal=compValNu , ROPE=NULL , cex=cex_plotPost , col="#0088aa8a")
    } else {
      plotPost( log10(mcmcMat[,"nu"]) , 
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
      plotPost( log10(mcmcMat[,"nu[1]"]) , 
                xlab=TeX(sprintf("$\\log_{10}\\nu_{%s,%s}$", subscript, groupNames[1])) , 
                compVal=compValNu[1] , ROPE=NULL , cex=cex_plotPost , 
                main=paste0("Group ", groupNames[1]) , col="#0088aa8a")
      plotPost( log10(mcmcMat[,"nu[2]"]) , 
                xlab=TeX(sprintf("$\\log_{10}\\nu_{%s,%s}$", subscript, groupNames[2])) , 
                compVal=compValNu[2] , ROPE=NULL , cex=cex_plotPost , 
                main=paste0("Group ", groupNames[2]) ,col="#0088aa8a")
    } else {
      plotPost( log10(mcmcMat[,"nu[1]"]) , 
                xlab=TeX(sprintf("$\\log_{10}\\nu_{%s}$", groupNames[1])) , 
                compVal=compValNu[1] , ROPE=NULL , cex=cex_plotPost , 
                main=paste0("Group ", groupNames[1]) , col="#0088aa8a")
      plotPost( log10(mcmcMat[,"nu[2]"]) , 
                xlab=TeX(sprintf("$\\log_{10}\\nu_{%s}$", groupNames[2])) , 
                compVal=compValNu[2] , ROPE=NULL , cex=cex_plotPost , 
                main=paste0("Group ", groupNames[2]) ,col="#0088aa8a")
    }
  }
  saveGraph( paste0(saveName,"-nu") , type=graphFileType )
  
  # Posterior effect size:
  if (plotEffsz) {
    if ( nG == 1 ) {
      # One group:
      postEffSz = ( mcmcMat[,"mu"] - nullValEff ) / mcmcMat[,"sigma"]
      openGraph(width=3.5,height=4)
      if (subsEffsz != "" && subscript != "") {
        plotPost( postEffSz , xlab=TeX(sprintf("$d_{%s,%s}$", subsEffsz, subscript)) , 
                  compVal=compValEff , ROPE=ropeEffSz , 
                  col="#0088aa8a" , cex=cex_plotPost )
      } 
      if (subsEffsz == "" && subscript != "") {
        plotPost( postEffSz , xlab=TeX(sprintf("$d_{%s}$", subscript)) , 
                  compVal=compValEff , ROPE=ropeEffSz , 
                  col="#0088aa8a" , cex=cex_plotPost )
      } 
      if (subsEffsz != "" && subscript == "") {
        plotPost( postEffSz , xlab=TeX(sprintf("$d_{%s}$", subsEffsz)) , 
                  compVal=compValEff , ROPE=ropeEffSz , 
                  col="#0088aa8a" , cex=cex_plotPost )
      } 
      if (subsEffsz == "" && subscript == "") {
        plotPost( postEffSz , xlab=TeX(sprintf("$d$")) , 
                  compVal=compValEff , ROPE=ropeEffSz , 
                  col="#0088aa8a" , cex=cex_plotPost )
      }
    }
    if ( nG == 2 ) {
      # Two groups:
      postEffSz = ( ( mcmcMat[,"mu[1]"] - mcmcMat[,"mu[2]"] ) 
                    / sqrt( ( mcmcMat[,"sigma[1]"]^2 + mcmcMat[,"sigma[2]"]^2 ) / 2 ) )
      openGraph(width=7,height=4)
      if (subsEffsz != "" && subscript != "") {
        plotPost( postEffSz , compVal=compValEff , ROPE=ropeEffSz , 
                  col="#0088aa8a" , cex=cex_plotPost,
                  xlab=TeX(sprintf("$d_{%s,%s}$", subsEffsz, subscript)) )
      }
      if (subsEffsz == "" && subscript != "") {
        plotPost( postEffSz , compVal=compValEff , ROPE=ropeEffSz , 
                  col="#0088aa8a" , cex=cex_plotPost,
                  xlab=TeX(sprintf("$d_{%s}$", subscript)) )
      }
      if (subsEffsz != "" && subscript == "") {
        plotPost( postEffSz , compVal=compValEff , ROPE=ropeEffSz , 
                  col="#0088aa8a" , cex=cex_plotPost,
                  xlab=TeX(sprintf("$d_{%s}$", subsEffsz)) )
      }
      if (subsEffsz == "" && subscript == "") {
        plotPost( postEffSz , compVal=compValEff , ROPE=ropeEffSz , 
                  col="#0088aa8a" , cex=cex_plotPost,
                  xlab=TeX(sprintf("$d$")) )
      }
    }
    saveGraph( paste0(saveName,"-EffSz") , type=graphFileType )
  }
  # Posterior predictive (histogram of data with a smattering of curves 
  # superimposed):
  openGraph(height=min(2.5+0.75,14),width=3.5*nG)
  par( mar=c(3.5,3.5,2.5,0.5) , mgp=c(2.0,0.7,0) , oma=c(0,0,3.5,0) ) # , xpd=NA )
  layout( matrix(1:(nG),nrow=1,ncol=nG,byrow=FALSE) )
  # Select thinned steps in chain for plotting of posterior predictive curves:
  nCurvesToPlot = 20
  for ( gIdx in 1:nG ) {
    thisY = y[ g==gIdx ]
    xLim = c( min(thisY)-0.1*(max(thisY)-min(thisY)) , 
              max(thisY)+0.1*(max(thisY)-min(thisY)) )
    xBreaks = seq( xLim[1] , xLim[2] , 
                   length=ceiling((xLim[2]-xLim[1])/(sd(thisY)/4)) )
    histInfo = hist(thisY, breaks=xBreaks, plot=FALSE)
    yMax = 1.2 * max( histInfo$density )
    # Histogram:
    if (nG == 1) {
      histInfo = hist( thisY , prob=TRUE , ylim=c(0,yMax) , 
                       xlim=xLim , breaks=xBreaks , xlab="Data" , ylab="" , 
                       cex.lab=1.0 , col="#f4e3d7ff" , border="white" , yaxt="n" ,
                       main="")
    } else {
      histInfo = hist( thisY , prob=TRUE , ylim=c(0,yMax) , 
                       xlim=xLim , breaks=xBreaks , xlab="Data" , ylab="" , 
                       cex.lab=1.0 , col="#f4e3d7ff" , border="white" , yaxt="n" ,
                       main=paste0("Group ", groupNames[gIdx]))
    }
    # Posterior predictive curves:
    stepIdxVec = seq( 1 , chainLength , floor(chainLength/nCurvesToPlot) )
    if (nG == 1) {
      for ( stepIdx in 1:length(stepIdxVec) ) {
        xVec = seq( xLim[1] , xLim[2] , length=501 )
        lines(xVec, dt( (xVec-mcmcMat[ stepIdxVec[stepIdx] , "mu" ]) /
                          mcmcMat[ stepIdxVec[stepIdx] , "sigma" ], 
                        df=mcmcMat[ stepIdxVec[stepIdx] , "nu" ] ) /
                mcmcMat[ stepIdxVec[stepIdx] , "sigma" ] , 
              type="l" , col="#0088aa8a" , lwd=1 )
      }
    }
    if (nG > 1) {
      for ( stepIdx in 1:length(stepIdxVec) ) {
        xVec = seq( xLim[1] , xLim[2] , length=501 )
        lines(xVec, dt( (xVec-mcmcMat[ stepIdxVec[stepIdx] , paste0("mu[",gIdx,"]") ]) /
                          mcmcMat[ stepIdxVec[stepIdx] , paste0("sigma[",gIdx,"]") ], 
                        df=mcmcMat[ stepIdxVec[stepIdx] , paste0("nu[",gIdx,"]") ] ) /
                mcmcMat[ stepIdxVec[stepIdx] , paste0("sigma[",gIdx,"]") ] , 
              type="l" , col="#0088aa8a" , lwd=1 )
      }
    }
  }
  mtext( text=TeX(sprintf("$\\bar{x}=%f\\;\\,\\tilde{x}=%f\\;\\,s=%f$", round(mean(y), 2), 
                          round(median(y), 2), round(sd(y), 2))) ,
         at=c(0.5) , cex=ifelse(nG>=2,1.5,1.0) , side=3 , outer=TRUE , 
         adj=c(0.5,0.5) , padj=c(-0.5,-0.5) )
  saveGraph( paste0(saveName,"-PostPred") , type=graphFileType )
  
  # Parameters pairwise, to see correlations:
  if (nG == 1) {
    openGraph(width=7*3/5,height=7*3/5)
    nPtToPlot = 1000
    plotIdx = floor(seq(1,NROW(mcmcMat),by=NROW(mcmcMat)/nPtToPlot))
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
      pairs( cbind( mcmcMat[,"mu"] , mcmcMat[,"sigma"] , 
                    log10(mcmcMat[,"nu"]) )[plotIdx,] ,
             labels=c( TeX(sprintf("$\\mu_{%s}$", subscript)) , 
                       TeX(sprintf("$\\tau_{%s}$", subscript)) , 
                       TeX(sprintf("$\\log_{10}\\nu_{%s}$", subscript)) ) , 
             lower.panel=panel.cor , col="#0088aa8a" )
    } else {
      pairs( cbind( mcmcMat[,"mu"] , mcmcMat[,"sigma"] , 
                    log10(mcmcMat[,"nu"]) )[plotIdx,] ,
             labels=c( TeX(sprintf("$\\mu$")) , 
                       TeX(sprintf("$\\tau$")) , 
                       TeX(sprintf("$\\log_{10}\\nu$")) ) , 
             lower.panel=panel.cor , col="#0088aa8a" )
    }
  }
  if (nG == 2) {
    openGraph(width=7,height=7)
    nPtToPlot = 1000
    plotIdx = floor(seq(1,length(mcmcMat[,"mu[1]"]),
                        by=length(mcmcMat[,"mu[1]"])/nPtToPlot))
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
      pairs( cbind( mcmcMat[,"mu[1]"] , mcmcMat[,"mu[2]"] , 
                    mcmcMat[,"sigma[1]"], mcmcMat[,"sigma[2]"] , 
                    log10(mcmcMat[,"nu[1]"]), log10(mcmcMat[,"nu[2]"]))[plotIdx,] ,
             labels=c( TeX(sprintf("$\\mu_{%s,%s}$", subscript, groupNames[1])) , 
                       TeX(sprintf("$\\mu_{%s,%s}$", subscript, groupNames[2])) , 
                       TeX(sprintf("$\\tau_{%s,%s}$", subscript, groupNames[1])) , 
                       TeX(sprintf("$\\tau_{%s,%s}$", subscript, groupNames[2])) , 
                       TeX(sprintf("$\\log_{10}\\nu_{%s,%s}$", subscript, groupNames[1])) , 
                       TeX(sprintf("$\\log_{10}\\nu_{%s,%s}$", subscript, groupNames[2])) ) , 
             lower.panel=panel.cor , col="#0088aa8a" )
    } else {
      pairs( cbind( mcmcMat[,"mu[1]"] , mcmcMat[,"mu[2]"] , 
                    mcmcMat[,"sigma[1]"], mcmcMat[,"sigma[2]"] , 
                    log10(mcmcMat[,"nu[1]"]), log10(mcmcMat[,"nu[2]"]))[plotIdx,] ,
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

#===============================================================================
# Functions to estimate power in the robust bayesian analyses of metric data.
#   The functions are adapted from Jags-Ydich-Xnom1subj-MbernBeta-Power.R and
#   BEST.R (from *Bayesian estimation supersedes the t test* [Kruschke, 2013]) 
#   scripts.

goalAchievedForSample = function( data, muNULL, effROPE, effHDImaxWid,
                                  mcmcLength=10000) {
  # Generate the MCMC chain:
  mcmcCoda = genMCMC( data=data , numSavedSteps=mcmcLength , saveName=NULL )
  mcmcMat = as.matrix(mcmcCoda,chains=TRUE)
  # Calculate effect size:
  mcmcMatEff = ( mcmcMat[,"mu"] - muNULL ) / mcmcMat[,"sigma"]
  
  # Check goal achievement. First, compute the HDI of the effect size:
  effHDI = HDIofMCMC( mcmcMatEff )
  
  # Define list for recording results:
  goalAchieved = list()
  
  # All the goals are related to the effect size posterior distribution.
  # Goal 1: Exclude ROPE around null value:
  goalAchieved = c( goalAchieved , 
                    "ExcludeROPE"=( effHDI[1] > effROPE[2] 
                                    | effHDI[2] < effROPE[1] ) )
  
  # Goal 2: HDI less than max width:
  goalAchieved = c( goalAchieved , 
                    "NarrowHDI"=( effHDI[2]-effHDI[1] < effHDImaxWid ) )
  
  # Goal 3: HDI all above the ROPE around null value:
  goalAchieved = c( goalAchieved , 
                    "HDIaboveROPE"=( effHDI[1] > effROPE[2] ) )
  
  # Return list of goal results:
  return(goalAchieved)
}

powerEstimation = function( mcmcChain, N, muNULL, 
                            effROPE, effHDImaxWid, 
                            mcmcLength=10000, nRep=1000 , 
                            saveName=NULL, recover=0) {
  # Description of arguments:
  #   - mcmcChain is a matrix with a MCMC chain.
  #   - N is the sample size.
  #   - muNULL is the mean null value.
  #   - effROPE is a two element vector, such as c(-1,1), specifying the limit
  #     of the ROPE on the effect size.
  #   - effHDImaxWid is the maximum desired width of the 95% HDI on the effect size.
  #   - numSavedSteps is the number of steps in the chains generated with simulated data.
  #   - nRep is the number of simulated experiments used to estimate the power.
  #   - saveName, if set, is the path to save the power analysis data.
  #   - recover indicates if data should be recovered (1) or not (0).
  
  # Select thinned steps in chain for posterior predictions:
  # We select nRep steps in the MCMC chain to generate simulated data for the 
  # power analysis. The steps are selected evenly from across the entire chain.
  chainLength = NROW( mcmcChain )
  stepIdxVec = seq( 1 , chainLength , floor(chainLength/nRep) )
  
  # For each selected step of the chain, we get the parameter values, create 
  # a simulated data set, run bayesian estimation to obtain the posterior 
  # distributions and then we check if the goals were achieved or not.
  nSim = 0
  if (recover == 1) {
    load( paste(saveName, "Power.Rdata", sep="") )
    nSim = nrow(goalTally)
  }
  while (nSim < length(stepIdxVec) ) {
    #for ( stepIdx in stepIdxVec ) {
    nSim = nSim + 1
    stepIdx = stepIdxVec[nSim]
    cat( "\n:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::\n" )
    cat( paste( "Power computation: Simulated Experiment" , nSim , "of" , 
                length(stepIdxVec) , ":\n\n" ) )
    
    # Get parameter values for this simulation:
    muVal = mcmcChain[stepIdx,"mu"]
    sigmaVal = mcmcChain[stepIdx,"sigma"]
    nuVal = mcmcChain[stepIdx,"nu"]
    # Generate simulated data:
    simulatedData = rt( N , df=nuVal ) * sigmaVal + muVal
    
    # Do bayesian analysis on simulated data:
    goalAchieved = goalAchievedForSample( simulatedData, muNULL, 
                                          effROPE, effHDImaxWid, mcmcLength )
    
    # Tally the results:
    # goalTally is a matrix to store the results of each simulated analysis.
    # It is created after the first iteration to consider the multiple goals
    # defined inside the goalAchievedForSample function.
    if (!exists("goalTally")) { # if goalTally does not exist, create it
      goalTally=matrix( nrow=0 , ncol=length(goalAchieved) ) 
    }
    goalTally = rbind( goalTally , goalAchieved )
    
    if ( !is.null(saveName) ) {  
      save( goalTally, file=paste(saveName, "Power.Rdata", sep="") )
    }
    
  }
  
  result_text = ""
  # Now we calculate the proportion that each goal was achieved:
  # For each goal...
  for ( goalIdx in 1:NCOL(goalTally) ) {
    # Extract the goal name for subsequent display:
    goalName = colnames(goalTally)[goalIdx]
    # Compute number of successes:
    goalHits = sum(unlist(goalTally[,goalIdx]))
    # Compute number of attempts:
    goalAttempts = NROW(goalTally)
    # Compute proportion of successes:
    goalEst = goalHits/goalAttempts
    # Compute HDI around proportion:
    goalEstHDI = HDIofICDF( qbeta ,
                            shape1=1+goalHits , 
                            shape2=1+goalAttempts-goalHits )
    # Display the result:
    powerResult = paste0( goalName,
                          ": Est.Power=" , round(goalEst,3) , 
                          "; Low Bound=" , round(goalEstHDI[1],3) ,
                          "; High Bound=" , round(goalEstHDI[2],3) )
    show( powerResult )
    result_text = paste(result_text, powerResult, "\n", sep="")
  }
  # Save the result:
  writeLines( result_text , con=paste(saveName, "PowerResult.txt", sep="") )
}