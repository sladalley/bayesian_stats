source("DBDA2E-utilities.R")
library('latex2exp')
normalise_data <- function(x, return_norm = TRUE, eps = 1e-6) {
  
  y_min <- min(x, na.rm = TRUE)
  y_max <- max(x, na.rm = TRUE)
  
  # min-max normalisation
  if (y_max == y_min) {
    x_norm <- rep(0.5, length(x))  # centre if constant signal
  } else {
    x_norm <- (x - y_min) / (y_max - y_min)
  }
  
  # push away from exact bounds [0,1] → [eps, 1-eps]
  x_norm <- x_norm * (1 - 2 * eps) + eps
  
  if (return_norm) {
    return(list(
      x_norm = x_norm,
      y_min = y_min,
      y_max = y_max
    ))
  } else {
    return(list(
      y_min = y_min,
      y_max = y_max
    ))
  }
}
# ====== MCMC chain ===========================================================
genMCMC_beta = function(datFrm, yName="y", gName="cond", muAlphaPrior = NULL, muBetaPrior = NULL, 
                   shapeKappa = NULL, rateKappa = NULL, numSavedSteps = 20000, saveName=NULL) { 
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
  if ( any( !is.finite(y) ) ) { stop("All y values alphast be finite.") }
  if (any(y <= 0)) {stop("Weibull requires strictly positive data!") }
  
  Ntotal = length(y)
  g = as.numeric(as.factor(datFrm[,gName])) # group data vector
  gLevels = levels(as.factor(datFrm[,gName]))
  nG = max(g) # number of groups
  


  y_norm = normalise_data(y)$x_norm
  
  y <- y_norm
  
  
  dataList = list(
     y = y,
    g = g,
    nG = nG,
    Ntotal = Ntotal
  )
  #-----------------------------------------------------------------------------
  # THE METRIC BETA MODEL.
  modelString = "
  model {
  
    for (i in 1:Ntotal) {
      y[i] ~ dbeta(alpha[g[i]], beta[g[i]])
    }
  
    for (gIdx in 1:nG) {
        mu[gIdx] ~ dbeta(1, 1)

      
       kappa[gIdx] ~ dgamma(2, 0.1)
     
  

      alpha[gIdx] <- mu[gIdx] * kappa[gIdx]
      beta[gIdx]  <- (1 - mu[gIdx]) * kappa[gIdx]
      
      #alpha[gIdx] ~ dunif(0,50)
      #beta[gIdx]  ~ dunif(0,50)
    }
  }
  "
  
  writeLines(modelString, con=paste0(saveName, "-model.txt"))
  
  # INTIALIZE THE CHAINS.
  # Initial values of MCMC chains based on data:
  #   The start of the MCMC chains is defined according to the data. 
  #   The estimated paremeters are: beta, alpha.
  #   The initial alpha and beta. are the ones obtained from the data and
  
  nChains = 10
  if (nG == 1) { 
    initsList = list(
      mu = mean(log(y)),
      sigma = sd(log(y))
    )
  }
  
  if (nG == 2) {
    baseInit = list( mu = c( mean(y[g==1]) , mean(y[g==2]) ) , 
                      kappa = c( 1 , 1) 
    )
  }
  
  initsList <- lapply(1:nChains,function(i){
    c(
      baseInit,
      list(
        .RNG.name = "base::Mersenne-Twister",
        .RNG.seed = 3000 + i
      )
      )
  })
  
  
  #-----------------------------------------------------------------------------
  # RUN THE CHAINS
  parameters = c("alpha", "beta") ## <- the parameters to be monitored
  thinSteps = 5
  adaptSteps = 500               # Number of steps to "tune" the samplers
  burnInSteps = 1000
  runjagsMethod=runjagsMethodDefault # from DBDA2E-utilities
 # nChains=nChainsDefault # from DBDA2E-utilities
  
  runJagsOut <- run.jags( method=runjagsMethod ,
                          model=paste0(saveName, "-model.txt") , 
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
    save( codaSamples , file=paste(saveName,"_", yName,"-Mcmc.Rdata",sep="") )
  }
  
  return(codaSamples)
}


# ====== Summary of data and diagnostics ======================================
smryMCMC_beta = function (codaSamples, nG, nullValEff=0, saveName=NULL, 
                           diagnostics=TRUE, graphFileType="png") {
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
  
  # ================= ONE GROUP =================
  if (nG == 1) {
    # alpha:
    summaryInfo = rbind(summaryInfo,
                        "alpha" = summarizePost(mcmcMat[,"alpha"])
    )
    
    summaryInfo = rbind(summaryInfo,
                        "beta" = summarizePost(mcmcMat[,"beta"])
    )
    
    
    
  }
  
  # ================= TWO GROUPS =================
  # Two groups:
  if (nG == 2) {
    # alpha parameter of group 1:
    summaryInfo = rbind( summaryInfo , 
                         "alpha[1]" = summarizePost( mcmcMat[,"alpha[1]"] , 
                                                  compVal=NULL , ROPE=NULL ) )
    # alpha parameter of group 1:
    summaryInfo = rbind( summaryInfo , 
                         "alpha[2]" = summarizePost( mcmcMat[,"alpha[2]"] , 
                                                  compVal=NULL , ROPE=NULL ) )
    
    # alpha parameter of group 1:
    summaryInfo = rbind( summaryInfo , 
                         "beta[1]" = summarizePost( mcmcMat[,"beta[1]"] , 
                                                     compVal=NULL , ROPE=NULL ) )
    
    # alpha parameter of group 1:
    summaryInfo = rbind( summaryInfo , 
                         "beta[2]" = summarizePost( mcmcMat[,"beta[2]"] , 
                                                     compVal=NULL , ROPE=NULL ) )
    
    # mean parameter of group 1:
    summaryInfo = rbind( summaryInfo , 
                         "mean[1]" = summarizePost( (mcmcMat[,"alpha[1]"]/(mcmcMat[,"alpha[1]"]+mcmcMat[,"beta[1]"])) , 
                                                    compVal=NULL , ROPE=NULL ) )
    
    # mean parameter of group 2:
    summaryInfo = rbind( summaryInfo , 
                         "mean[2]" = summarizePost( (mcmcMat[,"alpha[2]"]/(mcmcMat[,"alpha[2]"]+mcmcMat[,"beta[2]"])) , 
                                                    compVal=NULL , ROPE=NULL ) )
    
    # mean parameter of group 2:
    summaryInfo = rbind( summaryInfo , 
                         "kappa[1]" = summarizePost( ((mcmcMat[,"alpha[1]"]+mcmcMat[,"beta[1]"])) , 
                                                    compVal=NULL , ROPE=NULL ) )
    
    # mean parameter of group 2:
    summaryInfo = rbind( summaryInfo , 
                         "kappa[2]" = summarizePost( ((mcmcMat[,"alpha[2]"]+mcmcMat[,"beta[2]"])) , 
                                                    compVal=NULL , ROPE=NULL ) )
    
    # Difference between means, mu1 - mu2:
    summaryInfo = rbind( summaryInfo , 
                         "muDiff" = summarizePost( 
                           (mcmcMat[,"alpha[1]"]/(mcmcMat[,"alpha[1]"]+mcmcMat[,"beta[1]"])) - (mcmcMat[,"alpha[2]"]/(mcmcMat[,"alpha[2]"]+mcmcMat[,"beta[2]"])) , 
                           compVal=NULL , ROPE=NULL ) )
    
  }
  
  if ( !is.null(saveName) ) {
    write.csv( summaryInfo , file=paste0(saveName,"_", yName,"-SummaryInfo.csv") )
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
                saveName=paste0(diagnosticsPath, "/", sub(".*/", "", saveName,),"_", yName, "-") , 
                saveType=graphFileType )
    }
  }
}

# ====== Posterior distributions ==============================================
plotMCMC_beta = function( codaSamples, datFrm, yName="y", gName="cond", 
                             nullValEff=0, compValMu=NULL, compValMuDiff=NULL, 
                             compValSigma=NULL, compValSigmaDiff=NULL, compValNu=NULL, 
                             compValEff=NULL, ropeMu=NULL, ropeSigma = NULL, ropeEffSz=NULL, graphFileType="png", 
                             saveName=NULL, groupNames=c(1,2)) {
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
  #   - ropeMu: ROPE in the posterior mean (mu) plot.
  #   - ropeSigma: ROPE in the posterior s.d. (sigma) plot.
  #   - ropeEffSz: ROPE in the posterior effect size plot.
  #   - graphFileType: type of the image output files.
  #   - saveName: prefix of the output files.
  #   - groupNames: names of each group to be displayed in the plots.
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
  
  # Normalise Data
  y_norm <- numeric(length(y))
  eps <- 1e-6
  
  # Global min and max across ALL data
  y_min <- min(y)
  y_max <- max(y)
  
  # Apply scaling to all points
  y_scaled <- (y - y_min) / (y_max - y_min)
  
  # Push away from exact 0 and 1 for Beta
  y_norm <- (y_scaled * (1 - 2*eps)) + eps
  
  y <- y_norm
  
  cex_plotPost = 1.0
  # Posterior alpha and beta of one or two groups:
  if (nG == 1) {
    # One group:
    openGraph(height=2.5*nG,width=7)
    par( mar=c(3.5,3.5,3,1) , mgp=c(2.25,0.7,0) )
    layout( matrix(1:(2*nG),nrow=nG,byrow=TRUE) )
    alphaLim = range( mcmcMat[, grep("^alpha",colnames(mcmcMat)) ] )
    betaLim = range( mcmcMat[, grep("^beta",colnames(mcmcMat)) ] )
    plotPost( mcmcMat[,"alpha"] , xlab=TeX(sprintf("alpha")) , xlim=alphaLim , 
              compVal=compValMu , cex=cex_plotPost , col="#0088aa8a")
    plotPost( mcmcMat[,"beta"] , xlab=TeX(sprintf("beta")) , xlim=betaLim , 
              compVal=compValSigma , cex=cex_plotPost , col="#0088aa8a")
  }
  if (nG == 2) {
    # Two groups:
    openGraph(height=2.5*(nG+1),width=7)
    par( mar=c(3.5,3.5,3,1) , mgp=c(2.25,0.7,0) )
    layout( matrix(1:(2*(nG+1)),nrow=nG+1,byrow=TRUE) )
    alphaLim = range( mcmcMat[, grep("^alpha",colnames(mcmcMat)) ] )
    betaLim = range( mcmcMat[, grep("^beta",colnames(mcmcMat)) ] )
    
    plotPost( mcmcMat[,"alpha[1]"] , 
              xlab=TeX(sprintf("$\\alpha_{%s}$", groupNames[1])) , 
              xlim=alphaLim  , 
              col="#0088aa8a" , cex=cex_plotPost )
    plotPost( mcmcMat[,"beta[1]"] , 
              xlab=TeX(sprintf("$\\beta_{%s}$", groupNames[1])) , 
              xlim=betaLim , 
              col="#0088aa8a" , cex=cex_plotPost )
    plotPost( mcmcMat[,"alpha[2]"] , 
              xlab=TeX(sprintf("$\\alpha_{%s}$", groupNames[2])) , 
              xlim=alphaLim  , 
              col="#0088aa8a" , cex=cex_plotPost )
    plotPost( mcmcMat[,"beta[2]"] , 
              xlab=TeX(sprintf("$\\beta_{%s}$",  groupNames[2])) , 
              xlim=betaLim , 
              col="#0088aa8a" , cex=cex_plotPost )
  } 
  saveGraph( paste0(saveName,"_", yName,"-alphabeta") , type=graphFileType )
  
  # Posterior mu and kappa of one or two groups:
  if (nG == 1) {
    # One group:
    openGraph(height=2.5*nG,width=7)
    plotPost( mcmcMat[,"alpha"] , xlab=TeX(sprintf("alpha")) , xlim=alphaLim , 
              compVal=compValMu , cex=cex_plotPost , col="#0088aa8a")

  }
  if (nG == 2) {
    # Two groups:
    openGraph(height=2.5*(nG+1),width=7)
    par( mar=c(3.5,3.5,3,1) , mgp=c(2.25,0.7,0) )
    layout( matrix(1:(2*(nG+1)),nrow=nG+1,byrow=TRUE) )

    plotPost( (mcmcMat[,"alpha[1]"]/(mcmcMat[,"alpha[1]"]+mcmcMat[,"beta[1]"])) , 
              xlab=TeX(sprintf("$\\mu_{%s}$", groupNames[1])) , 
              col="#0088aa8a" , cex=cex_plotPost )
    plotPost( (mcmcMat[,"alpha[2]"]/(mcmcMat[,"alpha[2]"]+mcmcMat[,"beta[2]"])) , 
              xlab=TeX(sprintf("$\\mu_{%s}$", groupNames[2])) , 
              col="#0088aa8a" , cex=cex_plotPost )
    
    plotPost( ((mcmcMat[,"alpha[1]"]+mcmcMat[,"beta[1]"])) , 
              xlab=TeX(sprintf("$\\kappa_{%s}$", groupNames[1])) , 
              col="#0088aa8a" , cex=cex_plotPost )
    plotPost( ((mcmcMat[,"alpha[2]"]+mcmcMat[,"beta[2]"])) , 
              xlab=TeX(sprintf("$\\kappa_{%s}$", groupNames[2])) , 
              col="#0088aa8a" , cex=cex_plotPost )
    
    plotPost( (mcmcMat[,"alpha[1]"]/(mcmcMat[,"alpha[1]"]+mcmcMat[,"beta[1]"])) - (mcmcMat[,"alpha[2]"]/(mcmcMat[,"alpha[2]"]+mcmcMat[,"beta[2]"])) , compVal=compValMuDiff , ROPE=ropeMu,
              xlab=TeX(sprintf("$\\mu_{%s} - \\mu_{%s}$", groupNames[1],  groupNames[2])) , 
              cex=cex_plotPost , col="#0088aa8a" )

  } 
  saveGraph( paste0(saveName,"_", yName,"-mukappa") , type=graphFileType )
  
  
  
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
    yMax = 1.25 * max( histInfo$density )
    # Histogram:
    if (nG == 1) {
      histInfo = hist( thisY , prob=TRUE , ylim=c(0,yMax) , 
                       xlim=xLim , breaks=xBreaks , xlab="Data" , ylab="" , 
                       cex.lab=1.0 , col="#f4e3d7ff" , border="white" , yaxt="n" ,
                       main="")
    } else {
      histInfo = hist( thisY , prob=TRUE , ylim=c(0,yMax+0.2) , 
                       xlim=xLim , breaks=xBreaks , xlab="Data" , ylab=" " , 
                       cex.lab=1.0 , col="#f4e3d7ff" , border="white" ,
                       main=paste0("Group ", groupNames[gIdx]))
    }
    # Posterior predictive curves:
    stepIdxVec = seq( 1 , chainLength , floor(chainLength/nCurvesToPlot) )
    if (nG == 1) {
      for ( stepIdx in 1:length(stepIdxVec) ) {
        xVec = seq(eps, 1 - eps, length=501)
        alpha = mcmcMat[stepIdx, "alpha"]
        beta = mcmcMat[stepIdx, "beta"]

        
        dens = dbeta(xVec, alpha, beta)
        lines(xVec, dens, type="l" , col="#0088aa8a" , lwd=1 )
      }
    }
    if (nG > 1) {
      for ( stepIdx in 1:length(stepIdxVec) ) {
        xVec = seq(eps, 1 - eps, length=501)
        alpha = mcmcMat[stepIdxVec[stepIdx], paste0("alpha[", gIdx, "]")]
        beta = mcmcMat[stepIdxVec[stepIdx], paste0("beta[", gIdx, "]")]
        

        dens = dbeta(xVec, alpha, beta)
        lines(xVec, dens, type="l" , col="#0088aa8a" , lwd=1 )
      }
    }
    mtext( text=TeX(sprintf("$\\bar{x}=%f\\;\\,\\tilde{x}=%f\\;\\,s=%f$", 
                            round(mean(thisY), 2), 
                            round(median(thisY), 2), round(sd(thisY), 2))) ,
           at=ifelse(nG>=2,0.5*gIdx-0.2,0.5) , cex=1.0 , side=3 , outer=TRUE , 
           adj=c(0.5,0.5) , padj=c(-0.5,-0.5) )
  }
  saveGraph( paste0(saveName,"_", yName,"-PostPred") , type=graphFileType )
  
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
    
    pairs( cbind( mcmcMat[,"alpha"] , mcmcMat[,"beta"] )[plotIdx,] ,
           labels=c( TeX(sprintf("alpha")) , 
                     TeX(sprintf("beta"))),
           lower.panel=panel.cor , col="#0088aa8a" )
  }
  if (nG == 2) {
    openGraph(width=7,height=7)
    nPtToPlot = 1000
    plotIdx = floor(seq(1,length(mcmcMat[,"alpha[1]"]),
                        by=length(mcmcMat[,"alpha[1]"])/nPtToPlot))
    panel.cor = function(x, y, digits=2, prefix="", cex.cor, ...) {
      usr = par("usr"); on.exit(par(usr = usr))
      par(usr = c(0, 1, 0, 1))
      r = (cor(x, y))
      txt = format(c(r, 0.123456789), digits=digits)[1]
      txt = paste(prefix, txt, sep="")
      if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
      text(0.5, 0.5, txt, cex=1.25 ) # was cex=cex.cor*r
    }
    pairs( cbind( mcmcMat[,"alpha[1]"] , mcmcMat[,"alpha[2]"] , 
                  mcmcMat[,"beta[1]"], mcmcMat[,"beta[2]"])[plotIdx,] ,
           labels=c( TeX(sprintf("alpha %s", groupNames[1])) , 
                     TeX(sprintf("alpha %s", groupNames[2])) , 
                     TeX(sprintf("beta %s", groupNames[1])) , 
                     TeX(sprintf("beta %s", groupNames[2]))),
           lower.panel=panel.cor , col="#0088aa8a" )
    
  }
  saveGraph( file=paste(saveName,"_", yName,"-PostPairs",sep=""), type=graphFileType)
}

#===============================================================================
# Functions to estimate power in the robust bayesian analyses of metric data.
#   The functions are adapted from Jags-Ydich-Xnom1subj-MbernBeta-Power.R and
#   BEST.R (from *Bayesian estimation supersedes the t test* [Kruschke, 2013]) 
#   scripts.

#===============================================================================
# Functions to estimate power in the robust bayesian analyses of metric data.
#   The functions are adapted from Jags-Ydich-Xnom1subj-MbernBeta-Power.R and
#   BEST.R (from *Bayesian estimation supersedes the t test* [Kruschke, 2013]) 
#   scripts.

goalAchievedForSample = function( data, muNULL, effROPE, effHDImaxWid, muROPE,
                                  mcmcLength=10000, yName = 'y', gName = 'cond') {
  # Generate the MCMC chain:
  mcmcCoda = genMCMC_beta( datFrm = data , numSavedSteps=mcmcLength , saveName=NULL, yName = yName, gName = gName )
  mcmcMat = as.matrix(mcmcCoda,chains=TRUE)
  # Calculate effect size:
  HDImat = apply( mcmcMat , 2 , "HDIofMCMC" )
  mcmcMatDM =  (mcmcMat[,"alpha[1]"]/(mcmcMat[,"alpha[1]"]+mcmcMat[,"beta[1]"])) - (mcmcMat[,"alpha[2]"]/(mcmcMat[,"alpha[2]"]+mcmcMat[,"beta[2]"]))
  mcmcMatMu1 = (mcmcMat[,"alpha[1]"]/(mcmcMat[,"alpha[1]"]+mcmcMat[,"beta[1]"]))
  mcmcMatMu2 = (mcmcMat[,"alpha[2]"]/(mcmcMat[,"alpha[2]"]+mcmcMat[,"beta[2]"]))
  
  # Check goal achievement. First, compute the HDI of the effect size:
  effDM = HDIofMCMC( mcmcMatDM )
  effMu1 = HDIofMCMC(mcmcMatMu1)
  effMu2 = HDIofMCMC(mcmcMatMu2)
  
  
  # Define list for recording results:
  goalAchieved = list()
  
  effHDImaxWid = 0.5*(muROPE[2] - muROPE[1])
  effHDImaxWid = 0.02
  
  goalAchieved = c(goalAchieved,
                   "ExcludeROPE" = (effDM[1] > muROPE[2] | effDM[2] < muROPE[1])
  )
  
  goalAchieved = c(goalAchieved,
                   "InsideROPE" = (effDM[1] > muROPE[1] & effDM[2] < muROPE[2])
  )
  
  goalAchieved = c(goalAchieved,
                   "Decisive" = (
                     (effDM[1] > muROPE[2] | effDM[2] < muROPE[1]) |
                       (effDM[1] > muROPE[1] & effDM[2] < muROPE[2])
                   )
  )
  
  goalAchieved = c(goalAchieved,
                   "DMNarrowHDI" = (effDM[2] - effDM[1]) < effHDImaxWid
  )
  #print((effDM[2] - effDM[1]))
  #print(effHDImaxWid)
  goalAchieved = c(goalAchieved,
                   "PreciseDecision" =
    ((
      (effDM[1] > muROPE[2] | effDM[2] < muROPE[1]) |
        (effDM[1] > muROPE[1] & effDM[2] < muROPE[2])
    ) &
       (effDM[2] - effDM[1]) < effHDImaxWid)
  )
  
  
  
  
  
  # Return list of goal results:
  return(list(
    goalAchieved = goalAchieved,
    mcmcCoda = mcmcCoda
  ))
  
}

powerEstimation_beta = function( mcmcChain, N, muNULL, muROPE,
                                    effROPE, effHDImaxWid, 
                                    mcmcLength=10000, nRep=1000 , 
                                    saveName=NULL, recover=0, yName="y", gName="cond", groupNames) {
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
  N1 =N;
  N2 = N;
  nSim = 0
  splitted = strsplit(saveName, "/")
  savePath = ""
  
  for (i in 1:(length(splitted[[1]]) - 1)) {
    savePath = paste0(savePath, splitted[[1]][i], "/")
  }
  
  simPath = paste0(savePath, "mcmc_sims","_", yName)
  dir.create(simPath, showWarnings = FALSE, recursive = TRUE)
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
    #print(mcmcChain)
    #print(stepIdx)
    # Extract parameters for BOTH controllers
    alpha1 = mcmcChain[stepIdx,"alpha[1]"]
    alpha2 = mcmcChain[stepIdx,"alpha[2]"]
    
    beta1 = mcmcChain[stepIdx,"beta[1]"]
    beta2 = mcmcChain[stepIdx,"beta[2]"]
    
    
    # Simulate data for each controller
    y1 = rbeta(N1, alpha1, beta1)
    y2 = rbeta(N2, alpha2, beta2)
    
    
    
    # Compute trial-wise differences (paired case)
    # meanDiffSamples[i] <- mean(y1s) - mean(y2s)
    # Generate simulated data:
    dataMat <- matrix(0, ncol = 2, nrow = 0,
                      dimnames = list(NULL, c(yName,gName)))
    
    dataMat <- rbind(dataMat, cbind(y1, rep(groupNames[1], length(y1))))
    dataMat <- rbind(dataMat, cbind(y2, rep(groupNames[2], length(y2))))
    simulatedData = dataMat
    
    
    
    # print(simulatedData)
    # Do bayesian analysis on simulated data:
    goal = goalAchievedForSample( data.frame(simulatedData), muNULL, muROPE = muROPE,
                                          effROPE, effHDImaxWid, mcmcLength, yName = yName, gName = gName )
    goalAchieved = goal$goalAchieved
    mcmc = goal$mcmcCoda
  
    fileName = paste0(
      simPath, "/",
      sub(".*/", "", saveName), "_",
      yName, "_N", N,
      "_sim", nSim, ".Rdata"
    )
    save(mcmc, file = fileName)
    
    # Tally the results:
    # goalTally is a matrix to store the results of each simulated analysis.
    # It is created after the first iteration to consider the multiple goals
    # defined inside the goalAchievedForSample function.
    if (!exists("goalTally")) { # if goalTally does not exist, create it
      goalTally=matrix( nrow=0 , ncol=length(goalAchieved) ) 
    }
    goalTally = rbind( goalTally , goalAchieved )
    
    if ( !is.null(saveName) ) {  
      save( goalTally, file = paste0(saveName, "_N", N, "_", yName, "_Power.Rdata") )
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
  writeLines( result_text , con=paste0(saveName, "_N", N, "_", yName, "PowerResult.txt", sep="") )
}
ppc_beta = function(
    mcmcChain,
    N,
    mcmcLength = 10000,
    nRep = 100,
    saveName = NULL,
    recover = 0,
    yName = "y",
    gName = "cond",
    groupNames,
    realdat = NULL
){
  
  # --------------------------------------------------
  # basic checks
  # --------------------------------------------------
  if (is.null(realdat)) {
    stop("You must provide real data for PPC.")
  }
  
  if (!is.null(saveName)) {
    openGraph(
      height=2.5*(length(groupNames)+1),
      width  = 3.5 * length(groupNames)
    )
  }
  
  # --------------------------------------------------
  # select MCMC indices
  # --------------------------------------------------
  chainLength = NROW(mcmcChain)
  stepIdxVec = seq(1, chainLength, length.out = nRep)
  
  sim_data_list <- vector("list", length(stepIdxVec))
  
  # --------------------------------------------------
  # simulate datasets
  # --------------------------------------------------
  N1 = N[1]
  N2 = N[2]
  P = 100
  
  for (nSim in seq_along(stepIdxVec)) {
    cat( "\n:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::\n" )
    cat( paste( "PPC computation: Simulated Experiment" , nSim , "of" , 
                length(stepIdxVec) , ":\n\n" ) )
    
    stepIdx = stepIdxVec[nSim]
    
    alpha1 = mcmcChain[stepIdx, "alpha[1]"]
    alpha2 = mcmcChain[stepIdx, "alpha[2]"]
    
    beta1  = mcmcChain[stepIdx, "beta[1]"]
    beta2  = mcmcChain[stepIdx, "beta[2]"]
    
    y1 = rbeta(N1, alpha1, beta1)
    y2 = rbeta(N2, alpha2, beta2)
    
    sim_data_list[[nSim]] <- data.frame(
      value = c(y1, y2),
      group = factor(
        c(rep(groupNames[1], length(y1)),
          rep(groupNames[2], length(y2))),
        levels = groupNames
      ),
      sim_id = nSim
    )
  }
  
  sim_df <- do.call(rbind, sim_data_list)
  sim_keep <- unique(sim_df$sim_id)[seq(1, nRep, length.out = P)]
  sim_df <- sim_df[sim_df$sim_id %in% sim_keep, ]
  
 
  
  # --------------------------------------------------
  # NORMALISE REAL DATA
  # --------------------------------------------------
  res_real <- normalise_data(realdat[[yName]])
  
  obs_df <- data.frame(
    value = res_real$x_norm,
    group = factor(realdat[[gName]], levels = groupNames)
  )
  
  # --------------------------------------------------
  # SIMULATED DATA SPLIT BY GROUP
  # --------------------------------------------------
  sim_split <- split(sim_df, sim_df$group)
  obs_split <- split(obs_df, obs_df$group)
  
  # --------------------------------------------------
  # SET UP MULTI-PANEL PLOT (like your screenshot)
  # --------------------------------------------------
  par(mfrow = c(2, length(groupNames)))  # side-by-side plots
  par(mar = c(4, 4, 2, 1))
  
  for (g in groupNames) {
    
    thisY <- obs_split[[g]]$value
    
    # =========================
    # COMMON BINNING / LIMITS
    # =========================
    xLim <- c(
      min(thisY) - 0.1 * (max(thisY) - min(thisY)),
      max(thisY) + 0.1 * (max(thisY) - min(thisY))
    )
    
    xBreaks <- seq(
      xLim[1],
      xLim[2],
      length.out = ceiling((xLim[2] - xLim[1]) / (sd(thisY) / 4))
    )
    
    histInfo = hist(thisY, breaks=xBreaks, plot=FALSE)
    dens_max <- max(density(thisY)$y)
    hist_max <-max( histInfo$density )
    yMax = 1.1*max(dens_max,hist_max)
    
    # ==========================================================
    # ROW 1: DENSITY OVERLAYS
    # ==========================================================
    
    plot(
      density(thisY),
      col = "#f4e3d7ff",
      lwd = 2,
      xlim = xLim,
      ylim = c(0, yMax),
      main = paste("Group", g, "Density PPC"),
      xlab = "Value",
      ylab = "Density",
      bty = "n"
    )
    
    # simulated densities
    for (i in sim_keep) {
      
      sim_vals <- sim_split[[g]]$value[sim_split[[g]]$sim_id == i]
      
      d <- density(sim_vals)
      
      lines(
        d$x,
        d$y,
        col = rgb(0.2, 0.2, 0.2, 0.06),
        lwd = 1
      )
    }
    
    # ==========================================================
    # ROW 2: HISTOGRAM ENVELOPE
    # ==========================================================
    


    # plot observed histogram
    hist(
      thisY,
      probability = TRUE,
      breaks = xBreaks,
      xlim = xLim,
      ylim = c(0, yMax),
      col = "#f4e3d7ff",
      border = "white",
      main = paste("Group", g, "Histogram PPC"),
      xlab = "Value",
      ylab = "Density"
    )
    
    # simulated histograms
    for (i in sim_keep) {
      
      sim_vals <- sim_split[[g]]$value[sim_split[[g]]$sim_id == i]
      
      h_sim <- hist(
        sim_vals,
        breaks = xBreaks,
        probability = TRUE,
        plot = FALSE
      )
      
      rect(
        xleft   = h_sim$breaks[-length(h_sim$breaks)],
        ybottom = 0,
        xright  = h_sim$breaks[-1],
        ytop    = h_sim$density,
        #col = rgb(0.2, 0.2, 0.2, 0.06),
        border = rgb(0.2, 0.2, 0.2, 0.06)
      )

    }
  }
  saveGraph( file=paste(saveName,"_", yName,"-ppc",sep=""), type=graphFileType)
  return(list(
    sim_df = sim_df,
    obs_df = obs_df
  ))
}