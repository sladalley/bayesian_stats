source("DBDA2E-utilities.R")
library('latex2exp')

# plotPost summary for objectives
plotSumO = function( paramSampleVec , cenTend=c("mode","median","mean") , 
                     compVal=NULL, ROPE=NULL, credMass=0.95, HDItextPlace=0.7, 
                     xlab=NULL , xlim=NULL , yaxt=NULL , ylab=NULL , 
                     main=NULL , cex=NULL , cex.lab=NULL ,
                     col=NULL , border=NULL , showCurve=FALSE , breaks=NULL , 
                     ... ) {
  # Override defaults of hist function, if not specified by user:
  # (additional arguments "..." are passed to the hist function)
  if ( is.null(xlab) ) xlab="Param. Val."
  if ( is.null(cex.lab) ) cex.lab=1.5
  if ( is.null(cex) ) cex=1.4
  if ( is.null(xlim) ) xlim=range( c( compVal , ROPE , paramSampleVec ) )
  if ( is.null(main) ) main=""
  if ( is.null(yaxt) ) yaxt="n"
  if ( is.null(ylab) ) ylab=""
  if ( is.null(col) ) col="skyblue"
  if ( is.null(border) ) border="white"
  
  # convert coda object to matrix:
  if ( class(paramSampleVec) == "mcmc.list" ) {
    paramSampleVec = as.matrix(paramSampleVec)
  }
  
  summaryColNames = c("ESS","mean","median","mode",
                      "hdiMass","hdiLow","hdiHigh",
                      "compVal","pGtCompVal",
                      "ROPElow","ROPEhigh","pLtROPE","pInROPE","pGtROPE")
  postSummary = matrix( NA , nrow=1 , ncol=length(summaryColNames) , 
                        dimnames=list( c( xlab ) , summaryColNames ) )
  
  # require(coda) # for effectiveSize function
  postSummary[,"ESS"] = effectiveSize(paramSampleVec)
  
  postSummary[,"mean"] = mean(paramSampleVec)
  postSummary[,"median"] = median(paramSampleVec)
  mcmcDensity = density(paramSampleVec)
  postSummary[,"mode"] = mcmcDensity$x[which.max(mcmcDensity$y)]
  
  HDI = HDIofMCMC( paramSampleVec , credMass )
  postSummary[,"hdiMass"]=credMass
  postSummary[,"hdiLow"]=HDI[1]
  postSummary[,"hdiHigh"]=HDI[2]
  
  # Plot histogram.
  cvCol = "darkgreen"
  ropeCol = "darkred"
  if ( is.null(breaks) ) {
    if ( max(paramSampleVec) > min(paramSampleVec) ) {
      breaks = c( seq( from=min(paramSampleVec) , to=max(paramSampleVec) ,
                       by=(HDI[2]-HDI[1])/18 ) , max(paramSampleVec) )
    } else {
      breaks=c(min(paramSampleVec)-1.0E-6,max(paramSampleVec)+1.0E-6)
      border="skyblue"
    }
  }
  if ( !showCurve ) {
    par(xpd=NA)
    # histinfo = hist( paramSampleVec , xlab=xlab , yaxt=yaxt , ylab=ylab ,
    #                  freq=F , border=border , col=col ,
    #                  xlim=xlim , main=main , cex=cex , cex.lab=cex.lab ,
    #                  breaks=breaks , ... )
    histinfo = hist( paramSampleVec , plot=F )
    densCurve = density( paramSampleVec , adjust=2 )
    plot( densCurve$x , densCurve$y , type="n" , lwd=5 , col=col , bty="n" ,
          xlim=xlim , xlab=xlab , yaxt=yaxt , ylab=ylab , 
          main=main , cex=cex , cex.lab=cex.lab , ... )
  }
  if ( showCurve ) {
    par(xpd=NA)
    histinfo = hist( paramSampleVec , plot=F )
    densCurve = density( paramSampleVec , adjust=2 )
    plot( densCurve$x , densCurve$y , type="l" , lwd=5 , col=col , bty="n" ,
          xlim=xlim , xlab=xlab , yaxt=yaxt , ylab=ylab ,
          main=main , cex=cex , cex.lab=cex.lab , ... )
  }
  cenTendHt = 0.2*max(histinfo$density)
  cvHt = 0.7*max(histinfo$density)
  ROPEtextHt = 0.7*max(histinfo$density)
  
  # Display central tendency:
  mn = mean(paramSampleVec)
  med = median(paramSampleVec)
  mcmcDensity = density(paramSampleVec)
  mo = mcmcDensity$x[which.max(mcmcDensity$y)]
  # if ( "mode" %in% cenTend ){ 
  #   text( mo , cenTendHt ,
  #         bquote(mo==.(signif(mo,3))) , adj=c(-0.15,0) , cex=cex )
  # }
  if ( "median" %in% cenTend ){ 
    text( med , cenTendHt ,
          bquote(md==.(signif(med,3))) , adj=c(.5,0) , cex=cex , col=cvCol )
  }
  # if ( "mean" %in% cenTend ){ 
  #   text( mn , cenTendHt ,
  #         bquote(m==.(signif(mn,3))) , adj=c(1.15,0) , cex=cex )
  # }
  # Display the comparison value.
  if ( !is.null( compVal ) ) {
    pGtCompVal = sum( paramSampleVec > compVal ) / length( paramSampleVec ) 
    pLtCompVal = 1 - pGtCompVal
    lines( c(compVal,compVal) , c(0.96*cvHt,0) , 
           lty="dashed" , lwd=2.2 , col=cvCol )
    
    arrow_len = (HDI[2]-HDI[1])/4
    # below compVal:
    lines( c(compVal,compVal - arrow_len) , 
           c(0.96*cvHt,0.96*cvHt) , lwd=2 , 
           col=cvCol)
    text ( compVal - arrow_len , 1.2*cvHt , 
           bquote( .(round(100*pLtCompVal,1)) * "%") , col=cvCol)
    polygon( c(compVal - arrow_len, compVal - arrow_len, compVal - (arrow_len + (HDI[2]-HDI[1])/18)), 
             c(0.91*cvHt, 1.01*cvHt, 0.96*cvHt), 
             col=cvCol, border=cvCol)
    
    # above compVal:
    lines( c(compVal,compVal + arrow_len) , 
           c(0.96*cvHt,0.96*cvHt) , lwd=2 , 
           col=cvCol)
    text ( compVal + arrow_len , 1.2*cvHt , 
           bquote( .(round(100*pGtCompVal,1)) * "%") , col=cvCol)
    polygon( c(compVal + arrow_len, compVal + arrow_len, compVal + (arrow_len + (HDI[2]-HDI[1])/18)), 
             c(0.91*cvHt, 1.01*cvHt, 0.96*cvHt), 
             col=cvCol, border=cvCol)
    
    #text( compVal , cvHt ,
    #      bquote( .(round(100*pLtCompVal,1)) * "% < " *
    #               .(signif(compVal,3)) * " < " * 
    #               .(round(100*pGtCompVal,1)) * "%" ) ,
    #      adj=c(pLtCompVal,0) , cex=0.8*cex , col=cvCol )
    
    postSummary[,"compVal"] = compVal
    postSummary[,"pGtCompVal"] = pGtCompVal
  }
  # Display the ROPE.
  if ( !is.null( ROPE ) ) {
    pInROPE = ( sum( paramSampleVec > ROPE[1] & paramSampleVec < ROPE[2] )
                / length( paramSampleVec ) )
    pGtROPE = ( sum( paramSampleVec >= ROPE[2] ) / length( paramSampleVec ) )
    pLtROPE = ( sum( paramSampleVec <= ROPE[1] ) / length( paramSampleVec ) )
    lines( c(ROPE[1],ROPE[1]) , c(0.96*ROPEtextHt,0) , lty="dashed" , lwd=2.2 ,
           col=ropeCol )
    lines( c(ROPE[2],ROPE[2]) , c(0.96*ROPEtextHt,0) , lty="dashed" , lwd=2.2 ,
           col=ropeCol)
    
    arrow_len = (HDI[2]-HDI[1])/4
    # between rope:
    lines( c(ROPE[1],ROPE[2]) , c(1.1*ROPEtextHt,1.1*ROPEtextHt) , lwd=2 , 
           col=ropeCol) 
    lines( c(ROPE[1],ROPE[1]) , c(1.05*ROPEtextHt,1.15*ROPEtextHt) , lwd=2 , 
           col=ropeCol)
    lines( c(ROPE[2],ROPE[2]) , c(1.05*ROPEtextHt,1.15*ROPEtextHt) , lwd=2 , 
           col=ropeCol)
    text( mean(ROPE) , 1.35*ROPEtextHt , 
          bquote( .(round(100*pInROPE,1)) * "%") , col=ropeCol)
    
    # below rope:
    lines( c(ROPE[1],ROPE[1] - arrow_len) , 
           c(0.96*ROPEtextHt,0.96*ROPEtextHt) , lwd=2 , 
           col=ropeCol)
    text ( ROPE[1] - arrow_len , 1.2*ROPEtextHt , 
           bquote( .(round(100*pLtROPE,1)) * "%") , col=ropeCol)
    polygon( c(ROPE[1] - arrow_len, ROPE[1] - arrow_len, ROPE[1] - (arrow_len + (HDI[2]-HDI[1])/18)), 
             c(0.91*ROPEtextHt, 1.01*ROPEtextHt, 0.96*ROPEtextHt), 
             col=ropeCol, border=ropeCol)
    
    # above rope:
    lines( c(ROPE[2],ROPE[2] + arrow_len) , 
           c(0.96*ROPEtextHt,0.96*ROPEtextHt) , lwd=2 , 
           col=ropeCol)
    text ( ROPE[2] + arrow_len , 1.2*ROPEtextHt , 
           bquote( .(round(100*pGtROPE,1)) * "%") , col=ropeCol)
    polygon( c(ROPE[2] + arrow_len, ROPE[2] + arrow_len, ROPE[2] + (arrow_len + (HDI[2]-HDI[1])/18)), 
             c(0.91*ROPEtextHt, 1.01*ROPEtextHt, 0.96*ROPEtextHt), 
             col=ropeCol, border=ropeCol)
    
    #text( mean(ROPE) , ROPEtextHt ,
    #     bquote( .(round(100*pLtROPE,1)) * "% < " * .(ROPE[1]) * " < " * 
    #               .(round(100*pInROPE,1)) * "% < " * .(ROPE[2]) * " < " * 
    #               .(round(100*pGtROPE,1)) * "%" ) ,
    #      adj=c(pLtROPE+.5*pInROPE,0) , cex=1 , col=ropeCol )
    
    postSummary[,"ROPElow"]=ROPE[1] 
    postSummary[,"ROPEhigh"]=ROPE[2] 
    postSummary[,"pLtROPE"]=pLtROPE
    postSummary[,"pInROPE"]=pInROPE
    postSummary[,"pGtROPE"]=pGtROPE
  }
  # Display the HDI.
  lines( HDI , c(0,0) , lwd=4 , lend=1 )
  # text( mean(HDI) , 0 , bquote(.(100*credMass) * "% HDI" ) ,
  #       adj=c(.5,-1.7) , cex=cex )
  #text( mean(HDI) , 0 , bquote(.(100*credMass) * "% HDI" ) ,
  #      adj=c(.5,-0.8) , cex=cex )
  text( HDI[1] , 0 , bquote(.(signif(HDI[1],3))) ,
        adj=c(HDItextPlace,-0.5) , cex=cex )
  text( HDI[2] , 0 , bquote(.(signif(HDI[2],3))) ,
        adj=c(1.0-HDItextPlace,-0.5) , cex=cex )
  par(xpd=F)
  #
  return( postSummary )
}

# plotPost summary for subjectives
plotSumS = function( paramSampleVec , cenTend=c("mode","median","mean") , 
                     compVal=NULL, ROPE=NULL, credMass=0.95, HDItextPlace=0.7, 
                     xlab=NULL , xlim=NULL , yaxt=NULL , ylab=NULL , 
                     main=NULL , cex=NULL , cex.lab=NULL ,
                     col=NULL , border=NULL , showCurve=FALSE , breaks=NULL , 
                     ... ) {
  # Override defaults of hist function, if not specified by user:
  # (additional arguments "..." are passed to the hist function)
  if ( is.null(xlab) ) xlab="Param. Val."
  if ( is.null(cex.lab) ) cex.lab=1.5
  if ( is.null(cex) ) cex=1.4
  if ( is.null(xlim) ) xlim=range( c( compVal , ROPE , paramSampleVec ) )
  if ( is.null(main) ) main=""
  if ( is.null(yaxt) ) yaxt="n"
  if ( is.null(ylab) ) ylab=""
  if ( is.null(col) ) col="skyblue"
  if ( is.null(border) ) border="white"
  
  # convert coda object to matrix:
  if ( class(paramSampleVec) == "mcmc.list" ) {
    paramSampleVec = as.matrix(paramSampleVec)
  }
  
  summaryColNames = c("ESS","mean","median","mode",
                      "hdiMass","hdiLow","hdiHigh",
                      "compVal","pGtCompVal",
                      "ROPElow","ROPEhigh","pLtROPE","pInROPE","pGtROPE")
  postSummary = matrix( NA , nrow=1 , ncol=length(summaryColNames) , 
                        dimnames=list( c( xlab ) , summaryColNames ) )
  
  # require(coda) # for effectiveSize function
  postSummary[,"ESS"] = effectiveSize(paramSampleVec)
  
  postSummary[,"mean"] = mean(paramSampleVec)
  postSummary[,"median"] = median(paramSampleVec)
  mcmcDensity = density(paramSampleVec)
  postSummary[,"mode"] = mcmcDensity$x[which.max(mcmcDensity$y)]
  
  HDI = HDIofMCMC( paramSampleVec , credMass )
  postSummary[,"hdiMass"]=credMass
  postSummary[,"hdiLow"]=HDI[1]
  postSummary[,"hdiHigh"]=HDI[2]
  
  # Plot histogram.
  cvCol = "darkgreen"
  ropeCol = "darkred"
  if ( is.null(breaks) ) {
    if ( max(paramSampleVec) > min(paramSampleVec) ) {
      breaks = c( seq( from=min(paramSampleVec) , to=max(paramSampleVec) ,
                       by=(HDI[2]-HDI[1])/18 ) , max(paramSampleVec) )
    } else {
      breaks=c(min(paramSampleVec)-1.0E-6,max(paramSampleVec)+1.0E-6)
      border="skyblue"
    }
  }
  if ( !showCurve ) {
    par(xpd=NA)
    # histinfo = hist( paramSampleVec , xlab=xlab , yaxt=yaxt , ylab=ylab ,
    #                  freq=F , border=border , col=col ,
    #                  xlim=xlim , main=main , cex=cex , cex.lab=cex.lab ,
    #                  breaks=breaks , ... )
    histinfo = hist( paramSampleVec , plot=F )
    densCurve = density( paramSampleVec , adjust=2 )
    plot( densCurve$x , densCurve$y , type="n" , lwd=5 , col=col , bty="n" ,
          xlim=xlim , xlab=xlab , yaxt=yaxt , ylab=ylab , 
          main=main , cex=cex , cex.lab=cex.lab , ... )
  }
  if ( showCurve ) {
    par(xpd=NA)
    histinfo = hist( paramSampleVec , plot=F )
    densCurve = density( paramSampleVec , adjust=2 )
    plot( densCurve$x , densCurve$y , type="l" , lwd=5 , col=col , bty="n" ,
          xlim=xlim , xlab=xlab , yaxt=yaxt , ylab=ylab ,
          main=main , cex=cex , cex.lab=cex.lab , ... )
  }
  cenTendHt = 0.2*max(histinfo$density)
  cvHt = 0.7*max(histinfo$density)
  ROPEtextHt = 0.7*max(histinfo$density)
  
  # Display central tendency:
  mn = mean(paramSampleVec)
  med = median(paramSampleVec)
  mcmcDensity = density(paramSampleVec)
  mo = mcmcDensity$x[which.max(mcmcDensity$y)]
  # if ( "mode" %in% cenTend ){ 
  #   text( mo , cenTendHt ,
  #         bquote(mo==.(signif(mo,3))) , adj=c(-0.15,0) , cex=cex )
  # }
  if ( "median" %in% cenTend ){ 
    text( med , cenTendHt ,
          bquote(md==.(signif(med,3))) , adj=c(.5,-0.4) , cex=cex , col=cvCol )
  }
  # if ( "mean" %in% cenTend ){ 
  #   text( mn , cenTendHt ,
  #         bquote(m==.(signif(mn,3))) , adj=c(1.15,0) , cex=cex )
  # }
  # Display the comparison value.
  if ( !is.null( compVal ) ) {
    pGtCompVal = sum( paramSampleVec > compVal ) / length( paramSampleVec ) 
    pLtCompVal = 1 - pGtCompVal
    lines( c(compVal,compVal) , c(1.16*cvHt,0) , 
           lty="dashed" , lwd=2.2 , col=cvCol )
    
    arrow_len = (HDI[2]-HDI[1])/4
    # below compVal:
    lines( c(compVal,compVal - arrow_len) , 
           c(1.16*cvHt,1.16*cvHt) , lwd=2 , 
           col=cvCol)
    text ( compVal - arrow_len , 1.6*cvHt , 
           bquote( .(round(100*pLtCompVal,1)) * "%") , col=cvCol)
    polygon( c(compVal - arrow_len, compVal - arrow_len, compVal - (arrow_len + (HDI[2]-HDI[1])/18)), 
             c(1.01*cvHt, 1.31*cvHt, 1.16*cvHt), 
             col=cvCol, border=cvCol)
    
    # above compVal:
    lines( c(compVal,compVal + arrow_len) , 
           c(1.16*cvHt,1.16*cvHt) , lwd=2 , 
           col=cvCol)
    text ( compVal + arrow_len , 1.6*cvHt , 
           bquote( .(round(100*pGtCompVal,1)) * "%") , col=cvCol)
    polygon( c(compVal + arrow_len, compVal + arrow_len, compVal + (arrow_len + (HDI[2]-HDI[1])/18)), 
             c(1.01*cvHt, 1.31*cvHt, 1.16*cvHt), 
             col=cvCol, border=cvCol)
    
    #text( compVal , cvHt ,
    #      bquote( .(round(100*pLtCompVal,1)) * "% < " *
    #               .(signif(compVal,3)) * " < " * 
    #               .(round(100*pGtCompVal,1)) * "%" ) ,
    #      adj=c(pLtCompVal,0) , cex=0.8*cex , col=cvCol )
    
    postSummary[,"compVal"] = compVal
    postSummary[,"pGtCompVal"] = pGtCompVal
  }
  # Display the ROPE.
  if ( !is.null( ROPE ) ) {
    pInROPE = ( sum( paramSampleVec > ROPE[1] & paramSampleVec < ROPE[2] )
                / length( paramSampleVec ) )
    pGtROPE = ( sum( paramSampleVec >= ROPE[2] ) / length( paramSampleVec ) )
    pLtROPE = ( sum( paramSampleVec <= ROPE[1] ) / length( paramSampleVec ) )
    lines( c(ROPE[1],ROPE[1]) , c(0.96*ROPEtextHt,0) , lty="dashed" , lwd=2.2 ,
           col=ropeCol )
    lines( c(ROPE[2],ROPE[2]) , c(0.96*ROPEtextHt,0) , lty="dashed" , lwd=2.2 ,
           col=ropeCol)
    
    arrow_len = (HDI[2]-HDI[1])/4
    # between rope:
    lines( c(ROPE[1],ROPE[2]) , c(1.1*ROPEtextHt,1.1*ROPEtextHt) , lwd=2 , 
           col=ropeCol) 
    lines( c(ROPE[1],ROPE[1]) , c(1.05*ROPEtextHt,1.15*ROPEtextHt) , lwd=2 , 
           col=ropeCol)
    lines( c(ROPE[2],ROPE[2]) , c(1.05*ROPEtextHt,1.15*ROPEtextHt) , lwd=2 , 
           col=ropeCol)
    text( mean(ROPE) , 1.23*ROPEtextHt , 
          bquote( .(round(100*pInROPE,1)) * "%") , col=ropeCol)
    
    # below rope:
    lines( c(ROPE[1],ROPE[1] - arrow_len) , 
           c(0.96*ROPEtextHt,0.96*ROPEtextHt) , lwd=2 , 
           col=ropeCol)
    text ( ROPE[1] - arrow_len , 1.08*ROPEtextHt , 
           bquote( .(round(100*pLtROPE,1)) * "%") , col=ropeCol)
    polygon( c(ROPE[1] - arrow_len, ROPE[1] - arrow_len, ROPE[1] - (arrow_len + (HDI[2]-HDI[1])/18)), 
             c(0.91*ROPEtextHt, 1.01*ROPEtextHt, 0.96*ROPEtextHt), 
             col=ropeCol, border=ropeCol)
    
    # above rope:
    lines( c(ROPE[2],ROPE[2] + arrow_len) , 
           c(0.96*ROPEtextHt,0.96*ROPEtextHt) , lwd=2 , 
           col=ropeCol)
    text ( ROPE[2] + arrow_len , 1.08*ROPEtextHt , 
           bquote( .(round(100*pGtROPE,1)) * "%") , col=ropeCol)
    polygon( c(ROPE[2] + arrow_len, ROPE[2] + arrow_len, ROPE[2] + (arrow_len + (HDI[2]-HDI[1])/18)), 
             c(0.91*ROPEtextHt, 1.01*ROPEtextHt, 0.96*ROPEtextHt), 
             col=ropeCol, border=ropeCol)
    
    #text( mean(ROPE) , ROPEtextHt ,
    #     bquote( .(round(100*pLtROPE,1)) * "% < " * .(ROPE[1]) * " < " * 
    #               .(round(100*pInROPE,1)) * "% < " * .(ROPE[2]) * " < " * 
    #               .(round(100*pGtROPE,1)) * "%" ) ,
    #      adj=c(pLtROPE+.5*pInROPE,0) , cex=1 , col=ropeCol )
    
    postSummary[,"ROPElow"]=ROPE[1] 
    postSummary[,"ROPEhigh"]=ROPE[2] 
    postSummary[,"pLtROPE"]=pLtROPE
    postSummary[,"pInROPE"]=pInROPE
    postSummary[,"pGtROPE"]=pGtROPE
  }
  # Display the HDI.
  lines( HDI , c(0,0) , lwd=4 , lend=1 )
  # text( mean(HDI) , 0 , bquote(.(100*credMass) * "% HDI" ) ,
  #       adj=c(.5,-1.7) , cex=cex )
  #text( mean(HDI) , 0 , bquote(.(100*credMass) * "% HDI" ) ,
  #      adj=c(.5,-0.8) , cex=cex )
  text( HDI[1] , 0 , bquote(.(signif(HDI[1],3))) ,
        adj=c(HDItextPlace,-0.5) , cex=cex )
  text( HDI[2] , 0 , bquote(.(signif(HDI[2],3))) ,
        adj=c(1.0-HDItextPlace,-0.5) , cex=cex )
  par(xpd=F)
  #
  return( postSummary )
}

plotObjectives = function( codaSamples, datFrm, yName="y", gName="cond", 
                           nullHypValMu=0, compValMu=NULL, compValMuDiff=NULL, 
                           compValSigma=NULL, compValSigmaDiff=NULL, compValNu=NULL, 
                           compValEff=NULL, ropeEffSz=NULL, graphFileType="png", 
                           saveName=NULL, groupNames=c(1,2), subscript="", subsEffsz="", 
                           plotOption=2, summ=FALSE) {
  # Display posterior information.
  # List of parameters:
  #   - codaSamples: codaSamples object with the MCMC chain.
  #   - datFrm: data set.
  #   - yName: name of the column with the data.
  #   - gName: name of the column with the group indexes. 
  #   - nullHypValMu: null hypothesis value for the mean (for single group).
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
  #   - plotOption: plots in one line (1) or two lines (2)
  #   - summary: if TRUE, plot only the summary of distributions
  mcmcMat = as.matrix(codaSamples,chains=TRUE)
  chainLength = NROW( mcmcMat )
  
  y = as.numeric(datFrm[,yName])
  Ntotal = length(y)
  g = as.numeric(as.factor(datFrm[,gName])) # group data vector
  gLevels = levels(as.factor(datFrm[,gName]))
  nG = max(g) # number of groups
  
  cex_plotPost = 1.0
  
  if (plotOption == 1 & summ == FALSE) {
    openGraph(height=2.5,width=13)
    par( mar=c(3.5,0,3,1) , mgp=c(2.25,0.7,0))
    layout( matrix(c(1,2,3,4,5), nrow=1, byrow=TRUE) )
  }
  if (plotOption == 1 & summ == TRUE) {
    openGraph(height=1.5,width=13)
    par( mar=c(3.5,0,3,1) , mgp=c(2.25,0.7,0))
    layout( matrix(c(1,2,3,4,5), nrow=1, byrow=TRUE) )
  }
  if (plotOption == 2) {
    openGraph(height=2.5*2,width=9)
    par( mar=c(3.5,1,3,1) , mgp=c(2.25,0.7,0) )
    layout( matrix(c(1,1,2,2,3,3,
                     0,4,4,5,5,0), nrow=2, byrow=TRUE) )
  }

  # Posterior mu, sigma, and nu of one group:
  muLim = range( mcmcMat[, grep("^mu",colnames(mcmcMat)) ] )
  sigmaLim = range( mcmcMat[, grep("^sigma",colnames(mcmcMat)) ] )
  if (summ) {
    plotSumO( mcmcMat[,"mu"] , xlab=TeX(sprintf("$\\mu_{%s}$", subscript)) , xlim=muLim , 
              compVal=compValMu , cex=cex_plotPost , col="#0088aa8a",
              main="Mean")
    plotSumO( mcmcMat[,"sigma"] , xlab=TeX(sprintf("$\\tau_{%s}$", subscript)) , xlim=sigmaLim , 
              compVal=compValSigma , cex=cex_plotPost , col="#0088aa8a", 
              main="Scale")
    plotSumO( log10(mcmcMat[,"nu"]) , 
              xlab=TeX(sprintf("$\\log_{10}\\nu_{%s}$", subscript)) , 
              compVal=compValNu , ROPE=NULL , cex=cex_plotPost , col="#0088aa8a",
              main="Normality")
  } else {
    plotPost( mcmcMat[,"mu"] , xlab=TeX(sprintf("$\\mu_{%s}$", subscript)) , xlim=muLim , 
              compVal=compValMu , cex=cex_plotPost , col="#0088aa8a",
              main="Mean")
    plotPost( mcmcMat[,"sigma"] , xlab=TeX(sprintf("$\\tau_{%s}$", subscript)) , xlim=sigmaLim , 
              compVal=compValSigma , cex=cex_plotPost , col="#0088aa8a", 
              main="Scale")
    plotPost( log10(mcmcMat[,"nu"]) , 
              xlab=TeX(sprintf("$\\log_{10}\\nu_{%s}$", subscript)) , 
              compVal=compValNu , ROPE=NULL , cex=cex_plotPost , col="#0088aa8a",
              main="Normality")
  }
  
  # Posterior effect size:
  postEffSz = ( mcmcMat[,"mu"] - nullHypValMu ) / mcmcMat[,"sigma"]
  if (summ) {
    plotSumO( postEffSz , xlab=TeX(sprintf("$d_{%s,%s}$", subsEffsz, subscript)) , 
              compVal=compValEff , ROPE=ropeEffSz , col="#0088aa8a" , 
              cex=cex_plotPost, main="Effect size" )
  } else {
    plotPost( postEffSz , xlab=TeX(sprintf("$d_{%s,%s}$", subsEffsz, subscript)) , 
              compVal=compValEff , ROPE=ropeEffSz , col="#0088aa8a" , 
              cex=cex_plotPost, main="Effect size" )
  }
  
  # Posterior predictive:
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
    histInfo = hist( thisY , prob=TRUE , ylim=c(0,yMax) , 
                     xlim=xLim , breaks=xBreaks , xlab="Data" , ylab="" , 
                     cex.lab=1.0 , col="#f4e3d7ff" , border="white" , yaxt="n" ,
                     main="Histogram and credible distributions")
    text(par('usr')[1]+(par('usr')[2] - par('usr')[1])/2,0.9*par('usr')[4],
         labels=TeX(sprintf("$\\bar{x}=%f\\;\\,\\tilde{x}=%f\\;\\,s=%f$", round(mean(y), 2), 
                            round(median(y), 2), round(sd(y), 2))))
    
    # Posterior predictive curves:
    stepIdxVec = seq( 1 , chainLength , floor(chainLength/nCurvesToPlot) )
    for ( stepIdx in 1:length(stepIdxVec) ) {
      xVec = seq( xLim[1] , xLim[2] , length=501 )
      lines(xVec, dt( (xVec-mcmcMat[ stepIdxVec[stepIdx] , "mu" ]) /
                        mcmcMat[ stepIdxVec[stepIdx] , "sigma" ], 
                      df=mcmcMat[ stepIdxVec[stepIdx] , "nu" ] ) /
              mcmcMat[ stepIdxVec[stepIdx] , "sigma" ] , 
            type="l" , col="#0088aa8a" , lwd=1 )
    }
  }
  
  if (summ) 
    saveGraph( file=paste0(saveName,"-PlotsSumm",plotOption), type=graphFileType)
  else
    saveGraph( file=paste0(saveName,"-Plots",plotOption), type=graphFileType)
}

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
          data = rbind(data, c(level, item, subj, gIdx))
          extra[item, level] = extra[item, level] - 1
          break
        }
      }
    }
  }
  
  return(list("data"=data, "extra"=extra_return))
}

plotSubjectives = function (ordCodaSamples, datFrm, yName="resp", qName="item", 
                            gName="cond", nullHypValMu=0, compValMu=NULL, 
                            compValMuDiff=NULL, compValSigma=NULL, compValSigmaDiff=NULL, 
                            compValNu=NULL, compValEff=NULL, ropeEffSz=NULL, 
                            minLevel, maxLevel, graphFileType="png", saveName=NULL, 
                            groupNames=c(1,2), subscript="", subsEffsz="",
                            extraInfo=NULL, summ=FALSE) {
  # Display posterior information.
  # List of parameters:
  #   - ordCodaSamples: codaSamples object with the MCMC chain.
  #   - datFrm: data set.
  #   - yName: name of the column with the answers (ordinal values).
  #   - qName: name of the column with the item indexes.
  #   - gName: name of the column with the group indexes. 
  #   - nullHypValMu: null hypothesis value for the mean (for single group).
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
  #   - summary: if TRUE, plot only the summary of distributions
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
  
  if (summ) {
    openGraph(height=2.3,width=12)
    par( mar=c(3.5,1,3,1) , mgp=c(2.25,0.7,0) )
    layout( matrix(c(1,3,5,7,
                     2,4,6,7), nrow=2, byrow=TRUE) )
  } else {
    openGraph(height=2.5*(nG+1),width=7)
    par( mar=c(3.5,3.5,3,1) , mgp=c(2.25,0.7,0) )
    layout( matrix(1:9,nrow=3,byrow=TRUE) )
  }
  
  if (summ) {
    # Posterior mu, sigma, nu, and differences:
    muLim = range( ordMcmcMat[, grep("^mu\\[",colnames(ordMcmcMat)) ] )
    sigmaLim = range( ordMcmcMat[, grep("^sigma\\[",colnames(ordMcmcMat)) ] )
    plotSumS( ordMcmcMat[,"mu[1]"] , 
              xlab=TeX(sprintf("$\\mu_{%s,%s}$", subscript, groupNames[1])) , 
              xlim=muLim , compVal=compValMu[1] , 
              col="#0088aa8a" , cex=cex_plotPost , 
              main=paste0("Mean ",groupNames[1]))
    plotSumS( ordMcmcMat[,"sigma[1]"] , 
              xlab=TeX(sprintf("$\\tau_{%s,%s}$", subscript, groupNames[1])) , 
              xlim=sigmaLim , compVal=compValSigma[1] , 
              col="#0088aa8a" , cex=cex_plotPost , 
              main=paste0("Scale ",groupNames[1]))
    plotSumS( ordMcmcMat[,"mu[2]"] , 
              xlab=TeX(sprintf("$\\mu_{%s,%s}$", subscript, groupNames[2])) , 
              xlim=muLim , compVal=compValMu[2] , 
              col="#0088aa8a" , cex=cex_plotPost,  
              main=paste0("Mean ",groupNames[2])) 
    plotSumS( ordMcmcMat[,"sigma[2]"] , 
              xlab=TeX(sprintf("$\\tau_{%s,%s}$", subscript, groupNames[2])) , 
              xlim=sigmaLim , compVal=compValSigma[2] , 
              col="#0088aa8a" , cex=cex_plotPost , 
              main=paste0("Scale ",groupNames[2]))
    plotSumS( ordMcmcMat[,"mu[2]"]-ordMcmcMat[,"mu[1]"] , compVal=compValMuDiff ,
              xlab=TeX(sprintf("$\\mu_{%s,%s} - \\mu_{%s,%s}$", subscript, groupNames[2], subscript, groupNames[1])) , 
              cex=cex_plotPost , col="#0088aa8a" , 
              main="Difference between means")
    plotSumS( ordMcmcMat[,"sigma[2]"]-ordMcmcMat[,"sigma[1]"] , compVal=compValSigmaDiff , 
              xlab=TeX(sprintf("$\\tau_{%s,%s} - \\tau_{%s,%s}$", subscript, groupNames[2], subscript, groupNames[1])) , 
              cex=cex_plotPost , col="#0088aa8a" , 
              main="Difference between scales")
    
    # Posterior effect size:
    postEffSz = ( ( ordMcmcMat[,"mu[2]"] - ordMcmcMat[,"mu[1]"] ) 
                  / sqrt( ( ordMcmcMat[,"sigma[1]"]^2 + ordMcmcMat[,"sigma[2]"]^2 ) / 2 ) )
      
    plotSumS( postEffSz , compVal=compValEff , ROPE=ropeEffSz , 
              col="#0088aa8a" , cex=cex_plotPost , cex.lab=1.2,
              xlab=TeX(sprintf("$d_{%s,%s}$", subsEffsz, subscript)) ,
              main="Effect size")
    
    saveGraph( paste0(saveName,"-PlotsSumm") , type=graphFileType )
  } else {
    # Posterior mu, sigma, nu, and differences:
    muLim = range( ordMcmcMat[, grep("^mu\\[",colnames(ordMcmcMat)) ] )
    sigmaLim = range( ordMcmcMat[, grep("^sigma\\[",colnames(ordMcmcMat)) ] )
    plotPost( ordMcmcMat[,"mu[1]"] , 
              xlab=TeX(sprintf("$\\mu_{%s,%s}$", subscript, groupNames[1])) , 
              xlim=muLim , compVal=compValMu[1] , 
              col="#0088aa8a" , cex=cex_plotPost , 
              main=paste0("Mean ",groupNames[1]))
    plotPost( ordMcmcMat[,"sigma[1]"] , 
              xlab=TeX(sprintf("$\\tau_{%s,%s}$", subscript, groupNames[1])) , 
              xlim=sigmaLim , compVal=compValSigma[1] , 
              col="#0088aa8a" , cex=cex_plotPost , 
              main=paste0("Scale ",groupNames[1]))
    plotPost( log10(ordMcmcMat[,"nu[1]"]) , 
              xlab=TeX(sprintf("$\\log_{10}\\nu_{%s,%s}$", subscript, groupNames[1])) , 
              compVal=compValNu[1] , ROPE=NULL , cex=cex_plotPost , col="#0088aa8a" , 
              main=paste0("Normality ",groupNames[1]))
    plotPost( ordMcmcMat[,"mu[2]"] , 
              xlab=TeX(sprintf("$\\mu_{%s,%s}$", subscript, groupNames[2])) , 
              xlim=muLim , compVal=compValMu[2] , 
              col="#0088aa8a" , cex=cex_plotPost,  
              main=paste0("Mean ",groupNames[2])) 
    plotPost( ordMcmcMat[,"sigma[2]"] , 
              xlab=TeX(sprintf("$\\tau_{%s,%s}$", subscript, groupNames[2])) , 
              xlim=sigmaLim , compVal=compValSigma[2] , 
              col="#0088aa8a" , cex=cex_plotPost , 
              main=paste0("Scale ",groupNames[2]))
    plotPost( log10(ordMcmcMat[,"nu[2]"]) , 
              xlab=TeX(sprintf("$\\log_{10}\\nu_{%s,%s}$", subscript, groupNames[2])) , 
              compVal=compValNu[2] , ROPE=NULL , cex=cex_plotPost , col="#0088aa8a" , 
              main=paste0("Normality ",groupNames[2]))
    plotPost( ordMcmcMat[,"mu[2]"]-ordMcmcMat[,"mu[1]"] , compVal=compValMuDiff ,
              xlab=TeX(sprintf("$\\mu_{%s,%s} - \\mu_{%s,%s}$", subscript, groupNames[2], subscript, groupNames[1])) , 
              cex=cex_plotPost , col="#0088aa8a" , 
              main="Difference between means")
    plotPost( ordMcmcMat[,"sigma[2]"]-ordMcmcMat[,"sigma[1]"] , compVal=compValSigmaDiff , 
              xlab=TeX(sprintf("$\\tau_{%s,%s} - \\tau_{%s,%s}$", subscript, groupNames[2], subscript, groupNames[1])) , 
              cex=cex_plotPost , col="#0088aa8a" , 
              main="Difference between scales")
    
    # Posterior effect size:
    postEffSz = ( ( ordMcmcMat[,"mu[2]"] - ordMcmcMat[,"mu[1]"] ) 
                  / sqrt( ( ordMcmcMat[,"sigma[1]"]^2 + ordMcmcMat[,"sigma[2]"]^2 ) / 2 ) )
    
    plotPost( postEffSz , compVal=compValEff , ROPE=ropeEffSz , 
              col="#0088aa8a" , cex=cex_plotPost , cex.lab=1.2,
              xlab=TeX(sprintf("$d_{%s,%s}$", subsEffsz, subscript)) ,
              main="Effect size")
    saveGraph( paste0(saveName,"-Plots") , type=graphFileType )
  }
}

subBarplots = function (datFrm, yName="resp", qName="item", minLevel=1, 
                        groupNames=c(1,2), maxLevel=5, leg=TRUE) {
  
  # Getting number of items and groups from data:
  q = as.numeric(as.factor(datFrm[,qName]))
  g = as.numeric(as.factor(datFrm[,gName]))
  nQ = max(q) # Number of items.
  nG = max(g) # Number of groups.
  
  # Creating histogram with original data:
  y = as.numeric(datFrm[,yName])
  dataScaleDensMax = 0.7
  
  for ( gIdx in 1:nG ) {
    m = matrix(rep(0,nQ*maxLevel), nrow=nQ, ncol=maxLevel, byrow=TRUE)
    for ( qIdx in 1:nQ ) {
      for ( lIdx in 1:maxLevel ) {
        m[qIdx, lIdx] = sum(y[ g==gIdx & q==qIdx ]== lIdx)
      }
    }
    marromE ="#916e6eff"
    marromC = "#916e6ea4"
    azulE="#0088aaff"
    azulC= "#0087a864"
    amarelo="#ffe680ff"
    rosa="#f4e3d7ff"
    colors = c(amarelo, azulE, azulC, rosa, marromC, marromE)
    barplot(m, names.arg=seq(minLevel, maxLevel), yaxt='n', axis.lty=1 ,
            xlab=paste0("Responses ", groupNames[gIdx]) , cex.lab=1.2 , cex.axis = 3.0 ,
            col=colors, border=c(NA, NA, NA, NA), space=c(0.4,0.4,0.4,0.4,0.4),
            width=c(5.5,5.5,5.5,5.5,5.5) )
  }
  if (leg) {
    legend("top", legend = seq(1,nQ), horiz = TRUE, seg.len=1, 
           fill = colors, cex = 2.0, xpd = TRUE, 
           box.lty = 0, title="Item") 
  }
  
}
