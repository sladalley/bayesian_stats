graphics.off()
rm(list=ls())

source("metric_robust_functions.R")
source("beta_functions.R")

# Generating a data set using a uniform distribution just to test the script:
datFrm = read.csv(file="iros_data_7.csv")

datFrm = subset(datFrm, 
                Solved == 1 & 
                  IsStable == 1 & 
                  ConstraintsOK == 1 & JointsOK == 1 & JointsVelOK == 1 & SimTime <= 110 &( ( Controller == "FBL-QP") |( Controller == "QP")))

gName = "Controller"
yName = "SimTime"

# Defining directory for results and prefix of the files names:
fileNameRoot = "iros_beta_ppc_ce"
dir.create(fileNameRoot)
fileNameRoot = paste0(fileNameRoot, "/", fileNameRoot)
graphFileType = "pdf"

analysis = "independent"
power = "true"
saveName = fileNameRoot
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
  #ropeMu = c(-0.05, 0.05)  
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
  
  #ropeMu = c(-0.05*mean(y[ g==2 ]), 0.05*mean(y[ g==2 ]))
  mean_raw = median(y[ g==2 ])
  print(mean(y[g==1]))
  print(mean(y[g==2]))

  #print(max(y))
  #print(min(y))
  
  y_norm <- numeric(length(y))
  print(length(y[g==1]))
  print(length(y[g==2]))
  eps <- 1e-6
  
  # Global min and max across ALL data
  y_min <- min(y)
  y_max <- max(y)
  
  # Apply scaling to all points
  y_scaled <- (y - y_min) / (y_max - y_min)
  mean_norm <- (mean_raw - y_min) / (y_max - y_min)
  
  # Push away from exact 0 and 1 for Beta
  y_norm <- (y_scaled * (1 - 2*eps)) + eps
  
  y <- y_norm
  ropeMu <- c(-0.05*mean_norm, 0.05*mean_norm)
  
  muPrior = c( mean(y[ g==1 ]), mean(y[ g==2 ]) ) # mean of means priors (normal distribution)
  muSdPrior = 100*c( sd(y[ g==1 ]), sd(y[ g==2 ]) ) # s.d. of means priors (normal distribution)
  sigmaPriorLow = c( 0, 0 ) # low limits of scales priors (uniform distribution)
  sigmaPriorHigh = c( sd(y[ g==1 ]), sd(y[ g==2 ]) ) # high limits of scales priors (uniform distribution)
  
 # ropeMu = c(-0.05*mean(y[ g==2 ]), 0.05*mean(y[ g==2 ]))  
  
  #print(ropeMu)
  
  c_mu <- 5   # prior strength (equivalent role to your 100 scaling, but stable)
  muAlphaPrior <- muPrior * c_mu
  muBetaPrior  <- (1 - muPrior) * c_mu
  
  c_mu <- 5  # prior strength (like your 100 scaling, but stable)
  
  muAlphaPrior <- muPrior * c_mu
  muBetaPrior  <- (1 - muPrior) * c_mu
  kappa_hat <- (muPrior * (1 - muPrior)) / (muSdPrior^2) - 1
  kappa_hat <- pmax(kappa_hat, 0.1)
  
  shapeKappa <- 0.5
  rateKappa  <- shapeKappa / kappa_hat
  
  print("Control Effort")
  # Generating MCMC chain CE:
  codaSamples = genMCMC_beta(datFrm, muAlphaPrior = muAlphaPrior, muBetaPrior = muBetaPrior, 
                             rateKappa = rateKappa, 
                             shapeKappa = shapeKappa,
                              saveName=fileNameRoot, yName=yName, gName = gName)
  save( codaSamples , file=paste(saveName,"_", yName,"-Mcmc_real.Rdata",sep="") )
  
  # Summary and diagnostics of the MCMC chain:
  smryMCMC_beta(codaSamples, nG=nG, saveName=fileNameRoot, 
                 diagnostics=diagnostics, graphFileType=graphFileType)
  
  # Plotting the posterior distributions:
  plotMCMC_beta(codaSamples, datFrm, compValMu=compValMu, 
                 compValMuDiff=compValMuDiff, compValSigma=compValSigma,
                 compValSigmaDiff=compValSigmaDiff, compValNu=compValNu, 
                 compValEff=compValEff, ropeEffSz=ropeEffSz,ropeSigma = ropeSigma, ropeMu=ropeMu,
                 graphFileType=graphFileType, saveName=fileNameRoot, 
                 groupNames=groupNames, yName=yName, gName = gName)
  mcmcMat_ce = as.matrix(codaSamples)
  res = ppc_beta(mcmcChain = mcmcMat_ce, N = c(length(y[g==1]), length(y[g==2])), yName = yName, gName = gName, groupNames=groupNames, saveName = fileNameRoot, realdat = datFrm)
#  graphics.off()
 #  powerEstimation_beta(mcmcChain = mcmcMat_ce, muNULL = meanPrior, effROPE = ropeEffSz, effHDImaxWid = 0.2, N = 600, yName = yName, gName = gName, groupNames=groupNames, muROPE = ropeMu, saveName = fileNameRoot)
 #  powerEstimation_beta(mcmcChain = mcmcMat_ce, muNULL = meanPrior, effROPE = ropeEffSz, effHDImaxWid = 0.2, N = 800, yName = yName, gName = gName, groupNames=groupNames, muROPE = ropeMu, saveName = fileNameRoot)
 #  powerEstimation_beta(mcmcChain = mcmcMat_ce, muNULL = meanPrior, effROPE = ropeEffSz, effHDImaxWid = 0.2, N = 1000, yName = yName, gName = gName, groupNames=groupNames, muROPE = ropeMu, saveName = fileNameRoot)
 # 
 # 
 #  print("Sim_Time")
 #  # Generating MCMC chain ST:
 #  yName = "SimTime"
 #  # Getting number of groups from data:
 #  g = as.numeric(as.factor(datFrm[,gName]))
 #  nG = max(g) # number of groups
 # 
 #  # Setting parameters for prior distributions:
 #  y = as.numeric(datFrm[,yName])
 # 
 #  #ropeMu = c(-0.05*mean(y[ g==2 ]), 0.05*mean(y[ g==2 ]))
 #  mean_raw = median(y[ g==2 ])
 #  #print(mean_raw)
 # 
 # # print(max(y))
 # # print(min(y))
 # 
 #  y_norm <- numeric(length(y))
 #  eps <- 1e-6
 # 
 #  # Global min and max across ALL data
 #  y_min <- min(y)
 #  y_max <- max(y)
 # 
 #  # Apply scaling to all points
 #  y_scaled <- (y - y_min) / (y_max - y_min)
 #  mean_norm <- (mean_raw - y_min) / (y_max - y_min)
 # 
 #  # Push away from exact 0 and 1 for Beta
 #  y_norm <- (y_scaled * (1 - 2*eps)) + eps
 # 
 #  y <- y_norm
 #  ropeMu <- c(-0.05*mean_norm, 0.05*mean_norm)
 # 
 #  muPrior = c( mean(y[ g==1 ]), mean(y[ g==2 ]) ) # mean of means priors (normal distribution)
 #  muSdPrior = 100*c( sd(y[ g==1 ]), sd(y[ g==2 ]) ) # s.d. of means priors (normal distribution)
 #  sigmaPriorLow = c( 0, 0 ) # low limits of scales priors (uniform distribution)
 #  sigmaPriorHigh = c( sd(y[ g==1 ]), sd(y[ g==2 ]) ) # high limits of scales priors (uniform distribution)
 # 
 #  # ropeMu = c(-0.05*mean(y[ g==2 ]), 0.05*mean(y[ g==2 ]))
 # 
 #  #print(ropeMu)
 # 
 #  c_mu <- 5   # prior strength (equivalent role to your 100 scaling, but stable)
 #  
 #  muAlphaPrior <- muPrior * c_mu
 #  muBetaPrior  <- (1 - muPrior) * c_mu
 # 
 #  c_mu <- 5  # prior strength (like your 100 scaling, but stable)
 # 
 #  muAlphaPrior <- muPrior * c_mu
 #  muBetaPrior  <- (1 - muPrior) * c_mu
 #  kappa_hat <- (muPrior * (1 - muPrior)) / (muSdPrior^2) - 1
 #  kappa_hat <- pmax(kappa_hat, 0.1)
 # 
 #  shapeKappa <- 0.5
 #  rateKappa  <- shapeKappa / kappa_hat
 # 
 # 
 #  # Generating MCMC chain ST:
 #  codaSamples = genMCMC_beta(datFrm, muAlphaPrior = muAlphaPrior, muBetaPrior = muBetaPrior,
 #                             rateKappa = rateKappa,
 #                             shapeKappa = shapeKappa,
 #                             saveName=fileNameRoot, yName=yName, gName = gName)
 #  save(codaSamples , file=paste(saveName,"_", yName,"-Mcmc_real.Rdata",sep="") )
 # 
 #  # Summary and diagnostics of the MCMC chain:
 #  smryMCMC_beta(codaSamples, nG=nG, saveName=fileNameRoot,
 #                diagnostics=diagnostics, graphFileType=graphFileType)
 # 
 #  # Plotting the posterior distributions:
 #  plotMCMC_beta(codaSamples, datFrm, compValMu=compValMu,
 #                compValMuDiff=compValMuDiff, compValSigma=compValSigma,
 #                compValSigmaDiff=compValSigmaDiff, compValNu=compValNu,
 #                compValEff=compValEff, ropeEffSz=ropeEffSz,ropeSigma = ropeSigma, ropeMu=ropeMu,
 #                graphFileType=graphFileType, saveName=fileNameRoot,
 #                groupNames=groupNames, yName=yName, gName = gName)
 #  mcmcMat_st = as.matrix(codaSamples)
 #  graphics.off()
 # 
 #  powerEstimation_beta(mcmcChain = mcmcMat_st, muNULL = meanPrior, effROPE = ropeEffSz, effHDImaxWid = 0.2, N = 600, yName = yName, groupNames=groupNames, muROPE = ropeMu, saveName = fileNameRoot)
 #  powerEstimation_beta(mcmcChain = mcmcMat_st, muNULL = meanPrior, effROPE = ropeEffSz, effHDImaxWid = 0.2, N = 800, yName = yName, groupNames=groupNames, muROPE = ropeMu, saveName = fileNameRoot)
 #  powerEstimation_beta(mcmcChain = mcmcMat_st, muNULL = meanPrior, effROPE = ropeEffSz, effHDImaxWid = 0.2, N = 1000, yName = yName, groupNames=groupNames, muROPE = ropeMu, saveName = fileNameRoot)
 #  

  
  
}