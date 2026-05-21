graphics.off()
rm(list=ls())
library(tidyr)
library(dplyr)

source("metric_robust_functions.R")
source("lognormal_functions.R")


# Generating a data set using a uniform distribution just to test the script:
# Generating a data set using a uniform distribution just to test the script:
datFrm_full = read.csv(file="data_csv/iros_data_7.csv")





datFrm = subset(datFrm_full, 
                Solved == 1 & 
                  IsStable == 1 & 
                  ConstraintsOK == 1 & JointsOK == 1 & JointsVelOK == 1 &( ( Controller == "FBL-QP") |( Controller == "QP"))  &SimTime<=300)


# =========================
# HISTOGRAMS
# =========================

gName = "Controller"
yName = "RMSE"

g = as.numeric(as.factor(datFrm[,gName]))


nG = max(g) # number of groups

if (yName == "SuccessProbability")
{
  tab <- datFrm %>%
    group_by(Scene, Controller) %>%
    summarise(
      Success :=  n(),
      .groups = "drop"
    ) %>%
    group_by(Scene) %>%
    filter(n() == nG) %>%   
    ungroup()
} else {

tab <- datFrm %>%
  group_by(Scene, Controller) %>%
  summarise(
    "{yName}" := mean(.data[[yName]], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(Scene) %>%
  filter(n() == nG) %>%   
  ungroup()

dat_summary <- datFrm %>%
  group_by(Scene, Controller) %>%
  summarise(
    mean = mean(.data[[yName]], na.rm = TRUE),
    sd   = sd(.data[[yName]], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from  = Controller,
    values_from = c(mean, sd)
  )

names(tab) <- gsub("\\.", "-", names(tab))
names(dat_summary) <- gsub("\\.", "-", names(dat_summary))

}

names(tab) <- gsub("\\.", "-", names(tab))
datFrm <- as.data.frame(tab)




# Defining directory for results and prefix of the files names:
fileNameRoot = "paired_rmse"
dir.create(fileNameRoot)
fileNameRoot = paste0(fileNameRoot, "/", fileNameRoot)
graphFileType = "pdf"
Scenes = unique(datFrm$Scene)


summary_stats <- datFrm %>%
  group_by(Controller) %>%
  summarise(
    median = median(.data[[yName]], na.rm = TRUE),
    mean   = mean(.data[[yName]], na.rm = TRUE),
    sd = mean(.data[[yName]], na.rm = TRUE)
  )






analysis = "paired"
# Choosing type of analysis ("paired" or "independent"):
name = "FBL-QP - QP"
# ======= Paired/one-group analysis ===========================================
if (analysis == "paired") {
  
  # Generating data frame with differences between groups:
  diff = data.frame(temp=numeric(), Controller=character(), Scene=integer())
  names(diff)[1] <- yName
  for (Scene in Scenes) {

    y1 = datFrm[datFrm$Scene == Scene & datFrm$Controller == "FBL-QP", yName]
    y2 = datFrm[datFrm$Scene == Scene & datFrm$Controller == "QP", yName]
    diff[nrow(diff) + 1,] <- list(y1 - y2, name, Scene)
  }
  datFrm = diff
  
  summary_stats$diff_mean = mean(datFrm[,yName])
  summary_stats$diff_median = median(datFrm[,yName])
  summary_stats$diff_sd = sd(datFrm[,yName])
  
  ropeMu = c(-0.1*abs(summary_stats$diff_median[summary_stats$Controller=="QP"]),0.1*abs(summary_stats$diff_median[summary_stats$Controller=="QP"]))
  
  
  # Setting the comparison values and ROPE intervals:
  compValMu = NULL # comparison value to show in plots of means
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
                        saveName=fileNameRoot, yName=yName, gName = gName)
 
  
  # Summary and diagnostics of the MCMC chain:
  smryMCMC(codaSamples, nG=nG, nullValEff=nullValEff, saveName=fileNameRoot, 
           diagnostics=diagnostics, graphFileType=graphFileType, 
           computeEffsz=computeEffsz)
  
  # Plotting the posterior distributions:
  plotMCMC(codaSamples, datFrm, nullValEff=nullValEff, compValMu=compValMu,
           compValSigma=compValSigma,
           compValNu=compValNu, 
           compValEff=compValEff, ropeEffSz=ropeEffSz, 
           graphFileType=graphFileType, saveName=fileNameRoot, ropeMu=ropeMu,
           yName=yName, gName = gName,
           subscript=subscript, subsEffsz=subsEffsz, plotEffsz=computeEffsz)
}
# ======= Independent/two-groups analysis =====================================
if (analysis == "independent") {
  # Setting the comparison values and ROPE intervals:
  compValMu = NULL # comparison value to show in plots of means
  compValMuDiff = NULL # comparison value to show in plot of means difference
  compValSigma = NULL # comparison value to show in plots of scales
  compValSigmaDiff = NULL # comparison value to show in plot of scales difference
  compValNu = NULL # comparison value to show in plots of normality
  compValEff = NULL # comparison value to show in plot of effect size
  ropeEffSz = c(-0.1, 0.1) # ROPE around null value of effect size
  ropeMu = NULL
  ropeSigma = NULL
  
  diagnostics = TRUE # generate diagnostics plots or not
  computeEffsz = TRUE # calculate and plot effect size or not
  groupNames = c("FBL-QP","QP") # names of the groups
  subscript = "" # subscript to identify measured variable
  subsEffsz = "" # subscript to identify the effect size expression
  
  MuRopeValue = 0.1*abs(summary_stats$median[summary_stats$Controller == "QP"])
  
  ropeMu = c(-MuRopeValue, MuRopeValue)
  
  # Getting number of groups from data:
  g = as.numeric(as.factor(datFrm[,gName]))
  nG = max(g) # number of groups
  
  # Setting parameters for prior distributions:
  y = as.numeric(datFrm[,yName])
  muPrior = c( mean(y[ g==1 ]), mean(y[ g==2 ]) ) # mean of means priors (normal distribution)
  muSdPrior = 100*c( sd(y[ g==1 ]), sd(y[ g==2 ]) ) # s.d. of means priors (normal distribution)
  sigmaPriorLow = c( sd(y[ g==1 ]), sd(y[ g==2 ]) )/1000 # low limits of scales priors (uniform distribution)
  sigmaPriorHigh = 1000*c( sd(y[ g==1 ]), sd(y[ g==2 ]) ) # high limits of scales priors (uniform distribution)
  
  #ropeMu = c(-0.05*mean(y[ g==2 ]), 0.05*mean(y[ g==2 ]))
  
 
  
  # Generating MCMC chain:
  codaSamples = genMCMC(datFrm, muPrior=muPrior, muSdPrior=muSdPrior, 
                        sigmaPriorLow=sigmaPriorLow, 
                        sigmaPriorHigh=sigmaPriorHigh,
                        saveName=fileNameRoot, yName=yName, gName = gName)
  
  # Summary and diagnostics of the MCMC chain:
  smryMCMC(codaSamples, nG=nG, saveName=fileNameRoot, 
           diagnostics=diagnostics, graphFileType=graphFileType, 
           computeEffsz=computeEffsz)
  
  # Plotting the posterior distributions:
  plotMCMC(codaSamples, datFrm, compValMu=compValMu, 
           compValMuDiff=compValMuDiff, compValSigma=compValSigma,
           compValSigmaDiff=compValSigmaDiff, compValNu=compValNu, 
           compValEff=compValEff, ropeEffSz=ropeEffSz,ropeSigma = ropeSigma, ropeMu=ropeMu,
           graphFileType=graphFileType, saveName=fileNameRoot, 
           groupNames=groupNames, yName=yName, gName = gName,
           subscript=subscript, subsEffsz=subsEffsz, plotEffsz=computeEffsz)
}
save(summary_stats, file = paste0(fileNameRoot, "_summary_stats.RData"))
# =============================================================================
