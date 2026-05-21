graphics.off()
rm(list=ls())
source("metric_robust_functions.R")
source("DBDA2E-utilities.R")
library(fitdistrplus)
library(gofgamma)
library(metRology)



library(latex2exp)
library(nortest)
#library(DescTools)
library(MASS)
library(actuar)     # Burr, log-logistic
library(statmod)    # inverse Gaussian if needed later

# =========================
# LOAD + FILTER DATA
# =========================
datFrm_full = read.csv(file="data_csv/iros_data_7.csv")


metrics <- c("SimTime", "RMSE", "ControlEffort")
eff_thr <- 10000
time_thr <- 50
min_n   <- 5
sd <- 0.1



datFrm = subset(datFrm_full, 
                Solved == 1 & 
                  IsStable == 1 & 
                  ConstraintsOK == 1 & JointsOK == 1 & JointsVelOK == 1 &( ( Controller == "FBL-QP") |( Controller == "QP")) &SimTime<=300)

# datFrm = subset(datFrm_full, 
#                 Solved == 1 & 
#                   IsStable == 1 & 
#                   ConstraintsOK == 1 & JointsOK == 1 & JointsVelOK == 1 &( ( Controller == "FBL-QP") |( Controller == "QP")) &SimTime<=300)

dat_bad <- subset(datFrm_full,
                  !(Solved == 1 &
                      IsStable == 1 &
                      ConstraintsOK == 1 &
                      JointsOK == 1 &
                      JointsVelOK == 1 &
                      
                      SimTime <= 300) &( ( Controller == "FBL-QP") |( Controller == "QP"))
)




diagnostic <- datFrm %>%
  group_by(Scene, Controller) %>%
  summarise(
    n = sum(!is.na(ControlEffort)),
    
    ControlEffort_mean = mean(ControlEffort, na.rm = TRUE),
    ControlEffort_sd = sd(ControlEffort, na.rm = TRUE),
    RMSE_mean = mean(RMSE, na.rm = TRUE),
    RMSE_sd = sd(RMSE, na.rm = TRUE),
    SimTime_mean       = mean(SimTime, na.rm = TRUE),
    SimTime_sd       = sd(SimTime, na.rm = TRUE),
    
    
    
    eff_flag  = any(ControlEffort_mean > eff_thr, na.rm = TRUE),
    time_flag = any(SimTime_mean > time_thr, na.rm = TRUE),
    low_sample_flag = n < min_n,
    g = 10,
    sd_flag = ((SimTime_sd/sqrt(g))/SimTime_mean > 0.1) || ((ControlEffort_sd/sqrt(g))/ControlEffort_mean > 0.1) || ((RMSE_sd/sqrt(g))/RMSE_mean > 0.1),
    
    .groups = "drop"
  ) %>%
  filter(
    eff_flag | time_flag | low_sample_flag | sd_flag
  )
dat_summary <- datFrm %>%
  group_by(Scene, Controller) %>%
  summarise(
    across(
      all_of(metrics),
      list(
        mean = ~mean(.x, na.rm = TRUE),
        sd   = ~sd(.x, na.rm = TRUE)
      )
    ),
    .groups = "drop"
  )

# =========================
# Find Scene Outliers
# =========================
threshold <- 10  # set your limit

result <- datFrm %>%
  group_by(Scene, Controller) %>%
  summarise(
    ControlEffort = mean(ControlEffort, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(ControlEffort > threshold)
# =========================
# HISTOGRAMS
# =========================
par(mfrow = c(3, 2))
colour = "red"

# Sim Time Histograms

tab_st <- with(datFrm, tapply(SimTime, list(Scene,Controller),mean))
tab_st <-data.frame(scene = rownames(tab_st),tab_st)
tab_st <-na.omit(tab_st)
thisY <- tab_st$FBL.QP
#thisY  = (datFrm$SimTime[datFrm$Controller == "FBL-QP"])
xLim = c( min(thisY)-0.1*(max(thisY)-min(thisY)) ,  max(thisY)+0.1*(max(thisY)-min(thisY)) )

xBreaks = seq( xLim[1] , xLim[2] ,
                          length=ceiling((xLim[2]-xLim[1])/(sd(thisY)/4)) )
histInfo = hist(thisY, breaks=xBreaks, plot=FALSE)
yMax = 1.2 * max( histInfo$density )
histInfo = hist( thisY , prob=TRUE , ylim=c(0,yMax) ,
                            xlim=xLim , breaks=xBreaks , xlab="Sim Time" , ylab="" ,
                            cex.lab=1.0 , col=colour , border="white" , yaxt="n",
                            main=paste0( "FBL-QP"))

binMids <- 0.5 * (histInfo$breaks[-1] + histInfo$breaks[-length(histInfo$breaks)])
counts <- histInfo$counts
idx <- counts > 0


text(x = binMids[idx],
     y = histInfo$density[idx],
     labels = counts[idx],
     pos = 3,
     cex = 0.8)

tab_st <- with(datFrm, tapply(SimTime, list(Scene,Controller),mean))
tab_st <-data.frame(scene = rownames(tab_st),tab_st)
tab_st <-na.omit(tab_st)
thisY <- tab_st$QP
 #thisY  = (datFrm$SimTime[datFrm$Controller == "QP"])
 xLim = c( min(thisY)-0.1*(max(thisY)-min(thisY)) ,  max(thisY)+0.1*(max(thisY)-min(thisY)) )
 xBreaks = seq( xLim[1] , xLim[2] ,
                length=ceiling((xLim[2]-xLim[1])/(sd(thisY)/4)) )
 histInfo = hist(thisY, breaks=xBreaks, plot=FALSE)
 yMax = 1.2 * max( histInfo$density )

 histInfo = hist( thisY , prob=TRUE , ylim=c(0,yMax) ,
                  xlim=xLim , breaks=xBreaks , xlab="Sim Time" , ylab="" ,
                  cex.lab=1.0 , col=colour , border="white" , yaxt="n",
                  main=paste0( "QP"))
 binMids <- 0.5 * (histInfo$breaks[-1] + histInfo$breaks[-length(histInfo$breaks)])
 counts <- histInfo$counts
 idx <- counts > 0


 text(x = binMids[idx],
      y = histInfo$density[idx],
      labels = counts[idx],
      pos = 3,
      cex = 0.8)
 
 
 tab_st <- with(datFrm, tapply(ControlEffort, list(Scene,Controller),mean))
 tab_st <-data.frame(scene = rownames(tab_st),tab_st)
 tab_st <-na.omit(tab_st)
 thisY <- tab_st$FBL.QP
 #thisY  = (datFrm$ControlEffort[datFrm$Controller == "FBL-QP"])
 xLim = c( min(thisY)-0.1*(max(thisY)-min(thisY)) ,  max(thisY)+0.1*(max(thisY)-min(thisY)) )
 
 xBreaks = seq( xLim[1] , xLim[2] ,
                length=ceiling((xLim[2]-xLim[1])/(sd(thisY)/4)) )
 histInfo = hist(thisY, breaks=xBreaks, plot=FALSE)
 yMax = 1.2 * max( histInfo$density )
 histInfo = hist( thisY , prob=TRUE , ylim=c(0,yMax) ,
                  xlim=xLim , breaks=xBreaks , xlab="Control Effort" , ylab="" ,
                  cex.lab=1.0 , col=colour , border="white" , yaxt="n",
                  main=paste0( "FBL-QP"))
 
 binMids <- 0.5 * (histInfo$breaks[-1] + histInfo$breaks[-length(histInfo$breaks)])
 counts <- histInfo$counts
 idx <- counts > 0
 
 
 text(x = binMids[idx],
      y = histInfo$density[idx],
      labels = counts[idx],
      pos = 3,
      cex = 0.8)
 
 tab_st <- with(datFrm, tapply(ControlEffort, list(Scene,Controller),mean))
 tab_st <-data.frame(scene = rownames(tab_st),tab_st)
 tab_st <-na.omit(tab_st)
 thisY <- tab_st$QP
# thisY  = (datFrm$ControlEffort[datFrm$Controller == "QP"])
 xLim = c( min(thisY)-0.1*(max(thisY)-min(thisY)) ,  max(thisY)+0.1*(max(thisY)-min(thisY)) )
 xBreaks = seq( xLim[1] , xLim[2] ,
                length=ceiling((xLim[2]-xLim[1])/(sd(thisY)/4)) )
 histInfo = hist(thisY, breaks=xBreaks, plot=FALSE)
 yMax = 1.2 * max( histInfo$density )
 
 histInfo = hist( thisY , prob=TRUE , ylim=c(0,yMax) ,
                  xlim=xLim , breaks=xBreaks , xlab="Control Effort" , ylab="" ,
                  cex.lab=1.0 , col=colour , border="white" , yaxt="n",
                  main=paste0( "FBL-QP - QP"))
 binMids <- 0.5 * (histInfo$breaks[-1] + histInfo$breaks[-length(histInfo$breaks)])
 counts <- histInfo$counts
 idx <- counts > 0
 
 
 text(x = binMids[idx],
      y = histInfo$density[idx],
      labels = counts[idx],
      pos = 3,
      cex = 0.8)
 
 tab_st <- with(datFrm, tapply(RMSE, list(Scene,Controller),mean))
 tab_st <-data.frame(scene = rownames(tab_st),tab_st)
 tab_st <-na.omit(tab_st)
 thisY <- tab_st$FBL.QP
# thisY  = (datFrm$RMSE[datFrm$Controller == "FBL-QP"])
 xLim = c( min(thisY)-0.1*(max(thisY)-min(thisY)) ,  max(thisY)+0.1*(max(thisY)-min(thisY)) )
 
 xBreaks = seq( xLim[1] , xLim[2] ,
                length=ceiling((xLim[2]-xLim[1])/(sd(thisY)/4)) )
 histInfo = hist(thisY, breaks=xBreaks, plot=FALSE)
 yMax = 1.2 * max( histInfo$density )
 histInfo = hist( thisY , prob=TRUE , ylim=c(0,yMax) ,
                  xlim=xLim , breaks=xBreaks , xlab="RMSE" , ylab="" ,
                  cex.lab=1.0 , col=colour , border="white" , yaxt="n",
                  main=paste0( "FBL-QP"))
 
 binMids <- 0.5 * (histInfo$breaks[-1] + histInfo$breaks[-length(histInfo$breaks)])
 counts <- histInfo$counts
 idx <- counts > 0
 
 
 text(x = binMids[idx],
      y = histInfo$density[idx],
      labels = counts[idx],
      pos = 3,
      cex = 0.8)
 
 tab_st <- with(datFrm, tapply(RMSE, list(Scene,Controller),mean))
 tab_st <-data.frame(scene = rownames(tab_st),tab_st)
 tab_st <-na.omit(tab_st)
 thisY <- tab_st$QP
# thisY  = (datFrm$RMSE[datFrm$Controller == "QP"])
 xLim = c( min(thisY)-0.1*(max(thisY)-min(thisY)) ,  max(thisY)+0.1*(max(thisY)-min(thisY)) )
 xBreaks = seq( xLim[1] , xLim[2] ,
                length=ceiling((xLim[2]-xLim[1])/(sd(thisY)/4)) )
 histInfo = hist(thisY, breaks=xBreaks, plot=FALSE)
 yMax = 1.2 * max( histInfo$density )
 
 histInfo = hist( thisY , prob=TRUE , ylim=c(0,yMax) ,
                  xlim=xLim , breaks=xBreaks , xlab="RMSE" , ylab="" ,
                  cex.lab=1.0 , col=colour , border="white" , yaxt="n",
                  main=paste0( " QP"))
 binMids <- 0.5 * (histInfo$breaks[-1] + histInfo$breaks[-length(histInfo$breaks)])
 counts <- histInfo$counts
 idx <- counts > 0
 
 
 text(x = binMids[idx],
      y = histInfo$density[idx],
      labels = counts[idx],
      pos = 3,
      cex = 0.8)

# # thisY_fbl_qp = (datFrm$SimTime[datFrm$Controller == "FBL-QP"]) 
# # thisY_qp = (datFrm$SimTime[datFrm$Controller == "QP"]) 
# # 
# # 
# # thisY = (datFrm$ControlEffort[datFrm$Controller == "FBL-QP"])
# # print("FBL QP Length")
# # print(length(thisY))
# # xLim = c( min(thisY)-0.1*(max(thisY)-min(thisY)) , 
# #           max(thisY)+0.1*(max(thisY)-min(thisY)) )
# # xBreaks = seq( xLim[1] , xLim[2] , 
# #                length=ceiling((xLim[2]-xLim[1])/(sd(thisY)/4)) )
# # histInfo = hist(thisY, breaks=xBreaks, plot=FALSE)
# # yMax = 1.2 * max( histInfo$density )
# # 
# # histInfo = hist( thisY , prob=TRUE , ylim=c(0,yMax) , 
# #                  xlim=xLim , breaks=xBreaks , xlab="Control Effort" , ylab="" , 
# #                  cex.lab=1.0 , col=colour , border="white" , yaxt="n" ,
# #                  main=paste0( "FBL-QP"))
# # #lines(density(thisY), col = "blue", lwd = 2)
# # 
# # thisY = (datFrm$ControlEffort[datFrm$Controller == "QP"])
# # print("QP Length")
# # print(length(thisY))
# # xLim = c( min(thisY)-0.1*(max(thisY)-min(thisY)) , 
# #           max(thisY)+0.1*(max(thisY)-min(thisY)) )
# # xBreaks = seq( xLim[1] , xLim[2] , 
# #                length=ceiling((xLim[2]-xLim[1])/(sd(thisY)/4)) )
# # histInfo = hist(thisY, breaks=xBreaks, plot=FALSE)
# # yMax = 1.2 * max( histInfo$density )
# # 
# # histInfo = hist( thisY , prob=TRUE , ylim=c(0,yMax) , 
# #                  xlim=xLim , breaks=xBreaks , xlab="Control Effort" , ylab="" , 
# #                  cex.lab=1.0 , col=colour , border="white" , yaxt="n" ,
# #                  main=paste0( "QP"))
# # #lines(density(thisY), col = "blue", lwd = 2)
# # 
# # thisY = ((datFrm$SimTime[datFrm$Controller == "FBL-QP"]))
# # xLim = c( min(thisY)-0.1*(max(thisY)-min(thisY)) , 
# #           max(thisY)+0.1*(max(thisY)-min(thisY)) )
# # xBreaks = seq( xLim[1] , xLim[2] , 
# #                length=ceiling((xLim[2]-xLim[1])/(sd(thisY)/4)) )
# # histInfo = hist(thisY, breaks=xBreaks, plot=FALSE)
# # yMax = 1.2 * max( histInfo$density )
# # 
# # histInfo = hist( thisY , prob=TRUE , ylim=c(0,yMax) , 
# #                  xlim=xLim , breaks=xBreaks , xlab="Convergence Time" , ylab="" , 
# #                  cex.lab=1.0 , col=colour , border="white" , yaxt="n" ,
# #                  main=paste0( "FBL-QP"))
# # #lines(density(thisY), col = "blue", lwd = 2)
# # 
# # thisY = (datFrm$SimTime[datFrm$Controller == "QP"])
# # xLim = c( min(thisY)-0.1*(max(thisY)-min(thisY)) , 
# #           max(thisY)+0.1*(max(thisY)-min(thisY)) )
# # xBreaks = seq( xLim[1] , xLim[2] , 
# #                length=ceiling((xLim[2]-xLim[1])/(sd(thisY)/4)) )
# # histInfo = hist(thisY, breaks=xBreaks, plot=FALSE)
# # yMax = 1.2 * max( histInfo$density )
# # 
# # histInfo = hist( thisY , prob=TRUE , ylim=c(0,yMax) , 
# #                  xlim=xLim , breaks=xBreaks , xlab="Convergence Time" , ylab="" , 
# #                  cex.lab=1.0 , col=colour , border="white" , yaxt="n" ,
# #                  main=paste0( "QP"))
# # #lines(density(thisY), col = "blue", lwd = 2)
# # 
# # thisY = datFrm$RMSE[datFrm$Controller == "FBL-QP"]
# # xLim = c( min(thisY)-0.1*(max(thisY)-min(thisY)) , 
# #           max(thisY)+0.1*(max(thisY)-min(thisY)) )
# # xBreaks = seq( xLim[1] , xLim[2] , 
# #                length=ceiling((xLim[2]-xLim[1])/(sd(thisY)/4)) )
# # histInfo = hist(thisY, breaks=xBreaks, plot=FALSE)
# # yMax = 1.2 * max( histInfo$density )
# # 
# # histInfo = hist( thisY , prob=TRUE , ylim=c(0,yMax) , 
# #                  xlim=xLim , breaks=xBreaks , xlab="RMSE" , ylab="" , 
# #                  cex.lab=1.0 , col=colour , border="white" , yaxt="n" ,
# #                  main=paste0( "FBL-QP"))
# # #lines(density(thisY), col = "blue", lwd = 2)
# # 
# # thisY = datFrm$RMSE[datFrm$Controller == "QP"]
# # xLim = c( min(thisY)-0.1*(max(thisY)-min(thisY)) , 
# #           max(thisY)+0.1*(max(thisY)-min(thisY)) )
# # xBreaks = seq( xLim[1] , xLim[2] , 
# #                length=ceiling((xLim[2]-xLim[1])/(sd(thisY)/4)) )
# # histInfo = hist(thisY, breaks=xBreaks, plot=FALSE)
# # yMax = 1.2 * max( histInfo$density )
# # 
# # histInfo = hist( thisY , prob=TRUE , ylim=c(0,yMax) , 
# #                  xlim=xLim , breaks=xBreaks , xlab="RMSE" , ylab="" , 
# #                  cex.lab=1.0 , col=colour , border="white" , yaxt="n" ,
# #                  main=paste0( "QP"))
# #lines(density(thisY), col = "blue", lwd = 2)
# 
# #hist(datFrm$ControlEffort[datFrm$Controller == "FBL-QP"],
# #     main = "FBL-QP", xlab = "Control Effort", col = "skyblue")
# 
# #hist(datFrm$ControlEffort[datFrm$Controller == "QP"],
# #     main = "QP", xlab = "Control Effort", col = "salmon")
# 
# #hist(datFrm$SimTime[datFrm$Controller == "FBL-QP"],
# #     main = "FBL-QP", xlab = "Sim Time", col = "skyblue")
# 
# #hist(datFrm$SimTime[datFrm$Controller == "QP"],
# #     main = "QP", xlab = "Sim Time", col = "salmon")
# 
# #hist(datFrm$RMSE[datFrm$Controller == "FBL-QP"],
# #     main = "FBL-QP", xlab = "RMSE", col = "skyblue")
# 
# #hist(datFrm$RMSE[datFrm$Controller == "QP"],
# #     main = "QP", xlab = "RMSE", col = "salmon")
# 
# # =========================
# # DATA (RMSE - FBL-QP)
# # =========================
# x_raw <- datFrm$SimTime[datFrm$Controller == "FBL-QP"]
# x <- x_raw[is.finite(x_raw)]
# 
# # =========================
# # RESULTS TABLE
# # =========================
# results <- data.frame(
#   Distribution = character(),
#   AD_Statistic = numeric(),
#   p_value = numeric(),
#   stringsAsFactors = FALSE
# )
# 
# # =========================
# # SAFE AD WRAPPER
# # =========================
# safe_ad <- function(x, nullfun, params, name) {
#   
#   out <- tryCatch({
#     do.call(
#       AndersonDarlingTest,
#       c(list(x = x, null = nullfun, nullname = name), params)
#     )
#   }, warning = function(w) NULL,
#   error = function(e) NULL)
#   
#   if (is.null(out) || is.na(out$statistic)) {
#     return(data.frame(
#       Distribution = name,
#       AD_Statistic = NA,
#       p_value = NA
#     ))
#   }
#   
#   data.frame(
#     Distribution = name,
#     AD_Statistic = as.numeric(out$statistic),
#     p_value = out$p.value
#   )
# }
# 
# safe_ad_generic <- function(x, cdf, params, name) {
#   
#   out <- tryCatch({
#     do.call(
#       AndersonDarlingTest,
#       c(list(x = x, null = cdf, nullname = name), params)
#     )
#   }, warning = function(w) NULL,
#   error = function(e) NULL)
#   
#   if (is.null(out) || is.na(out$statistic)) {
#     return(data.frame(Distribution=name,
#                       AD_Statistic=NA,
#                       p_value=NA))
#   }
#   
#   data.frame(Distribution=name,
#              AD_Statistic=as.numeric(out$statistic),
#              p_value=out$p.value)
# }
# 
# # =========================
# # 1. NORMAL
# # =========================
# results <- rbind(results,
#                  safe_ad(
#                    x,
#                    pnorm,
#                    list(mean = mean(x), sd = sd(x)),
#                    "Normal"
#                  )
# )
# 
# # =========================
# # 2. LOGNORMAL (strict positive)
# # =========================
# x_pos <- x[is.finite(x) & x > 0]
# mean <- 6.246304 
# sigma <- 0.900862
# 
# if (length(x_pos) > 5) {
#   
#   results <- rbind(results,
#                    safe_ad(
#                      x_pos,
#                      plnorm,
#                      list(meanlog = mean,
#                           sdlog = sigma),
#                      "Lognormal"
#                    )
#   )
# }
# 
# # =========================
# # 3. EXPONENTIAL
# # =========================
# if (length(x_pos) > 5) {
#   
#   rate <- 1 / mean(x_pos)
#   
#   results <- rbind(results,
#                    safe_ad(
#                      x_pos,
#                      pexp,
#                      list(rate = rate),
#                      "Exponential"
#                    )
#   )
# }
# 
# # =========================
# # 4. WEIBULL (FULLY ROBUST FIX)
# # =========================
# x_w <- x[is.finite(x) & x > 0]
# 
# # remove extreme numerical instability
# if (length(x_w) > 10) {
#   shape <-0.933276436
#   scale <- 826.0508994
#   
#   x_w <- x_w[x_w > quantile(x_w, 0.01)]
#   
#   if (sd(x_w) > 0) {
#     
#     fit_w <- tryCatch({
#       MASS::fitdistr(x_w, "weibull")
#     }, error = function(e) NULL,
#     warning = function(w) NULL)
#     
#     # fallback if fitdistr fails
#     if (is.null(fit_w) ||
#         any(!is.finite(fit_w$estimate))) {
#       
#       #shape <- 1.2 / sd(log(x_w))
#       #scale <- mean(x_w) / gamma(1 + 1/shape)
#       
#     } else {
#       #shape <- fit_w$estimate["shape"]
#       #scale <- fit_w$estimate["scale"]
#     }
#     #shape <- mean(x)
#     #scale  <- sd(x)
#     # final safety check before AD test
#     if (is.finite(shape) && is.finite(scale) &&
#         shape > 0 && scale > 0) {
#       
#       results <- rbind(results,
#                        safe_ad(
#                          x_w,
#                          pweibull,
#                          list(shape = shape, scale = scale),
#                          "Weibull"
#                        )
#       )
#       
#     } else {
#       results <- rbind(results,
#                        data.frame(Distribution="Weibull",
#                                   AD_Statistic=NA,
#                                   p_value=NA)
#       )
#     }
#     
#   } else {
#     results <- rbind(results,
#                      data.frame(Distribution="Weibull",
#                                 AD_Statistic=NA,
#                                 p_value=NA)
#     )
#   }
#   
# } else {
#   results <- rbind(results,
#                    data.frame(Distribution="Weibull",
#                               AD_Statistic=NA,
#                               p_value=NA)
#   )
# }
# 
# # =========================
# # 5. UNIFORM (FIXED PROPERLY)
# # =========================
# x_u <- x[is.finite(x)]
# 
# unif_min <- min(x_u)
# unif_max <- max(x_u)
# 
# results <- rbind(results,
#                  safe_ad(
#                    x_u,
#                    punif,
#                    list(min = unif_min, max = unif_max),
#                    "Uniform"
#                  )
# )
# # =========================
# # 5. T (FIXED PROPERLY)
# # =========================
# mu <- 488.6456771
# sigma <- 187.70203215
# 
# # safe df estimate
# df <- 1.162036
# 
# results <- rbind(results,
#                  safe_ad_generic(
#                    x,
#                    function(q, df, mu, sigma) pt((q - mu)/sigma, df),
#                    list(df = df, mu = mu, sigma = sigma),
#                    "Student-t"
#                  )
# )
# # =========================
# # 5. Log Logistic
# # =========================
# x_p <- x[is.finite(x) & x > 0]
# 
# if (length(x_p) > 5) {
#   
#   shape <- 1.2 / sd(log(x_p))
#   scale <- exp(mean(log(x_p)))
#   
#   results <- rbind(results,
#                    safe_ad_generic(
#                      x_p,
#                      function(q, shape, scale) pllogis(q, shape = shape, scale = scale),
#                      list(shape = shape, scale = scale),
#                      "Log-logistic"
#                    )
#   )
# }
# # =========================
# # FINAL SORTED OUTPUT
# # =========================
# x_raw <- (datFrm$SimTime[datFrm$Controller == "FBL-QP"])
# x_min <- min(x_raw)
# x_max <- max(x_raw)
# x_beta <- (x_raw - x_min) / (x_max - x_min)
# eps <- 1e-6
# x_beta <- (x_beta * (1 - 2*eps)) + eps
# results <- results[order(results$AD_Statistic, na.last = TRUE), ]
# print(results)
# # cf = descdist(x_beta, discrete=FALSE, boot=1000)
# fit.lnorm = fitdist(x_beta, "lnorm")
# fit.weibull = fitdist(x_beta, "weibull", lower = c(0, 0))
# fit.norm = fitdist(x_beta,"norm")
# fit.t = fitdist(x_beta,"t.scaled",
#                     start=list(df=5,mean=mean(x_raw),sd=sd(x_raw)))
# fit.gamma = fitdist(x_beta, "gamma", lower = c(0, 0), start = list(scale = 1, shape = 1))
# fit.beta = fitdist(x_beta, "beta", lower = c(0, 0),start = list(shape1 = 0.62, shape2 = 6.6))
# 
# 
# print(summary(fit.norm))
# print(summary(fit.lnorm))
# print(summary(fit.weibull))
# print(summary(fit.t))
# print(summary(fit.gamma))
# print(summary(fit.beta))
# 
# list_dist = list(fit.norm, fit.lnorm, fit.weibull, fit.t, fit.gamma, fit.beta)
# 
# # plot.legend <- c("normal", "lognormal", "weibull", "t", "gamma", "beta")
# # denscomp(list_dist, legendtext = plot.legend)
# # qqcomp(list_dist, legendtext = plot.legend)
# # cdfcomp(list_dist, legendtext = plot.legend)
# # ppcomp(list_dist, legendtext = plot.legend)
# 
# print(gofstat(list_dist))
# 
# #plot(fit.beta, las = 1)
# 
# #print(x_log)
# #print(x)
# #plot(fit.lnorm)
