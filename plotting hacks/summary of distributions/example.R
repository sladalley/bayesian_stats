graphics.off()
rm(list=ls())

source("plotFunctions.R")
load("example_metric_data.Rdata" )

# Setting null values and ROPE intervals:
nullHypValMu = 0.0 # Null hypothesis value to calculate effect size.
compValMu = 0.0 # Comparison value to show in mu's posterior plot.
ropeEffSz = c(-0.1, 0.1) # ROPE for the effect size.

# Defining directory for results and prefix of the files names:
fileNameRoot = "test"
dir.create(fileNameRoot)
fileNameRoot = paste0(fileNameRoot, "/", fileNameRoot)
graphFileType = "pdf"

# Plotting the posterior distributions:
plotObjectives(codaSamples, datFrm, yName=yName, gName=gName,
               nullHypValMu=nullHypValMu, compValMu=compValMu, ropeEffSz=ropeEffSz,
               graphFileType=graphFileType, saveName=fileNameRoot ,
               subscript="\\Delta t", subsEffsz = "obj", plotOption = 1, summ=TRUE)

graphics.off()