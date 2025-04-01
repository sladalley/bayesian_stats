# Functions for exploratory data analysis

# Checking that required packages are installed:
want = c("ggpubr")
have = want %in% rownames(installed.packages())
if ( any(!have) ) { install.packages( want[!have] ) }

# ======== Functions from Kruschke's DBDA2E-utilities.R =======================

openGraph = function( width=7 , height=7 , mag=1.0 , ... ) {
  if ( .Platform$OS.type != "windows" ) { # Mac OS, Linux
    tryInfo = try( X11( width=width*mag , height=height*mag , type="cairo" , 
                        ... ) )
    if ( class(tryInfo)=="try-error" ) {
      lineInput = readline("WARNING: Previous graphics windows will be closed because of too many open windows.\nTO CONTINUE, PRESS <ENTER> IN R CONSOLE.\n")
      graphics.off() 
      X11( width=width*mag , height=height*mag , type="cairo" , ... )
    }
  } else { # Windows OS
    tryInfo = try( windows( width=width*mag , height=height*mag , ... ) )
    if ( class(tryInfo)=="try-error" ) {
      lineInput = readline("WARNING: Previous graphics windows will be closed because of too many open windows.\nTO CONTINUE, PRESS <ENTER> IN R CONSOLE.\n")
      graphics.off() 
      windows( width=width*mag , height=height*mag , ... )    
    }
  }
}

saveGraph = function( file="saveGraphOutput" , type="pdf" , ... ) {
  if ( .Platform$OS.type != "windows" ) { # Mac OS, Linux
    if ( any( type == c("png","jpeg","jpg","tiff","bmp")) ) {
      sptype = type
      if ( type == "jpg" ) { sptype = "jpeg" }
      savePlot( file=paste0(file,".",type) , type=sptype , ... )     
    }
    if ( type == "pdf" ) {
      dev.copy2pdf(file=paste0(file,".",type) , ... )
    }
    if ( type == "eps" ) {
      dev.copy2eps(file=paste0(file,".",type) , ... )
    }
  } else { # Windows OS
    file=paste0(file,".",type) 
    savePlot( file=file , type=type , ... )
  }
}
# =============================================================================

# Calculate general statistics of data. 
# It prints the results and creates a CSV file if a file name is given.
summary_info = function(datFrm, yName, gName=NULL, 
                        groupNames=NULL, saveName=NULL) {
  # List of parameters:
  #   - datFrm: data set.
  #   - yName: name of the column with the data.
  #   - gName: name of the column with the group indexes.
  #   - groupNames: names of the groups.
  #   - saveName: prefix of the output file.
  
  # Getting number of groups from data:
  if (is.null(gName)) {
    nG = 1
  } else {
    g = as.numeric(as.factor(datFrm[,gName])) # group data vector
    nG = max(g) # number of groups
  }
  
  if (nG == 1) {
    # Single group
    
    # Obtaining summary data and creating data frame with it:
    summaryInfo = summary(datFrm[,yName])
    df = data.frame(matrix(summaryInfo, nrow=1))
    colnames(df) = c("Min.","1st Qu.","Median", "Mean", "3rd Qu.", "Max.")
    
    # Adding standard deviation information:
    df = cbind(df, data.frame(Sd = sd(datFrm[,yName])))
    
    print(df)
    if (!is.null(saveName))
      write.csv( df , file=paste(saveName,"-SummaryInfo.csv",sep="") )
    
  } else {
    # Multiple groups
    
    # Defining groups names, if not given as parameter:
    if (is.null(groupNames))
      groupNames = levels(factor(g))
    
    df = NULL
    rowNames = NULL
    for (i in seq(1,nG)) {
      # Obtaining summary info and standard deviation:
      df = rbind(df, c(summary(datFrm[g==i, yName]), sd(datFrm[g==i, yName])))
    }
    # Setting names of columns and rows:
    colnames(df) = c("Min.","1st Qu.","Median", "Mean", "3rd Qu.", "Max.", "Sd.")
    row.names(df) = groupNames
    
    print(df)
    if (!is.null(saveName))
      write.csv( df , file=paste(saveName,"-SummaryInfo.csv",sep="") )
  }
}

# Plot boxplots with the data.
# It saves an image file if a file name is given.
plot_boxplots = function(datFrm, yName, gName=NULL, sName=NULL, 
                         groupNames=NULL, connect=FALSE, plotTitle=NULL, 
                         plotXlab=NULL, plotYlab=NULL, colors=NULL, 
                         graphFileType="png", saveName=NULL) {
  # List of parameters:
  #   - datFrm: data set.
  #   - yName: name of the column with the data.
  #   - gName: name of the column with the group indexes.
  #   - sName: name of the column with the subject indexes.
  #   - groupNames: names of each group.
  #   - connect: if data from each subject should be connected or not.
  #   - plotTitle: title of the boxplot.
  #   - plotXlab: x label in the boxplot (groups/conditions).
  #   - plotYlab: y label in the boxplot (variable/value).
  #   - graphFileType: type of the image output files.
  #   - saveName: prefix of the output files.
  
  # Getting number of groups from data:
  if (is.null(gName)) {
    nG = 1
  } else {
    g = as.numeric(as.factor(datFrm[,gName])) # group data vector
    nG = max(g) # number of groups
  }
  
  if (!is.null(gName)) {
    s = as.numeric(datFrm[,sName]) # subjects data vector
    nS = max(s) # number of subjects
  }
  
  library(ggpubr)
  if (nG == 1) {
    # Single group
    
    # Defining plot colors, if not given as parameter:
    if (is.null(colors))
      colors = "#0088aa8a"
    
    # Plotting boxplots
    openGraph(width=7,height=7)
    plot <- ggboxplot(datFrm[,yName], fill=colors, title=plotTitle, 
                      xlab=plotXlab,ylab=plotYlab) + theme(axis.text.x = element_blank())
    print(plot)
    
  } else {
    # Multiple groups
    
    # Defining plot colors, if not given as parameter:
    if (is.null(colors)) {
      colors = rep("#0088aa8a", nG)
    } else {
      if (length(a) == 1) 
        colors = rep(colors, nG)
    }
    
    # Defining group labels using data, if not given as parameter:
    if (is.null(groupNames))
      groupNames = levels(factor(g))
    
    # Plotting boxplots:
    openGraph(width=7,height=7)
    plot <- ggboxplot(datFrm, x=gName, y=yName, fill = gName,
                      palette = colors, title=plotTitle, 
                      xlab=plotXlab, ylab=plotYlab) + 
      scale_x_discrete(labels=groupNames) + theme(legend.position="none")
    
    # Connecting data points from each subject:
    if (connect) {
      for (sIdx in 1:nS) {
        plot <- plot + geom_line(data=data.frame(x=as.numeric(levels(factor(g))), 
                                         y=as.numeric(datFrm[,yName])[s==sIdx]), 
                         mapping=aes(x=x , y=y), linewidth=0.8, color="gray") + 
          geom_point(color="gray", size=0.7)
      }
    }
    
    print(plot)
  }
  
  if ( !is.null(saveName) )
    saveGraph(file=paste(saveName, "-Boxplot", sep=""), type=graphFileType)
}