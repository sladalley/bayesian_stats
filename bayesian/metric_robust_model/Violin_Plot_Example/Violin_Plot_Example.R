
# ======================================================
# Code below is for my violin plot function. I've tried 
# to tidy it up and add relevant comments to explain 
# things as we go. There's also an example of how to
# call the function using example data that should be
# distributed with this file.
#
# If anything isn't working or isn't clear, let me know.
#
#                 ~Dan (daniel.derwent@manchester.ac.uk)
# ======================================================
rm(list=ls())
library(ggplot2)
library(HDInterval)
library(patchwork)

fmt3 <- function(x) formatC(x, format = "fg", digits = 3)
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

# ------------------------------------------------------
# Function takes:
#     df -----------> a data.frame containing the data to be plotted. Has columns:
#                  |-> value ---------> Numeric values to plot
#                  |-> horizon_depth -> The horizon depth
#                  |-> group ---------> If the data is bimodal then this is "A" or "B". If its unimodal then this is "0".
#                  |-> weighting -----> The numeric weighting of this component (if bimodal)
#     save_dir -----> a directory to save the plot in as a string
#     rope_size ----> size of the ROPE to use.
#     fig_width ----> width of the figure
#     fig_height ---> height of the figure
#     plot_title ---> Title of the plot (written at the top)
#     y_axis_label -> Label on the vertical axis
#     x_axis_label -> Label on the horizontal axis
# ------------------------------------------------------
plot_violins = function(
    df,
    save_dir,
    file_name="ViolinPlots-Combined",
    rope_size,
    rope_label = NULL,
    fig_width=2,
    fig_height=4.5,
    plot_title = "Hooray, Violins!",
    y_axis_label = "Goodness",
    x_axis_label = "Violins"
){
  # ------------------------------------------------------
  # Set up some variables to be used later
  # ------------------------------------------------------
  # Rope label information
  ref_label   = "ROPE"
  ref_height  = rope_size
  
  # Colour of the violins
  violin_fill = "skyblue"
  
  # Lengths for the horizontal bars denoting the median and HDI
  bar_length_med <- 0.50
  bar_length_hdi <- 0.30
  
  # Thickness of those bars
  bar_width <-0.5
  
  # Scaler for the text on the labels (e.g., HDI information)
  Text_scaler <-0.9
  
  # ROPE label alignment (left or right)
  #Rope_label_pos = 0 # uncomment for left aligned
  #Rope_label_just = -0.05
  Rope_label_pos = Inf # uncomment for right aligned
  Rope_label_just = 1.05
  
  # Y axis limits
  #y_axis_min <- -80 # uncomment to hard code
  #y_axis_max <- 1000 
  y_axis_min <- min(df$value, -rope_size)
  y_axis_max <- max(df$value, rope_size)
  
  # ------------------------------------------------------
  # Prepare data for plotting
  # ------------------------------------------------------
  # Update labels so that if the dataset is bimodal then the label includes which mode this plot is for

  df$label <- paste0(" ")
  df$label <- factor(df$label, levels = unique(df$label))
  
  # Extract relevant information about each violin into a summary data frame
  summary_df <- do.call(
    rbind,
    lapply(names(split(df$value, df$label)), function(lbl) {
      
      x <- df$value[df$label == lbl]
      
      h <- hdi(x, credMass = 0.95)
      print(h)
      
      data.frame(
        label     = lbl,
        median    = median(x),
        hdi_lower = h[1],
        hdi_upper = h[2],
        top       = max(x),
        bottom    = min(x),
        hdi_label_vertical_offset = 0,
        hdi_label_horizontal_offset = 0,
        median_label_vertical_offset = 0
      )
    })
  )

  # Assign numeric x positions for discrete labels
  summary_df$x <- as.numeric(
    factor(summary_df$label, levels = levels(df$label))
  )
  
  # ------------------------------------------------------
  # Overrides (if any) would be defined here.
  # ------------------------------------------------------
  
  # Overrides can include:
  #     - Moving the HDI label horizontally (hdi_label_horizontal_offset)
  #     - Moving the HDI label vertically (hdi_label_vertical_offset)
  #     - Moving the Median label vertically (median_label_vertical_offset)
  
  # Overrides must be declared by identifying the relevant violin using its label. Examples:
  #summary_df$hdi_label_horizontal_offset[summary_df$label == "k=1"] <- 0.15
  #summary_df$hdi_label_vertical_offset[summary_df$label == "k=3"] <- 100
  #summary_df$median_label_vertical_offset[summary_df$label == "k=4"] <- 700
  
  # ------------------------------------------------------
  # Create the plot
  # ------------------------------------------------------
  rope_text <- if (!is.null(rope_label) && rope_label != "") {
    rope_label
  } else {
    fmt3(rope_size)
  }
  
  label_pos <- paste0(ref_label, " (+", rope_text, ")")
  label_neg <- paste0(ref_label, " (-", rope_text, ")")
  # Add the violins
  p <- ggplot(df, aes(x = label, y = value)) +
    geom_violin(trim = FALSE, fill = violin_fill, colour = "black") +
    
  # Add the median horizontal bars
  geom_segment(
    data = summary_df,
    aes(
      x    = x - bar_length_med,
      xend = x + bar_length_med,
      y    = median,
      yend = median
    ),
    linewidth = bar_width,
    colour = "darkgreen",
    inherit.aes = FALSE
  ) +
    
  # Add the HDI horizontal bars
  geom_segment(
    data = summary_df,
    aes(
      x    = x - bar_length_hdi,
      xend = x + bar_length_hdi,
      y    = hdi_lower,
      yend = hdi_lower
    ),
    linewidth = (0.8/0.9)*bar_width,
    colour = "black",
    inherit.aes = FALSE
  ) +
    
    geom_segment(
      data = summary_df,
      aes(
        x    = x - bar_length_hdi,
        xend = x + bar_length_hdi,
        y    = hdi_upper,
        yend = hdi_upper
      ),
      linewidth = (0.8/0.9)*bar_width,
      colour = "black",
      inherit.aes = FALSE
    ) +
    
    # Add a horizontal line at y=0
    geom_hline(yintercept = 0, linewidth = (0.6/0.9)*bar_width) +
    
    # Add the ROPE lines
    geom_hline(
      yintercept = c(-ref_height, ref_height),
      linetype = "dashed",
      colour = "darkred",
      linewidth = (0.7/0.9)*bar_width
    ) +
    
    # Add the ROPE labels
    annotate(
      "text",
      x = Rope_label_pos,
      y = ref_height,
      label = label_pos,
      colour = "darkred",
      hjust = Rope_label_just,
      vjust = -0.5,
      size = 2.5*Text_scaler
    ) +
    annotate(
      "text",
      x = Rope_label_pos,
      y = -ref_height,
      label = label_neg,
      colour = "darkred",
      hjust = Rope_label_just,
      vjust = 1.5,
      size = 2.5*Text_scaler
    ) +
    
    # Add median annotation
    geom_label(
      data = summary_df,
      aes(
        x = label,
        y = bottom - (0.015 * (top - bottom))+median_label_vertical_offset,
        label = paste0("Md = ", fmt3(median))
      ),
      colour = "darkgreen",
      fill   = "white",
      alpha  = 0.75,
      label.size = 0,     # no border
      hjust = 0.5,
      vjust = 1.0,
      size  = 2.5*Text_scaler,
      inherit.aes = FALSE
    )+
    
    # Add HDI annotations
    geom_label(
      data = summary_df,
      aes(
        x = x+hdi_label_horizontal_offset,
        y = top+hdi_label_vertical_offset,
        label = paste0(
          "HDI: [", fmt3(hdi_lower), ", ", fmt3(hdi_upper),"]"
        )
      ),
      colour = "black",
      fill   = "white",
      alpha  = 0.75,
      label.size = 0,     # no border
      vjust = -0.25,#was -0.75
      size  = 2.5*Text_scaler,
      inherit.aes = FALSE
    )+
    
    # Set y axis properties
    scale_y_continuous(n.breaks = 10) +
    coord_cartesian(ylim = c(y_axis_min, y_axis_max)) + # WW Steps Viva
    
    # Apply classic theme defaults
    theme_classic() +
    
    # Apply axis labels
    labs(
      title = plot_title,
      y = y_axis_label,
      x = x_axis_label
    )+
    
    # Adjust theme settings
    theme(
      # Adjust text sizes
      text = element_text(size = 10 * Text_scaler),
      plot.title = element_text(
        size = 11 * Text_scaler,
        hjust = 0.5,
        margin = margin(b = 4)
      ),
      axis.title.x = element_text(size = 9 * Text_scaler),
      axis.title.y = element_text(size = 9 * Text_scaler),
      axis.text.x  = element_text(size = 8 * Text_scaler, hjust = 0.5),
      axis.text.y  = element_text(size = 8 * Text_scaler),
      
      # Adjust spacing and margins
      panel.spacing.x = unit(0.4, "lines"),
      plot.margin = margin(3, 3, 3, 3),
      
      # Adjust gridlines
      panel.grid.major.y = element_line(colour = "grey65", linewidth = (0.6/0.9)*bar_width),
      panel.grid.minor.y = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    )
  
  # ------------------------------------------------------
  # Save the resulting plot
  # ------------------------------------------------------
  
  ggsave(
    filename = paste0(save_dir, "/", file_name,".pdf"),
    plot = p,
    device = "pdf",
    width = fig_width,
    height = fig_height
  )
  
  # Return the ggplot object, in case you want to modify it any further, but don't print anything to the terminal.
  invisible(p)
}


# ------------------------------------------------------
#     Example of calling this function
# ------------------------------------------------------
# This example loads some data from my thesis experiments that compares
# my planner to the state of the art. My planner's data is bimodal,
# so both modes are plotted as separate violins. 

# It works by loading the data in, calculating the difference in means,
# constructing the data frame, and then calling the function.

# # Load the data in
# #load("Violin_Plot_Example/Wall_Walking-Corrected Contact Changes-RHCP_2-Bimodal-codaSamples.RData")
# load("iros_test_t_rmse/iros_test_t_rmse-Mcmc.Rdata")
# rmse <- codaSamples
# mcmcMatList_rmse = as.matrix(rmse,chains=TRUE)
# 
# 
#  load("iros_test_beta_ce_test/iros_test_beta_ce_test-Mcmc.Rdata")
#  ce <- codaSamples
#  mcmcMatList_ce = as.matrix(ce,chains=TRUE)
# # 
# # 
#  load("iros_test_beta_st_test/iros_test_beta_st_test-Mcmc.Rdata")
#  st <- codaSamples
#  mcmcMatList_st = as.matrix(st,chains=TRUE)
# 
# datFrm = read.csv(file="iros_data_6.csv")
# 
# # Keep only successful runs
# datFrm = subset(datFrm, 
#                 Solved == 1 & 
#                   IsStable == 1 & 
#                   ConstraintsOK == 1 & SimTime <=300)
# 
# 
# 
# 
# y_rmse = as.numeric(datFrm[,"RMSE"])
# 
# y_controleffort = as.numeric(datFrm[,"ControlEffort"])
# y_ce_median = median(y_controleffort)
# res_ce = normalise_data(y_controleffort)
# y_controleffort= res_ce$x_norm
# 
# 
# 
# y_simtime = (as.numeric(datFrm[,"SimTime"]))
# y_st_median = median(y_simtime)
# res_st = normalise_data(y_simtime)
# y_simtime =res_st$x_norm
# 
# gName = "Controller"
# # Getting number of groups from data:
# g = as.numeric(as.factor(datFrm[,gName]))
# nG = max(g) # number of groups
# 
# rope_rmse = 0.05*mean(y_rmse[ g==2 ])
# rope_ce = 0.05*((y_ce_median - res_ce$y_min) / (res_ce$y_max - res_ce$y_min))
# rope_st = 0.05*((y_st_median - res_st$y_min) / (res_st$y_max - res_st$y_min))
# 
# 
# so_qp_rmse = mcmcMatList_rmse[, "mu[2]"]
# fbl_qp_rmse = mcmcMatList_rmse[, "mu[1]"]
# 
# so_qp_ce = (mcmcMatList_ce[,"alpha[2]"]/(mcmcMatList_ce[,"alpha[2]"]+mcmcMatList_ce[,"beta[2]"]))
# fbl_qp_ce = (mcmcMatList_ce[,"alpha[1]"]/(mcmcMatList_ce[,"alpha[1]"]+mcmcMatList_ce[,"beta[1]"]))
# 
# 
# so_qp_st = (mcmcMatList_st[,"alpha[2]"]/(mcmcMatList_st[,"alpha[2]"]+mcmcMatList_st[,"beta[2]"]))
# fbl_qp_st = (mcmcMatList_st[,"alpha[1]"]/(mcmcMatList_st[,"alpha[1]"]+mcmcMatList_st[,"beta[1]"]))


load("paired_rmse/paired_rmse-Mcmc.Rdata")
rmse <- codaSamples
mcmcMat_rmse = as.matrix(rmse,chains=TRUE)


load("paired_controleffort/paired_controleffort-Mcmc.Rdata")
ce <- codaSamples
mcmcMat_ce = as.matrix(ce,chains=TRUE)

load("paired_simtime/paired_simtime-Mcmc.Rdata")
st <- codaSamples
mcmcMat_st = as.matrix(st,chains=TRUE)

load("paired_simtime/paired_simtime_summary_stats.RData")
st_median = 0.1*abs(summary_stats$diff_mad[summary_stats$Controller=="QP"])

load("paired_rmse/paired_rmse_summary_stats.RData")
rmse_median = 0.1*abs(summary_stats$diff_mad[summary_stats$Controller=="QP"])

load("paired_controleffort/paired_controleffort_summary_stats.RData")
ce_median = 0.1*abs(summary_stats$diff_mad[summary_stats$Controller=="QP"])

# Calculate difference in means
mean_diff_rmse = mcmcMat_rmse[,"mu"]
mean_diff_st = mcmcMat_st[,"mu"]
mean_diff_ce = mcmcMat_ce[,"mu"]


#print(mean_diff)


# Define data.frame
violin_df <- data.frame(
  # value = numeric(0),
#  horizon_depth = integer(0),
  group = character(0)
#  weighting = numeric(0)
)

# Add data for the first component
rmse_df <- rbind(
  violin_df,
  data.frame(
    value = mean_diff_rmse,
    group = 0

  )
)

st_df <- rbind(
  violin_df,
  data.frame(
    value = mean_diff_st,
    group = 0
    
  )
)

ce_df <- rbind(
  violin_df,
  data.frame(
    value = mean_diff_ce,
    group = 0
    
  )
)


groupNames=c("FBL-QP", "QP")
x_axis_label = xlab=TeX(sprintf("$\\mu_{diff}"))
x_axis_label_md = xlab=TeX(sprintf("$Md_{diff}"))
# Call the function
p_rmse = plot_violins(rmse_df, # data frame
             getwd(), # save directory (the working directory)
             plot_title = "RMSE",
             "p_rmse", y_axis_label = " ",
             x_axis_label = x_axis_label,
             rope_size = rmse_median, rope_label = "10%") # ROPE 

p_st = plot_violins(st_df, # data frame
                      getwd(), # save directory (the working directory)
                      plot_title = "Simulation Time",
                      "p_st", y_axis_label = " ",
                      x_axis_label = x_axis_label_md,
                      rope_size = st_median, rope_label = "10%") # ROPE 

p_ce= plot_violins(ce_df, # data frame
                      getwd(), # save directory (the working directory)
                      plot_title = "Control Effort",
                      "p_ce", y_axis_label = " ",
                      x_axis_label = x_axis_label,
                      rope_size = ce_median, rope_label = "10%") # ROPE 


combined <- p_rmse | p_ce | p_st
ggsave("combined.pdf", combined, width = 6, height = 4.5)