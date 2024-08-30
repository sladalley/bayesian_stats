This folder contains scripts for Bayesian estimation of parameters. Currently, the main reference for this type of analysis is John Kruschke's work, such as the book _Doing Bayesian Data Analysis_ (Kruschke, 2015).

Descriptions:

**metric_robust_model/**
- metric_robust_functions.R: functions implementing the model described in Chapter 16 of Kruschke (2015). It is for the analysis of metric variable on one or two groups, and uses the robust approach, i.e., it considers a t-distribution instead of a normal distribution.
- metric_robust_example.R: an example on how to use the functions to generate the results *(this is the script you’ll be able to run)*.
- DBDA2E-utilities.R: adapted version of the file with the same name accompanying Kruschke’s book (see [resources](../resources)). It contains a set of methods used for the generation of results.
- data_metric_robust_example.RData: example data set to illustrate the robust model for a metric variable. It is a sample with metric data of 26 subjects and two conditions (groups).
---
References:

- Kruschke, John K. 2015. Doing Bayesian Data Analysis: A Tutorial with R, JAGS, and Stan. Burlington, MA: Academic Press / Elsevier.
