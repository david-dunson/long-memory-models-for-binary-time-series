## R code for implementing the FRAP proposed in Chakraborty, Ovaskainen & Dunson (2022, Annals of Applied Statistics)

This code can be employed to generate posetrior samples for the model parameters in the FRAP model. The model is designed to handle binary 
time series data with long range dependence. The model includes a latent process generating the binary series and has a addtitive 
decomposition in terms of a fixed but unknown trend component and fractional Brownian motion. For more details of the model and the 
corresponding compuatational details refer to the paper: https://projecteuclid.org/journals/annals-of-applied-statistics/volume-16/issue-3/Bayesian-semiparametric-long-memory-models-for-discretized-event-data/10.1214/21-AOAS1546.full

The code depends on the R packages "doParallel" and "tmvtnorm", so make sure you install them from CRAN before running it. To run the Windows version, download the frap_windows.R  and git_demo_run. R file into your current folder in R. Then run the git_demo_run.R file. 
