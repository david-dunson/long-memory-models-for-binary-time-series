## R code for implementing the FRAP proposed in Chakraborty, Ovaskainen & Dunson (2020+)

This code can be employed to generate posetrior samples for the model parameters in FRAP model. The model is designed to handle binary 
time series data with long range dependence. The model proposes a latent process generating the binary series and has a addtitive 
decomposition in terms of a fixed but unknown trend component and fractional Brownian motion. For more details of the model and the 
corresponding compuatational details please visit my [website](https://antik015.github.io/publications/).

The code depends on the R packages "doParallel" and "tmvtnorm", so make sure you install them from CRAN before running it.
