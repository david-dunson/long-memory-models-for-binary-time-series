## Git demo run ##

## Trend function 

f1 = function(x) 
{
  r = sin(4*pi*x/90)
  return(r)
}

## Dependencies

cov_form = function(k, H)
{
  g = 0.5*((k+1)^(2*H) - 2*k^(2*H) + (k-1)^(2*H))
  return(g)
}

Sigma = function(hurst, p)
{
  H = hurst
  x = c(1, cov_form((1:(p-1)), H))
  S = toeplitz(x)
  return(S)
}


## Generate data

H = 0.75
tau = 0.1
n_replicates = 15
ntrain = 90
xtrain = 1:ntrain

Sigma_H_true = Sigma(hurst = H, p = ntrain)
chol_Sigma_H_true = chol(Sigma_H_true)
f_true = c(0, f1(xtrain))

Z = matrix(0, ntrain, n_replicates) ## Matrix of binary data - each column represents one day recording

tau = 0.1
wtrue = matrix(0, ntrain, n_replicates)
for(l in 1:n_replicates)
{
  wtrue[,l] = diff(f_true)  + tau*chol_Sigma_H_true%*%rnorm(ntrain)
  Z[wtrue[,l]>0, l] = 1
}

source("frap_windows.R") ## Download the "frap_windows.R" file, save it in the current folder then run this line

## MCMC implementation (short)
nmcmc = 500
burnin = 100
thin = 1
a_tau = 0.01
b_tau = 0.01
nugget = 0.001
prop_sd = 0.1
prop_sd_H = 0.2
res = frap_windows(Z, nmcmc, burnin, thin, a_tau, b_tau, prop_sd, prop_sd_H, nugget)
Hout = res$Hout ## Posterior samples of Hurst exponent
hist(Hout, main = "Histogram of posterior samples of Hurst exponent")
m_probs_out = res$m_probs_train
m_probs_intervals = t(apply(m_probs_out, 1, quantile, c(0.025, 0.975)))
m_probs_mean = rowMeans(m_probs_out)

plot(xtrain, m_probs_mean, type = "l", ylim = c(0,1), main = "Estimated marginal probabilities over time", xlab = "Time", ylab = "Marginal probability")
lines(m_probs_intervals[,1], col = "red")
lines(m_probs_intervals[,2], col = "red")
