frap_windows = function(Z, nmcmc, burnin, thin, a_tau, b_tau, prop_sd, prop_sd_H, nugget, smooth.posterior = F)
{
  ## Function for fitting the fractional probit model (FRAP) ##
  ### Inputs ###
  
  ## Z = matrix of binary event indicators with n rows and n_replicates columns. 
  ##     Each column represent a day of recording of n minutes.
  
  ## nmcmc = Number of MCMC samples to draw from the posterior
  ## burnin = Number of burnin samples to throw away from nmcmc
  ## thin = thinning parameter of the chain
  
  ## a_tau = shape parameter of the inverse-Gamma prior on the precision parameter \tau.
  ## b_tau = rate parameter of the inverse-Gamma prior on the precision parameter \tau.
  
  ## n_cores = Number of cores on which the program would run. Empirically fast results are obtained for 4 cores.
  
  ## prop_sd = Proposal standard deviation of the random walk Metropolis Hastings for the Gaussian process
  ##           smoothness parameters.
  
  ## prop_sd_H = Proposal standard deviation of the random walk Metropolis Hastings for the Hurst coefficient.
  
  ## nugget = Value of the nugget parameter to be added to the diagonal of the prior covariance matrix.
  
  ## smooth.posterior = A logical variable indicating whether the posterior samples of the marginal 
  ##                    posterior are to be smoothed or not. Default is set to False.
  
  ### Outputs ###
  
  ## xtest = Vector of starting points of test intervals. These are 0.5, 1.5, ...
  ## m_probs = Posterior samples of the marginal probability of observing at least one event on test intervals 
  ##           (0.5, 1.5), (1.5, 2.5) ... A matrix of n rows and (nmcmc - burnin)/thin columns.
  
  ## m_probs_train = Posterior samples of the marginal probability of observing at least one event on
  ##                 training intervals (0, 1), (1, 2), .... Dimension is same s m_probs.
  
  ## Hout = Vector of posterior samples of the Hurst coefficient of length (nmcmc - burnin)/thin.
  
  ## tauout = Vector of posterior samples of the precision parameter \tau of length (nmcmc - burnin)/thin.
  
  ## phiout = Vector of posterior samples of the lengthscale parameter of squared exponential kernel.
  ## sigmaout = Vector of posterior samples of the amplitude parameter of squared exponential kernel.
  
  
  
  ## Define internal functions ##
  
  ## fGN covariance kernel##
  cov_form = function(k, H)
  {
    g = 0.5*((k+1)^(2*H) - 2*k^(2*H) + (k-1)^(2*H))
    return(g)
  }
  
  ## fGN covariance matrix
  Sigma = function(hurst, p)
  {
    H = hurst
    x = c(1, cov_form((1:(p-1)), H))
    S = toeplitz(x)
    return(S)
  }
  
  ## Squared exponential covariance kernel
  sqe = function(d, phi = 1)
  {
    ifelse(d>0, exp(-(0.5*d^2)/(phi)^2), 1.0)
  }
  
  ## Matern covariance
  
  matern = function(d , phi, nu = 1.5)
  {
    ifelse(d>0, (sqrt(2*nu)*d*phi)^nu/(2^(nu-1)*gamma(nu))*besselK(x=d*phi*sqrt(2*nu), nu=nu), 1.0)
  }
  
  
  ## Replicated version of rtmvnorm
  #my_rtmvnorm = function(j, mean, sigma, lower, upper)
  #{
  #  x = rtmvnorm(1, lower = lower[,j], upper = upper[,j], mean = as.vector(mean), sigma = sigma, algorithm = "gibbs", burn.in.samples = 50)
  #  return(x)
  #}
  
  ## Data dimensions
  ntrain = nrow(Z)
  n_replicates = ncol(Z)
  xtrain = 1:ntrain
  xtest = seq(0.5, ntrain+0.5, by = 1)
  ntest = length(xtest)
  n = ntrain + ntest
  
  ## Difference matrices 
  A = matrix(0, ntrain, ntrain)
  A[1,] = c(1, rep(0, ntrain - 1))
  for(i in 2:ntrain)
  {
    A[i,i] = 1
    A[i,i-1] = -1
  }
  A_star = solve(A)
  
  ## Distance matrix ##
  D = as.matrix(dist(xtrain,diag=TRUE,upper=TRUE))
  D_test = as.matrix(dist(xtest, diag = T, upper = T))
  D_all = as.matrix(dist(c(xtrain, xtest), diag = T, upper = T))
  
  ## Initial values ##
  
  # fGN parameters
  logit_H = 0
  Sigma_H = Sigma(exp(logit_H)/(1+exp(logit_H)), ntrain)
  Sigma_H_inv = chol2inv(chol(Sigma_H))
  
  # Noise variance
  log_tau = 0
  tau = exp(log_tau)
  
  # Trend function differences
  g = rep(0, ntrain)
  
  # GP covarince
  
  log_phi = log(1)
  phi = exp(log_phi)
  log_sigma = 0
  sigma = exp(log_sigma)
  
  cov_mat_all = tau^2*exp(log_sigma)^2*sqe(D_all, exp(log_phi)) + nugget*diag(n)
  cov_mat = cov_mat_all[1:ntrain, 1:ntrain]
  cov_mat_test = cov_mat_all[((ntrain+1):n),((ntrain+1):n)]
  cov_mat_test_train = cov_mat_all[(ntrain+1):n, 1:ntrain]
  
  cov_mat_inv = chol2inv(chol(cov_mat)) 
  
  cov_mat_diff = A%*%cov_mat%*%t(A)
  cov_mat_diff_inv = chol2inv(chol(cov_mat_diff))
  
  kernel_params = c(log_sigma, log_phi)
  kernel_params_store = matrix(kernel_params, nrow = 1, ncol = 2)
  
  
  ## Constraint matrices
  
  L = matrix(0, ntrain, n_replicates)
  U = matrix(0, ntrain, n_replicates)
  for(i in 1:ntrain)
  {
    for(j in 1:n_replicates)
    {
      if(Z[i,j] == 1)
      {
        L[i,j] = 0
        U[i,j] = Inf
      }
      else if(Z[i,j] == 0)
      {
        L[i,j] = -Inf
        U[i,j] = 0
      }
    }
  }
  
  ## MCMC storage
  effsamp = (nmcmc - burnin)/thin
  m_probs = matrix(0, ntest-1, effsamp)
  m_probs_train = matrix(0, ntrain, effsamp)
  tauout = rep(0, effsamp)
  sigmaout = rep(0, effsamp)
  phiout = rep(0, effsamp)
  Hout = rep(0, effsamp)
  
  ## MH_paramaters
  MH_kernel = 0
  MH_Hurst = 0
  
  ## Latent variable matrix
  what = matrix(0, ntrain, n_replicates)
  
  ## Parallelization 
  #cl = makeForkCluster(n_cores)
  #registerDoParallel(cl)
  
  ## Start MCMC ##
  
  for(ii in 1:nmcmc)
  {
    ## Update latent variables
    
    #what = parSapply(cl, 1:n_replicates, my_rtmvnorm, lower = L, upper = U, mean = as.vector(g), sigma = tau^2*Sigma_H)
    for(j in 1:n_replicates)
    {
      what[,j] = rtmvnorm(1, lower = L[,j], upper = U[,j], mean = as.vector(g), sigma = tau^2*Sigma_H, algorithm = "gibbs", burn.in.samples = 50)
    }
    
    ## Update trend differences
    
    wbar = rowMeans(what, na.rm = T)
    total_cov = (n_replicates/tau^2)*Sigma_H_inv + cov_mat_diff_inv
    chol_total_cov = chol(total_cov)
    inv_total_cov = chol2inv(chol_total_cov)
    g = (n_replicates/tau^2)*inv_total_cov%*%Sigma_H_inv%*%wbar + solve(chol_total_cov, rnorm(ntrain))
    if(sum(is.na(g)) > 0)
    {
      g =A%*%ftrain
    }
    ftrain = A_star%*%g
    
    ## Update trend at test points
    
    y = A_star%*%wbar
    total_cov1 = (tau^2/n_replicates)*A_star%*%Sigma_H%*%t(A_star) + cov_mat
    inv_total_cov1 = chol2inv(chol(total_cov1))
    test_cov = cov_mat_test - cov_mat_test_train%*%inv_total_cov1%*%t(cov_mat_test_train) 
    test_cov = (test_cov + t(test_cov))/2
    test_mean = cov_mat_test_train%*%inv_total_cov1%*%y
    ftest = t(rmvnorm(1, test_mean, test_cov))
    
    ## Update Hurst parameter
    
    if(ii%%50 == 0 && MH_Hurst/ii > 0.3)
    {
      prop_sd_H = prop_sd_H*exp(min(0.01, ii^(-0.5)))
    }else if(ii%%50 == 0 && MH_Hurst/ii < 0.3){
      prop_sd_H = prop_sd_H*exp(-min(0.01, ii^(-0.5)))
    }
    
    logit_H_cand = min(rnorm(1, logit_H, prop_sd_H), 10) # for numerical stability
    Sigma_H_cand = Sigma(exp(logit_H_cand)/(1+exp(logit_H_cand)), ntrain)
    Sigma_H_cand_inv = chol2inv(chol(Sigma_H_cand))
    
    r1 = sum(dmvnorm(t(what), as.vector(g), tau^2*Sigma_H_cand, log = T), na.rm = T) + dnorm(logit_H_cand, log = T) + dnorm(logit_H, logit_H_cand, prop_sd_H, log = T)
    r2 = sum(dmvnorm(t(what), as.vector(g), tau^2*Sigma_H, log = T), na.rm = T) + dnorm(logit_H, log = T) + dnorm(logit_H_cand, logit_H, prop_sd_H, log = T)
    
    if(r1 - r2 > log(runif(1)))
    {
      logit_H = logit_H_cand
      Sigma_H = Sigma_H_cand
      Sigma_H_inv = Sigma_H_cand_inv
      MH_Hurst = MH_Hurst + 1
    }
    
    ## Update tau ##
    
    G = replicate(n_replicates, g, simplify = "matrix")
    S = 0.5*sum(diag(t(what- G)%*%Sigma_H_inv%*%(what - G)),na.rm = T) + (0.5/tau^2)*(t(g)%*%cov_mat_diff_inv%*%g)
    tau = max(sqrt(1/rgamma(1, ntrain*(n_replicates + 1)/2 + a_tau, scale = 1/(S + b_tau))), 0.1)
    #tau = sqrt(1/rgamma(1, ntrain*(n_replicates + 1)/2 + a_tau, scale = 1/(S + b_tau)))
    # ## Update GP kernel parameters 
    
    if(ii%%50 == 0 && MH_kernel/ii > 0.3)
    {
      prop_sd = prop_sd*exp(min(0.01, ii^(-0.5)))
    }else if(ii%%50 == 0 && MH_kernel/ii < 0.3){
      prop_sd = prop_sd*exp(-min(0.01, ii^(-0.5)))
    }
    prop_cov2 = prop_sd^2*diag(length(kernel_params))
    kernel_params_cand = rmvnorm(1, kernel_params, prop_cov2)
    
    d11 = sum(dnorm(kernel_params_cand, log = T)) + dmvnorm(kernel_params, kernel_params_cand, prop_cov2, log = T)
    d22 = sum(dnorm(kernel_params, log = T)) + dmvnorm(kernel_params_cand, kernel_params, prop_cov2, log = T)
    sigma_cand = exp(kernel_params_cand[1])
    phi_cand = exp(kernel_params_cand[2])
    cov_mat_cand = tau^2*sigma_cand^2*sqe(D, phi_cand) + nugget*diag(ntrain)
    
    l11 = dmvnorm(t(ftrain), rep(0, ntrain), cov_mat_cand, log = T) + d11
    l22 = dmvnorm(t(ftrain), rep(0, ntrain), cov_mat, log = T) + d22
    
    
    if(is.na(l11) + is.na(l22) > 0)
    {
      kernel_params = kernel_params
      sigma = sigma
      phi = phi
    }else if(is.na(l11) + is.na(l22) == 0){
      if(l11 - l22 > log(runif(1)))
      {
        kernel_params = kernel_params_cand
        sigma = sigma_cand
        phi = phi_cand
        cov_mat_all = tau^2*sigma^2*sqe(D_all, phi) + nugget*diag(n)
        cov_mat = cov_mat_all[1:ntrain, 1:ntrain]
        cov_mat_test = cov_mat_all[((ntrain+1):n),((ntrain+1):n)]
        cov_mat_test_train = cov_mat_all[(ntrain+1):n, 1:ntrain]
        cov_mat_inv = chol2inv(chol(cov_mat))
        cov_mat = cov_mat_cand
        cov_mat_diff = A%*%cov_mat%*%t(A)
        cov_mat_diff_inv = chol2inv(chol(cov_mat))
        MH_kernel = MH_kernel + 1
      }
    }
    kernel_params_store = rbind(kernel_params_store, kernel_params)
    
    ## Print ii
    
    # ## Checks 
    if(ii%%50 == 0)
    {
      print(paste("Number of MCMC iteration completed = ", ii))
      print(paste("Accetptance probability for Hurst coefficient = ", round(MH_Hurst/ii, 2)))
      print(paste("Accetptance probability for GP parameters = ", round(MH_kernel/ii, 2)))
    }
    
    ## Store values 
    if(ii > burnin && ii%%thin == 0)
    {
      if(smooth.posterior == T)
      {
        fsmooth = loess.smooth(1:(ntest-1), diff(ftest), evaluation = ntest - 1, span = 0.25)$y
        fsmooth_train = loess.smooth(1:ntrain, diff(c(0, ftrain)), evaluation = ntrain, span = 0.25)$y
      } else if(smooth.posterior == F)
      {
        fsmooth = diff(ftest)
        fsmooth_train = diff(c(0, ftrain))
      }
      m_probs[,(ii - burnin)/thin] = pnorm(fsmooth/tau)
      m_probs_train[,(ii - burnin)/thin] = pnorm(fsmooth_train/tau)
      Hout[(ii - burnin)/thin] = exp(logit_H)/(1+exp(logit_H))
      tauout[(ii - burnin)/thin] = tau
      phiout[(ii - burnin)/thin] = phi
      sigmaout[(ii - burnin)/thin] = sigma
    }
  }
  #stopCluster(cl)
  result = list("xtest" = xtest, "m_probs" = m_probs, "m_probs_train" = m_probs_train, "Hout" = Hout, "tauout" = tauout, "phiout" = phiout, "sigmaout" = sigmaout)
  return(result)
}