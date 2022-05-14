
# generate data, run the em algorithm and analyze the results for the CTHMM model
cthmm.func <- function(run.algorithm, debugging.mode, latex.table) {
  
  if (run.algorithm == TRUE) {
  
    # the mean and std of the number of visits
    H.params <- c(10,3)
    
    # the mean and std of the between visit time distribution
    between.visit.time.params <- c(10,3)
    
    # set both the true and initial model parameters
    ret.set.model.parameters <- set.model.parameters.func(latex.table)
    I <- round(ret.set.model.parameters$I)
    J <- round(ret.set.model.parameters$J)
    L <- round(ret.set.model.parameters$L)
    pi <- ret.set.model.parameters$pi
    lambda <- ret.set.model.parameters$lambda
    R <- ret.set.model.parameters$R
    Q <- ret.set.model.parameters$Q
    mu  <- ret.set.model.parameters$mu
    eta  <- ret.set.model.parameters$eta
    pi.init <- ret.set.model.parameters$pi.init
    Q.init <- ret.set.model.parameters$Q.init
    mu.init  <- ret.set.model.parameters$mu.init
    eta.init  <- ret.set.model.parameters$eta.init
    
    # run for different number of samples
    N.vals <- c(10, 100, 1000, 10000)
    counter.em.max <- 100
    pi.hat.vals <- array(numeric(), dim=c(length(N.vals), counter.em.max+1, I))
    Q.hat.vals <- array(numeric(), dim=c(length(N.vals), counter.em.max+1, I,I,L))
    mu.hat.vals <- array(numeric(), dim=c(length(N.vals), counter.em.max+1, I))
    eta.hat.vals <- array(numeric(), dim=c(length(N.vals), counter.em.max+1, J))
    for (N.itr in 1:length(N.vals)) {
      
      # set the number of samples
      N <- N.vals[N.itr]
  
      # generate synthetic data
      ret.data.generation <- data.generation(I, J, L, pi, lambda, R, Q, mu, eta, N, H.params, between.visit.time.params)
      H <- ret.data.generation$H
      a <- ret.data.generation$a
      tau.true <- ret.data.generation$tau.true
      tau.obs <- ret.data.generation$tau.obs
      z.true <- ret.data.generation$z.true
      z.obs <- ret.data.generation$z.obs
      y.obs <- ret.data.generation$y.obs
      u.obs <- ret.data.generation$u.obs
      
      # run the EM algorithm
      ret.EM.algorithm <- EM.algorithm.func(N, I, J, L, H, a, tau.obs, y.obs, u.obs, pi.init, Q.init, mu.init, eta.init, counter.em.max, debugging.mode)
      pi.hat.vals[N.itr,,] <- ret.EM.algorithm$pi.hat
      Q.hat.vals[N.itr,,,,] <- ret.EM.algorithm$Q.hat
      mu.hat.vals[N.itr,,] <- ret.EM.algorithm$mu.hat
      eta.hat.vals[N.itr,,] <- ret.EM.algorithm$eta.hat
    
    }
    
    # store the raw data (both observed and unobserved) in file (only for the largest N)
    store_raw(I, J, L, N.vals, H, a, tau.true, tau.obs, z.true, z.obs, y.obs, u.obs)
    
    # store the results in file (for all the values of N)
    store_results(pi, Q, mu, eta, pi.hat.vals, Q.hat.vals, mu.hat.vals, eta.hat.vals)
  
  } else {
    
    # load the raw data from file (only for the largest N)
    raw.cthmm <- load_raw()
    I <- round(raw.cthmm$I)
    J <- round(raw.cthmm$J)
    L <- round(raw.cthmm$L)
    N.vals <- raw.cthmm$N_vals
    H <- raw.cthmm$H
    a <- raw.cthmm$a
    tau.true <- raw.cthmm$tau_true
    tau.obs <- raw.cthmm$tau_obs
    z.true <- raw.cthmm$z_true
    z.obs <- raw.cthmm$z_obs
    y.obs <- raw.cthmm$y_obs
    u.obs <- raw.cthmm$u_obs

    # load the results from file (for all the values of N)
    res.cthmm <- load_results()
    pi <- res.cthmm$pi
    Q <- res.cthmm$Q
    mu <- res.cthmm$mu
    eta <- res.cthmm$eta
    pi.hat.vals <- res.cthmm$pi_hat_vals
    Q.hat.vals <- res.cthmm$Q_hat_vals
    mu.hat.vals <- res.cthmm$mu_hat_vals
    eta.hat.vals <- res.cthmm$eta_hat_vals
    
  }
  
  # use Python for visualization of the data
  visualize_data(I, J, L, N.vals, H, a, tau.true, tau.obs, z.true, z.obs, y.obs, u.obs)
  
  # use Python for visualization of the results
  visualize_results(pi, Q, mu, eta, pi.hat.vals, Q.hat.vals, mu.hat.vals, eta.hat.vals)
  
  ret <- list("I" = I, "J" = J, "L" = L, "N.vals" = N.vals, "H" = H, "a" = a, "tau.true" = tau.true, "tau.obs" = tau.obs, "z.true" = z.true, "z.obs" = z.obs, "y.obs" = y.obs, "u.obs" = u.obs, "pi" = pi, "Q" = Q, "mu" = mu, "eta" = eta, "pi.hat.vals" = pi.hat.vals, "Q.hat.vals" = Q.hat.vals, "mu.hat.vals" = mu.hat.vals, "eta.hat.vals" = eta.hat.vals) 
  
  return(ret)  
  
}

