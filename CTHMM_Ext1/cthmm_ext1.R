
# generate data, run the em algorithm and analyze the results for the CTHMM Ext1 model
cthmm.ext1.func <- function(run.algorithm, debugging.mode, latex.table) {
  
  if (run.algorithm == TRUE) {
  
    # the mean and std of the number of visits
    H.params <- c(10,3)
    
    # the mean and std of the between visit time distribution
    between.visit.time.params <- c(10,3)
    
    # initial probability of accurate examination as well as the mean and slope of its Sigmoid function
    accurate.examination.params <- c(0.1,8,1/5,2.7,1/5)
    
    # set both the true and initial model parameters
    ret.set.model.parameters <- cthmm.ext1.set.model.parameters.func(latex.table)
    I <- round(ret.set.model.parameters$I)
    J <- round(ret.set.model.parameters$J)
    L <- round(ret.set.model.parameters$L)
    pi <- ret.set.model.parameters$pi
    lambda <- ret.set.model.parameters$lambda
    R <- ret.set.model.parameters$R
    Q <- ret.set.model.parameters$Q
    mu  <- ret.set.model.parameters$mu
    eta  <- ret.set.model.parameters$eta
    eta.prime  <- ret.set.model.parameters$eta.prime
    pi.init <- ret.set.model.parameters$pi.init
    Q.init <- ret.set.model.parameters$Q.init
    mu.init  <- ret.set.model.parameters$mu.init
    eta.init  <- ret.set.model.parameters$eta.init
    eta.prime.init  <- ret.set.model.parameters$eta.prime.init
    
    # run for different number of samples
    N.vals <- c(50)
    counter.em.max <- 10
    pi.hat.vals <- array(numeric(), dim=c(length(N.vals), counter.em.max+1, I))
    Q.hat.vals <- array(numeric(), dim=c(length(N.vals), counter.em.max+1, I,I,L))
    mu.hat.vals <- array(numeric(), dim=c(length(N.vals), counter.em.max+1, I))
    eta.hat.vals <- array(numeric(), dim=c(length(N.vals), counter.em.max+1, J))
    eta.prime.hat.vals <- array(numeric(), dim=c(length(N.vals), counter.em.max+1, I))
    for (N.itr in 1:length(N.vals)) {
      
      # set the number of samples
      N <- N.vals[N.itr]
  
      # generate synthetic data
      ret.data.generation <- cthmm.ext1.data.generation(I, J, L, pi, lambda, R, Q, mu, eta, eta.prime, N, H.params, between.visit.time.params, accurate.examination.params)
      H <- ret.data.generation$H
      a <- ret.data.generation$a
      O <- ret.data.generation$O
      tau.true <- ret.data.generation$tau.true
      tau.obs <- ret.data.generation$tau.obs
      z.true <- ret.data.generation$z.true
      z.obs <- ret.data.generation$z.obs
      z.acc <- ret.data.generation$z.acc
      y.obs <- ret.data.generation$y.obs
      u.obs <- ret.data.generation$u.obs
      
      # run the EM algorithm
      ret.EM.algorithm <- cthmm.ext1.EM.algorithm.func(N, I, J, L, H, a, O, tau.obs, z.acc, y.obs, u.obs, pi.init, Q.init, mu.init, eta.init, eta.prime.init, counter.em.max, debugging.mode)
      pi.hat.vals[N.itr,,] <- ret.EM.algorithm$pi.hat
      Q.hat.vals[N.itr,,,,] <- ret.EM.algorithm$Q.hat
      mu.hat.vals[N.itr,,] <- ret.EM.algorithm$mu.hat
      eta.hat.vals[N.itr,,] <- ret.EM.algorithm$eta.hat
      eta.prime.hat.vals[N.itr,,] <- ret.EM.algorithm$eta.prime.hat
    
    }
    
    # store the raw data (both observed and unobserved) in file (only for the largest N)
    cthmm_ext1_store_raw(I, J, L, N.vals, H, a, O, tau.true, tau.obs, z.true, z.obs, z.acc, y.obs, u.obs)
    
    # store the results in file (for all the values of N)
    cthmm_ext1_store_results(pi, Q, mu, eta, eta.prime, pi.hat.vals, Q.hat.vals, mu.hat.vals, eta.hat.vals, eta.prime.hat.vals)
  
  } else {
    
    # load the raw data from file (only for the largest N)
    raw.cthmm.ext1 <- cthmm_ext1_load_raw()
    I <- round(raw.cthmm.ext1$I)
    J <- round(raw.cthmm.ext1$J)
    L <- round(raw.cthmm.ext1$L)
    N.vals <- raw.cthmm.ext1$N_vals
    H <- raw.cthmm.ext1$H
    a <- raw.cthmm.ext1$a
    O <- raw.cthmm.ext1$O
    tau.true <- raw.cthmm.ext1$tau_true
    tau.obs <- raw.cthmm.ext1$tau_obs
    z.true <- raw.cthmm.ext1$z_true
    z.obs <- raw.cthmm.ext1$z_obs
    z.acc <- raw.cthmm.ext1$z_acc
    y.obs <- raw.cthmm.ext1$y_obs
    u.obs <- raw.cthmm.ext1$u_obs

    # load the results from file (for all the values of N)
    res.cthmm.ext1 <- cthmm_ext1_load_results()
    pi <- res.cthmm.ext1$pi
    Q <- res.cthmm.ext1$Q
    mu <- res.cthmm.ext1$mu
    eta <- res.cthmm.ext1$eta
    eta.prime <- res.cthmm.ext1$eta_prime
    pi.hat.vals <- res.cthmm.ext1$pi_hat_vals
    Q.hat.vals <- res.cthmm.ext1$Q_hat_vals
    mu.hat.vals <- res.cthmm.ext1$mu_hat_vals
    eta.hat.vals <- res.cthmm.ext1$eta_hat_vals
    eta.prime.hat.vals <- res.cthmm.ext1$eta_prime_hat_vals
    
  }
  
  # use Python for visualization of the data
  cthmm_visualize_data(H, a, O, tau.true, tau.obs, z.true, z.obs, z.acc, y.obs, u.obs)
  
  # use Python for visualization of the results
  cthmm_visualize_results(I, J, L, N.vals, pi, Q, mu, eta, eta.prime, pi.hat.vals, Q.hat.vals, mu.hat.vals, eta.hat.vals, eta.prime.hat.vals)
  
  ret <- list("I" = I, "J" = J, "L" = L, "N.vals" = N.vals, "H" = H, "a" = a, "O" = O, "tau.true" = tau.true, "tau.obs" = tau.obs, "z.true" = z.true, "z.obs" = z.obs, "z.acc" = z.acc, "y.obs" = y.obs, "u.obs" = u.obs, "pi" = pi, "Q" = Q, "mu" = mu, "eta" = eta, "eta.prime" = eta.prime, "pi.hat.vals" = pi.hat.vals, "Q.hat.vals" = Q.hat.vals, "mu.hat.vals" = mu.hat.vals, "eta.hat.vals" = eta.hat.vals, "eta.prime.hat.vals" = eta.prime.hat.vals) 
  
  return(ret)  
  
}

