
# generate data, run the em algorithm and analyze the results for the DT-HMM model
dthmm.func <- function(run.algorithm, debugging.mode, latex.table) {
  
  # start the timer
  start.time <- Sys.time()
  
  # run the algorithm again
  if (run.algorithm == TRUE) {
    
    # the mean and std of the number of visits
    H.params <- c(10,3)
    
    # set both the true and initial model parameters
    ret.set.model.parameters <- dthmm.set.model.parameters.func(latex.table)
    I <- round(ret.set.model.parameters$I)
    J <- round(ret.set.model.parameters$J)
    L <- round(ret.set.model.parameters$L)
    pi <- ret.set.model.parameters$pi
    Q <- ret.set.model.parameters$Q
    mu  <- ret.set.model.parameters$mu
    eta  <- ret.set.model.parameters$eta
    pi.init <- ret.set.model.parameters$pi.init
    Q.init <- ret.set.model.parameters$Q.init
    mu.init  <- ret.set.model.parameters$mu.init
    eta.init  <- ret.set.model.parameters$eta.init
    
    # run for different number of data samples and Monte-Carlo samples (set M=0 for using the direct forward-backward approach instead of Gibbs sampling)
    N.vals <- c(1000, 10000)
    M.vals <- c(100, 0)
    counter.em.max <- 100
    pi.hat.vals <- array(numeric(), dim=c(length(N.vals), length(M.vals), counter.em.max+1, I))
    Q.hat.vals <- array(numeric(), dim=c(length(N.vals), length(M.vals), counter.em.max+1, I,I,L))
    mu.hat.vals <- array(numeric(), dim=c(length(N.vals), length(M.vals), counter.em.max+1, I))
    eta.hat.vals <- array(numeric(), dim=c(length(N.vals), length(M.vals), counter.em.max+1, J))
    for (N.itr in 1:length(N.vals)) {
      for (M.itr in 1:length(M.vals)) {
        
        # set the number of data samples and Monte-Carlo samples
        N <- N.vals[N.itr]
        M <- M.vals[M.itr]
    
        # the number of burn-in samples
        B <- round(0.1 * M)
        
        # generate synthetic data
        ret.data.generation <- dthmm.data.generation(I, J, L, pi, Q, mu, eta, N, H.params)
        H <- ret.data.generation$H
        z <- ret.data.generation$z
        y <- ret.data.generation$y
        u <- ret.data.generation$u
        
        # run the EM algorithm
        ret.EM.algorithm <- dthmm.EM.algorithm.func(N, M, B, I, J, L, H, y, u, pi.init, Q.init, mu.init, eta.init, counter.em.max, debugging.mode)
        pi.hat.vals[N.itr,M.itr,,] <- ret.EM.algorithm$pi.hat
        Q.hat.vals[N.itr,M.itr,,,,] <- ret.EM.algorithm$Q.hat
        mu.hat.vals[N.itr,M.itr,,] <- ret.EM.algorithm$mu.hat
        eta.hat.vals[N.itr,M.itr,,] <- ret.EM.algorithm$eta.hat
        
      }
    }
    
    # store the raw data (both observed and unobserved) in file (only for the last N and M)
    dthmm_store_raw(I, J, L, N.vals, M.vals, H, z, y, u)
    
    # store the results in file (for all the values of N)
    dthmm_store_results(pi, Q, mu, eta, pi.hat.vals, Q.hat.vals, mu.hat.vals, eta.hat.vals)
    
  } else {
    
    # load the raw data from file (only for the largest N)
    raw.dthmm <- dthmm_load_raw()
    I <- round(raw.dthmm$I)
    J <- round(raw.dthmm$J)
    L <- round(raw.dthmm$L)
    N.vals <- raw.dthmm$N_vals
    M.vals <- raw.dthmm$M_vals
    H <- raw.dthmm$H
    z <- raw.dthmm$z
    y <- raw.dthmm$y
    u <- raw.dthmm$u

    # load the results from file (for all the values of N)
    res.dthmm <- dthmm_load_results()
    pi <- res.dthmm$pi
    Q <- res.dthmm$Q
    mu <- res.dthmm$mu
    eta <- res.dthmm$eta
    pi.hat.vals <- res.dthmm$pi_hat_vals
    Q.hat.vals <- res.dthmm$Q_hat_vals
    mu.hat.vals <- res.dthmm$mu_hat_vals
    eta.hat.vals <- res.dthmm$eta_hat_vals
    
  }
  
  # use Python for visualization of the data
  dthmm_visualize_data(H, z, y, u)
  
  # use Python for visualization of the results
  dthmm_visualize_results(I, J, L, N.vals, M.vals, pi, Q, mu, eta, pi.hat.vals, Q.hat.vals, mu.hat.vals, eta.hat.vals)
  
  # stop the timer and measure the total run time
  end.time <- Sys.time()
  total.time <- end.time - start.time
  
  ret <- list("I" = I, "J" = J, "L" = L, "N.vals" = N.vals, "M.vals" = M.vals, "H" = H, "z" = z, "y" = y, "u" = u, "pi" = pi, "Q" = Q, "mu" = mu, "eta" = eta, "pi.hat.vals" = pi.hat.vals, "Q.hat.vals" = Q.hat.vals, "mu.hat.vals" = mu.hat.vals, "eta.hat.vals" = eta.hat.vals, "total.time" = total.time) 
  
  return(ret)  
  
}

