
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
      
      # use Python for visualization of the data
      # visualize_data(H, a, tau.true, tau.obs, z.true, z.obs, y.obs, u.obs)
      
      # run the EM algorithm
      ret.EM.algorithm <- EM.algorithm.func(N, I, J, L, H, a, tau.obs, y.obs, u.obs, pi.init, Q.init, mu.init, eta.init, counter.em.max, debugging.mode)
      pi.hat.vals[N.itr,,] <- ret.EM.algorithm$pi.hat
      Q.hat.vals[N.itr,,,,] <- ret.EM.algorithm$Q.hat
      mu.hat.vals[N.itr,,] <- ret.EM.algorithm$mu.hat
      eta.hat.vals[N.itr,,] <- ret.EM.algorithm$eta.hat
    
    }
    
    # store the results in file
    store_results(I, J, L, N.vals, pi, Q, mu, eta, pi.hat.vals, Q.hat.vals, mu.hat.vals, eta.hat.vals)
  
  } else {
    
    # load the results from file
    res.cthmm <- load_results()
    I <- round(ret.cthmm$I)
    J <- round(ret.cthmm$J)
    L <- round(ret.cthmm$L)
    N.vals <- ret.cthmm$N.vals
    pi <- ret.cthmm$pi
    Q <- ret.cthmm$Q
    mu <- ret.cthmm$mu
    eta <- ret.cthmm$eta
    pi.hat.vals <- ret.cthmm$pi.hat.vals
    Q.hat.vals <- ret.cthmm$Q.hat.vals
    mu.hat.vals <- ret.cthmm$mu.hat.vals
    eta.hat.vals <- ret.cthmm$eta.hat.vals
    
  }
  
  # use Python for visualization of the results
  visualize_results(I, J, L, N.vals, pi, Q, mu, eta, pi.hat.vals, Q.hat.vals, mu.hat.vals, eta.hat.vals)
  
  ret <- list("I" = I, "J" = J, "L" = L, "N.vals" = N.vals, "pi" = pi, "Q" = Q, "mu" = mu, "eta" = eta, "pi.hat.vals" = pi.hat.vals, "Q.hat.vals" = Q.hat.vals, "mu.hat.vals" = mu.hat.vals, "eta.hat.vals" = eta.hat.vals) 
  
  return(ret)  
  
}

