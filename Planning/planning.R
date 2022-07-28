
# generate data, estimate the parameters, and perform planning
planning.func <- function(run.algorithm, debugging.mode, latex.table) {
  
  if (run.algorithm == TRUE) {
  
    # the mean and std of the number of visits
    H.params <- c(10,3)
    
    # the mean and std of the between visit time distribution
    between.visit.time.params <- c(10,3)
    
    # set both the true and initial model parameters
    ret.set.model.parameters <- planning.set.model.parameters.func(latex.table)
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
    
    # planning parameters
    delta.tau.planning <- 10
    beta.planning <-c (1, 3, 8)
      
    # run for different number of samples
    N.train.vals <- c(10)
    N.test <- 1000
    counter.em.max <- 10
    pi.hat.vals <- array(numeric(), dim=c(length(N.train.vals), counter.em.max+1, I))
    Q.hat.vals <- array(numeric(), dim=c(length(N.train.vals), counter.em.max+1, I,I,L))
    mu.hat.vals <- array(numeric(), dim=c(length(N.train.vals), counter.em.max+1, I))
    eta.hat.vals <- array(numeric(), dim=c(length(N.train.vals), counter.em.max+1, J))
    for (N.train.itr in 1:length(N.train.vals)) {
      
      # set the number of samples
      N.train <- N.train.vals[N.train.itr]
  
      # generate synthetic data for training
      ret.data.generation.train <- cthmm.data.generation(I, J, L, pi, lambda, R, Q, mu, eta, N.train, H.params, between.visit.time.params)
      H.train <- ret.data.generation.train$H
      a.train <- ret.data.generation.train$a
      tau.true.train <- ret.data.generation.train$tau.true
      tau.obs.train <- ret.data.generation.train$tau.obs
      z.true.train <- ret.data.generation.train$z.true
      z.obs.train <- ret.data.generation.train$z.obs
      y.obs.train <- ret.data.generation.train$y.obs
      u.obs.train <- ret.data.generation.train$u.obs
      
      # generate synthetic data for testing
      ret.data.generation.test <- cthmm.data.generation(I, J, L, pi, lambda, R, Q, mu, eta, N.test, H.params, between.visit.time.params)
      H.test <- ret.data.generation.test$H
      a.test <- ret.data.generation.test$a
      tau.true.test <- ret.data.generation.test$tau.true
      tau.obs.test <- ret.data.generation.test$tau.obs
      z.true.test <- ret.data.generation.test$z.true
      z.obs.test <- ret.data.generation.test$z.obs
      y.obs.test <- ret.data.generation.test$y.obs
      u.obs.test <- ret.data.generation.test$u.obs
      
      # run the EM algorithm
      ret.EM.algorithm <- cthmm.EM.algorithm.func(N.train, I, J, L, H.train, a.train, tau.obs.train, y.obs.train, u.obs.train, pi.init, Q.init, mu.init, eta.init, counter.em.max, debugging.mode)
      pi.hat.vals[N.train.itr,,] <- ret.EM.algorithm$pi.hat
      Q.hat.vals[N.train.itr,,,,] <- ret.EM.algorithm$Q.hat
      mu.hat.vals[N.train.itr,,] <- ret.EM.algorithm$mu.hat
      eta.hat.vals[N.train.itr,,] <- ret.EM.algorithm$eta.hat
    
      
      # print(ret.EM.algorithm$pi.hat)
      
      # estimate the posterior of z in the last time period
      pi.hat <- pi.hat.vals[N.train.itr,,]
      Q.hat <- Q.hat.vals[N.train.itr,,,,]
      mu.hat <- mu.hat.vals[N.train.itr,,]
      eta.hat <- eta.hat.vals[N.train.itr,,]
      
      ###############
      ################################ replace the input probabilities with hat
      ################
      ret.posterior.estimation.lastobs <- planning.posterior.estimation.lastobs.func(N.test, I, J, L, H.test, a.test, y.obs.test, pi, mu)
      prob.posterior.lastobs <- ret.posterior.estimation.lastobs$prob.posterior
      ret.posterior.estimation.history <- planning.posterior.estimation.history.func(N.test, I, J, L, H.test, a.test, tau.obs.test, y.obs.test, u.obs.test, pi, Q, mu, eta)
      prob.posterior.history <- ret.posterior.estimation.history$prob.posterior
      
      ###############
      ################################ replace the input probabilities with hat
      ################
      ret.expected.disutility <- planning.expected.disutility.func(N.test, I, J, L, H.test, a.test, z.obs.test, prob.posterior.lastobs, prob.posterior.history, Q, delta.tau.planning, beta.planning)
      expected.disutility.optimal <- ret.expected.disutility$expected.disutility.optimal
      expected.disutility.lastobs <- ret.expected.disutility$expected.disutility.lastobs
      expected.disutility.history <- ret.expected.disutility$expected.disutility.history
      best.choice.optimal <- ret.expected.disutility$best.choice.optimal
      best.choice.lastobs <- ret.expected.disutility$best.choice.lastobs
      best.choice.history <- ret.expected.disutility$best.choice.history
      actual.disutility.optimal <- ret.expected.disutility$actual.disutility.optimal
      actual.disutility.lastobs <- ret.expected.disutility$actual.disutility.lastobs
      actual.disutility.history <- ret.expected.disutility$actual.disutility.history
    }
    
    # store the raw data (both observed and unobserved) in file (only for the largest N)
    # cthmm_store_raw(I, J, L, N.train.vals, H, a, tau.true, tau.obs, z.true, z.obs, y.obs, u.obs)
    
    # store the results in file (for all the values of N)
    # cthmm_store_results(pi, Q, mu, eta, pi.hat.vals, Q.hat.vals, mu.hat.vals, eta.hat.vals)
  
  } else {
    
    # load the raw data from file (only for the largest N)
    raw.cthmm <- cthmm_load_raw()
    I <- round(raw.cthmm$I)
    J <- round(raw.cthmm$J)
    L <- round(raw.cthmm$L)
    N.train.vals <- raw.cthmm$N_vals
    H <- raw.cthmm$H
    a <- raw.cthmm$a
    tau.true <- raw.cthmm$tau_true
    tau.obs <- raw.cthmm$tau_obs
    z.true <- raw.cthmm$z_true
    z.obs <- raw.cthmm$z_obs
    y.obs <- raw.cthmm$y_obs
    u.obs <- raw.cthmm$u_obs

    # load the results from file (for all the values of N)
    res.cthmm <- cthmm_load_results()
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
  # cthmm_visualize_data(H, a, tau.true, tau.obs, z.true, z.obs, y.obs, u.obs)
  
  # use Python for visualization of the results
  # cthmm_visualize_results(I, J, L, N.train.vals, pi, Q, mu, eta, pi.hat.vals, Q.hat.vals, mu.hat.vals, eta.hat.vals)
  
  ret <- list("I" = I, "J" = J, "L" = L, "N.train.vals" = N.train.vals, "N.test" = N.test, "H.train" = H.train, "H.test" = H.test, "a.train" = a.train, "a.test" = a.test, "tau.true.train" = tau.true.train, "tau.true.test" = tau.true.test, "tau.obs.train" = tau.obs.train, "tau.obs.test" = tau.obs.test, "z.true.train" = z.true.train, "z.true.test" = z.true.test, "z.obs.train" = z.obs.train, "z.obs.test" = z.obs.test, "y.obs.train" = y.obs.train, "y.obs.test" = y.obs.test, "u.obs.train" = u.obs.train, "u.obs.test" = u.obs.test, "pi" = pi, "Q" = Q, "mu" = mu, "eta" = eta, "pi.hat.vals" = pi.hat.vals, "Q.hat.vals" = Q.hat.vals, "mu.hat.vals" = mu.hat.vals, "eta.hat.vals" = eta.hat.vals, "prob.posterior.lastobs" = prob.posterior.lastobs, "prob.posterior.history" = prob.posterior.history, "expected.disutility.optimal" = expected.disutility.optimal, "expected.disutility.lastobs" = expected.disutility.lastobs, "expected.disutility.history" = expected.disutility.history, "best.choice.optimal" = best.choice.optimal, "best.choice.lastobs" = best.choice.lastobs, "best.choice.history" = best.choice.history, "actual.disutility.optimal" = actual.disutility.optimal, "actual.disutility.lastobs" = actual.disutility.lastobs, "actual.disutility.history" = actual.disutility.history) 
  
  return(ret)  
  
}

