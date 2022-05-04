
# generate data, run the em algorithm and analyze the results for the CTHMM model
cthmm.func <- function(debugging.mode) {
  
  # the number of observations
  N <- 100
  
  # the mean and std of the number of visits
  H.params <- c(10,3)
  
  # the mean and std of the between visit time distribution
  between.visit.time.params <- c(10,3)
  
  # set both the true and initial model parameters
  ret.set.model.parameters <- set.model.parameters.func()
  I <- ret.set.model.parameters$I
  J <- ret.set.model.parameters$J
  L <- ret.set.model.parameters$L
  pi <- ret.set.model.parameters$pi
  lambda <- ret.set.model.parameters$lambda
  R <- ret.set.model.parameters$R
  Q <- ret.set.model.parameters$Q
  mu  <- ret.set.model.parameters$mu
  eta  <- ret.set.model.parameters$eta
  pi.tilde <- ret.set.model.parameters$pi.tilde
  Q.tilde <- ret.set.model.parameters$Q.tilde
  mu.tilde  <- ret.set.model.parameters$mu.tilde
  eta.tilde  <- ret.set.model.parameters$eta.tilde
  
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
  
  # visualize the data
  # visualization.func(H, a, tau.true, tau.obs, z.true, z.obs, y.obs, u.obs)
  
  # run the EM algorithm
  counter.em.max <- 20
  ret.EM.algorithm <- EM.algorithm.func(N, I, J, L, H, a, tau.obs, y.obs, u.obs, pi.tilde, Q.tilde, mu.tilde, eta.tilde, counter.em.max, debugging.mode)
  pi.tilde <- ret.EM.algorithm$pi.tilde
  Q.tilde <- ret.EM.algorithm$Q.tilde
  mu.tilde <- ret.EM.algorithm$mu.tilde
  eta.tilde <- ret.EM.algorithm$eta.tilde
  
  # print the results
  print(pi.tilde)
  print(Q.tilde)
  print(mu.tilde)
  print(eta.tilde)
  
}

