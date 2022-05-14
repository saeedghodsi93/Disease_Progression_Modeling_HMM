
# EM algorithm
EM.algorithm.func <- function(N, I, J, L, H, a, tau.obs, y.obs, u.obs, pi.init, Q.init, mu.init, eta.init, counter.em.max, debugging.mode) {
  
  # arrays for keeping track of the parameters
  pi.hat <- array(numeric(), dim=c(counter.em.max+1,I))
  Q.hat <- array(numeric(), dim=c(counter.em.max+1,I,I,L))
  mu.hat <- array(numeric(), dim=c(counter.em.max+1,I))
  eta.hat <- array(numeric(), dim=c(counter.em.max+1,J))
  
  # save the initial parameters
  pi.hat[1,] <- pi.init
  Q.hat[1,,,] <- Q.init
  mu.hat[1,] <- mu.init
  eta.hat[1,] <- eta.init
  
  # EM iterations
  pi.tilde <- pi.init
  Q.tilde <- Q.init
  mu.tilde <- mu.init
  eta.tilde <- eta.init
  counter.em <- 1
  while (TRUE) {
    
    # add a progress bar
    progress.bar(counter.em, counter.em.max, TRUE, debugging.mode)
    
    # E-step 
    ret.E.step <- E.step.func(N, I, J, L, H, a, tau.obs, y.obs, u.obs, pi.tilde, Q.tilde, mu.tilde, eta.tilde)
    sufficient.pi <- ret.E.step$sufficient.pi
    sufficient.Q <- ret.E.step$sufficient.Q
    sufficient.mu <- ret.E.step$sufficient.mu
    sufficient.eta <- ret.E.step$sufficient.eta
    
    # M-step
    ret.M.step <- M.step.func(I, J, L, sufficient.pi, sufficient.Q, sufficient.mu, sufficient.eta)
    pi.tilde <- ret.M.step$pi.tilde
    Q.tilde <- ret.M.step$Q.tilde
    mu.tilde <- ret.M.step$mu.tilde
    eta.tilde <- ret.M.step$eta.tilde
    
    # save the current parameters
    pi.hat[counter.em+1,] <- pi.tilde
    Q.hat[counter.em+1,,,] <- Q.tilde
    mu.hat[counter.em+1,] <- mu.tilde
    eta.hat[counter.em+1,] <- eta.tilde
    
    # check the iterator and run the progress bar function
    counter.em <- counter.em + 1
    progress.bar(counter.em, counter.em.max, FALSE, debugging.mode)
    if (counter.em > counter.em.max) {
      break
    }
    
  }
  
  ret <- list("pi.hat" = pi.hat, "Q.hat" = Q.hat, "mu.hat" = mu.hat, "eta.hat" = eta.hat)
  
  return(ret)
  
}

