
# EM algorithm
EM.algorithm.func <- function(N, I, J, L, H, a, tau.obs, y.obs, u.obs, pi.tilde, Q.tilde, mu.tilde, eta.tilde, counter.em.max, debugging.mode) {
  
  # EM iterations
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
    
    # check the iterator and run the progress bar function
    counter.em <- counter.em + 1
    progress.bar(counter.em, counter.em.max, FALSE, debugging.mode)
    if (counter.em > counter.em.max) {
      break
    }
    
  }
  
  ret <- list("pi.tilde" = pi.tilde, "Q.tilde" = Q.tilde, "mu.tilde" = mu.tilde, "eta.tilde" = eta.tilde)
  
  return(ret)  
  
}

