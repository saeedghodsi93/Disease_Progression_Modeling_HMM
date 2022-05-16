
# M-step
dthmm.M.step.func <- function(I, J, L, sufficient.pi, sufficient.Q, sufficient.mu, sufficient.eta) {
  
  # update pi
  pi.tilde <- array(numeric(), dim=c(I))
  for (i in 1:I) {
    pi.tilde[i] <- sufficient.pi[i] / sum(sufficient.pi)
  }
  
  # update Q
  Q.tilde <- array(numeric(), dim=c(I,I,L))
  for (i in 1:I) {
    for (k in 1:I) {
      for (l in 0:(L-1)) {
        Q.tilde[i,k,l+1] <- sufficient.Q[i,k,l+1] / sum(sufficient.Q[i,,l+1])
      }
    }
  }
  
  # update mu
  mu.tilde <- array(numeric(), dim=c(I))
  for (i in 1:I) {
    mu.tilde[i] <- sum(c(0:(J-1)) * sufficient.mu[i,]) / ((J-1) * sum(sufficient.mu[i,]))
  }
  
  # update eta
  eta.tilde <- array(numeric(), dim=c(J))
  for (j in 0:(J-1)) {
    eta.tilde[j+1] <- sum(c(0:(L-1)) * sufficient.eta[j+1,]) / ((L-1) * sum(sufficient.eta[j+1,]))
  }
  
  # for debugging
  # print(pi.tilde)
  # print(Q.tilde)
  # print(mu.tilde)
  # print(eta.tilde)
  
  ret <- list("pi.tilde" = pi.tilde, "Q.tilde" = Q.tilde, "mu.tilde" = mu.tilde, "eta.tilde" = eta.tilde)
  
  return(ret)  
}

