
# M-step
cthmm.ext2.M.step.func <- function(I, J, L, sufficient.pi, sufficient.Q, sufficient.mu, sufficient.eta, sufficient.eta.prime) {
  
  # update pi
  pi.tilde <- array(numeric(), dim=c(I))
  for (i in 1:I) {
    pi.tilde[i] <- sufficient.pi[i] / sum(sufficient.pi)
  }
  
  # update Q
  Q.tilde <- array(numeric(), dim=c(I,I,L))
  for (i in 1:I) {
    for (l in 0:(L-1)) {
      temp <- 0
      for (k in 1:I) {
        if (k != i) {
          Q.tilde[i,k,l+1] <- sufficient.Q[i,k,l+1,1] / sufficient.Q[i,k,l+1,2]
          temp <- temp + Q.tilde[i,k,l+1]
        }
      }
      Q.tilde[i,i,l+1] <- -temp 
    }
  }
  
  # update mu
  mu.tilde <- array(numeric(), dim=c(I,L))
  for (i in 1:I) {
    for (l in 0:(L-1)) {
      mu.tilde[i,l+1] <- sum(c(0:(J-1)) * sufficient.mu[i,l+1,]) / ((J-1) * sum(sufficient.mu[i,l+1,]))
    }
  }
  
  # update eta
  eta.tilde <- array(numeric(), dim=c(J))
  for (j in 0:(J-1)) {
    eta.tilde[j+1] <- sum(c(0:(L-1)) * sufficient.eta[j+1,]) / ((L-1) * sum(sufficient.eta[j+1,]))
  }
  
  # update eta.prime
  eta.prime.tilde <- array(numeric(), dim=c(I))
  for (i in 1:I) {
    eta.prime.tilde[i] <- sum(c(0:(L-1)) * sufficient.eta.prime[i,]) / ((L-1) * sum(sufficient.eta.prime[i,]))
  }
  
  # for debugging
  # print(pi.tilde)
  # print(Q.tilde)
  # print(mu.tilde)
  # print(eta.tilde)
  # print(eta.prime.tilde)
  
  ret <- list("pi.tilde" = pi.tilde, "Q.tilde" = Q.tilde, "mu.tilde" = mu.tilde, "eta.tilde" = eta.tilde, "eta.prime.tilde" = eta.prime.tilde)
  
  return(ret)  
}

