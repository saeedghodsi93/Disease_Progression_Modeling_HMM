
# E-step using the forward-backward algorithm or the Gibbs sampling approach
dthmm.E.step.func <- function(N, M, B, I, J, L, H, y, u, pi.tilde, Q.tilde, mu.tilde, eta.tilde) {
  
  sufficient.pi <- array(0, c(I))
  sufficient.Q <- array(0, c(I,I,L))
  sufficient.mu <- array(0, c(I,J))
  sufficient.eta <- array(0, c(J,L))
  for (n in 1:N) {
    H.n <- H[n]
    
    # calculate the posterior probabilities of the end states (use direct forward-backward if M=0)
    if (M == 0) {
      ret.forward.backward <- dthmm.forward.backward.func(I, J, L, y, u, pi.tilde, Q.tilde, mu.tilde, eta.tilde, H.n, n)
      gamma <- ret.forward.backward$gamma
      nu <- ret.forward.backward$nu
    } else {
      ret.gibbs.sampling <- dthmm.gibbs.sampling.func(I, J, L, y, u, pi.tilde, Q.tilde, mu.tilde, eta.tilde, H.n, n, M, B)
      gamma <- ret.gibbs.sampling$gamma
      nu <- ret.gibbs.sampling$nu
    }
    
    # update the sufficient statistics
    for (i in 1:I) {
      sufficient.pi[i] <- sufficient.pi[i] + gamma[1,i]
    }
    for (i in 1:I) {
      for (k in 1:I) {
        for (l in 0:(L-1)) {
          for (t in 1:(H.n-1)) {
            if (u[n,t]==l) {
              sufficient.Q[i,k,l+1] <- sufficient.Q[i,k,l+1] + nu[t,i,k]
            }
          }
        }
      }
    }
    for (i in 1:I) {
      for (j in 0:(J-1)) {
        for (t in 1:H.n) {
          sufficient.mu[i,j+1] <- sufficient.mu[i,j+1] + gamma[t,i] * sum(y[n,t]==j)
        }
      }
    }
    for (j in 0:(J-1)) {
      for (l in 0:(L-1)) {
        for (t in 1:H.n) {
          sufficient.eta[j+1,l+1] <- sufficient.eta[j+1,l+1] + sum(y[n,t]==j) * sum(u[n,t]==l)
        }
      }
    }
    
  }
  
  ret <- list("sufficient.pi" = sufficient.pi, "sufficient.Q" = sufficient.Q, "sufficient.mu" = sufficient.mu, "sufficient.eta" = sufficient.eta)
  
  return(ret)
  
}

