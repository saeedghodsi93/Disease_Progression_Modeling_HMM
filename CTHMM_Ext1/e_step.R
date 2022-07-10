
# E-step using the forward-backward algorithm
cthmm.ext1.E.step.func <- function(N, I, J, L, H, a, O, tau.obs, z.acc, y.obs, u.obs, pi.tilde, Q.tilde, mu.tilde, eta.tilde, eta.prime.tilde) {
  
  sufficient.pi <- array(0, c(I))
  sufficient.Q <- array(0, c(I,I,L,2))
  sufficient.mu <- array(0, c(I,J))
  sufficient.eta <- array(0, c(J,L))
  sufficient.eta.prime <- array(0, c(J,L))
  for (n in 1:N) {
    H.n <- H[n]
    
    # calculate the posterior probabilities of the end states
    ret.forward.backward <- cthmm.ext1.forward.backward.func(I, J, L, tau.obs, z.acc, y.obs, u.obs, O, pi.tilde, Q.tilde, mu.tilde, eta.tilde, eta.prime.tilde, H.n, n)
    gamma <- ret.forward.backward$gamma
    nu <- ret.forward.backward$nu
      
    # calculate the expected number of transitions and sojourn times in each period
    expectation.transition.mat <- array(0, c(H.n-1,I,I))
    expectation.sojourn.time <- array(0, c(H.n-1,I))
    for (t in 1:(H.n-1)) {
      
      # calculate the end state conditioned expected number of transitions and sojourn times in each period
      ret.end.conditioned.expectations <- cthmm.ext1.end.conditioned.expectations.func(I, Q.tilde, u.obs[n,t]+1, tau.obs[n,t+1]-tau.obs[n,t])
      transition.mat <- ret.end.conditioned.expectations$transition.mat
      sojourn.time <- ret.end.conditioned.expectations$sojourn.time
      
      # calculate the expected number of transitions
      for (i in 1:I) {
        for (k in 1:I) {
          for (i.tilde in 1:I) {
            for (k.tilde in 1:I) {
              expectation.transition.mat[t,i,k] <- expectation.transition.mat[t,i,k] + transition.mat[i,k,i.tilde,k.tilde] * nu[t,i.tilde,k.tilde]
            }
          }
        }
      }
      
      # calculate the expected sojourn times
      for (i in 1:I) {
        for (i.tilde in 1:I) {
          for (k.tilde in 1:I) {
            expectation.sojourn.time[t,i] <- expectation.sojourn.time[t,i] + sojourn.time[i,i.tilde,k.tilde] * nu[t,i.tilde,k.tilde]
          }
        }
      }
      
    }
    
    # update the sufficient statistics
    for (i in 1:I) {
      sufficient.pi[i] <- sufficient.pi[i] + gamma[1,i]
    }
    for (i in 1:I) {
      for (k in 1:I) {
        for (l in 0:(L-1)) {
          for (t in 1:(H.n-1)) {
            if (u.obs[n,t]==l) {
              sufficient.Q[i,k,l+1,1] <- sufficient.Q[i,k,l+1,1] + expectation.transition.mat[t,i,k]
              sufficient.Q[i,k,l+1,2] <- sufficient.Q[i,k,l+1,2] + expectation.sojourn.time[t,i]
            }
          }
        }
      }
    }
    for (i in 1:I) {
      for (j in 0:(J-1)) {
        for (t in 1:H.n) {
          sufficient.mu[i,j+1] <- sufficient.mu[i,j+1] + gamma[t,i] * sum(y.obs[n,t]==j)
        }
      }
    }
    for (j in 0:(J-1)) {
      for (l in 0:(L-1)) {
        for (t in 1:H.n) {
          sufficient.eta[j+1,l+1] <- sufficient.eta[j+1,l+1] + sum(y.obs[n,t]==j) * sum(u.obs[n,t]==l)
        }
      }
    }
    
  }
  
  ret <- list("sufficient.pi" = sufficient.pi, "sufficient.Q" = sufficient.Q, "sufficient.mu" = sufficient.mu, "sufficient.eta" = sufficient.eta)
  
  return(ret)
  
}

