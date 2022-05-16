
# emission model
dthmm.emission.prob.func <- function(J, mu.in, j.in) {
  
  ret <- choose(J-1, j.in) * (mu.in ^ j.in) * ((1-mu.in) ^ (J-1 - j.in))
  
  return(ret)
  
}


# intervention model
dthmm.intervention.prob.func <- function(L, eta.in, l.in) {
  
  ret <- choose(L-1, l.in) * (eta.in ^ l.in) * ((1-eta.in) ^ (L-1 - l.in))
  
  return(ret)
  
}


# the forward-backward algorithm for calculating the end-state posterior probabilities
dthmm.forward.backward.func <- function(I, J, L, y, u, pi.tilde, Q.tilde, mu.tilde, eta.tilde, H.n, n) {
  
  # forward probabilities
  alpha <- array(numeric(), c(H.n,I))
  for (i in 1:I) {
    alpha[1,i] <- pi.tilde[i] * dthmm.emission.prob.func(J, mu.tilde[i], y[n,1]) * dthmm.intervention.prob.func(L, eta.tilde[y[n,1]+1], u[n,1])
  }
  for (t in 1:(H.n-1)) {
    for (i in 1:I) {
      alpha[t+1,i] <- 0
      for (k in 1:I) {
        alpha[t+1,i] <- alpha[t+1,i] + Q.tilde[k,i, u[n,t]+1] * dthmm.emission.prob.func(J, mu.tilde[i], y[n,t+1]) * dthmm.intervention.prob.func(L, eta.tilde[y[n,t+1]+1], u[n,t+1]) * alpha[t,k]
      }
    }
  }
  
  # backward probabilities
  beta <- array(numeric(), c(H.n-1,I))
  for (i in 1:I) {
    beta[H.n-1,i] <- 0
    for (k in 1:I) {
      beta[H.n-1,i] <- beta[H.n-1,i] + Q.tilde[i,k, u[n,H.n-1]+1] * dthmm.emission.prob.func(J, mu.tilde[k], y[n,H.n]) * dthmm.intervention.prob.func(L, eta.tilde[y[n,H.n]+1], u[n,H.n])
    }
  }
  for (t in (H.n-1):2) {
    for (i in 1:I) {
      beta[t-1,i] <- 0
      for (k in 1:I) {
        beta[t-1,i] <- beta[t-1,i] + Q.tilde[i,k, u[n,t-1]+1] * dthmm.emission.prob.func(J, mu.tilde[k], y[n,t]) * dthmm.intervention.prob.func(L, eta.tilde[y[n,t]+1], u[n,t]) * beta[t,k]
      } 
    }
  }
  
  # posterior probabilities
  gamma <- array(numeric(), c(H.n,I))
  for (i in 1:I) {
    gamma[H.n,i] <- alpha[H.n,i] / sum(alpha[H.n,])
  }
  for (t in 1:(H.n-1)) {
    temp <- array(numeric(), c(I))
    for (i in 1:I) {
      temp[i] <- alpha[t,i] * beta[t,i]
    }
    for (i in 1:I) {
      gamma[t,i] <- temp[i] / sum(temp)
    }
  }
  nu <- array(numeric(), c(H.n-1,I,I))
  temp <- array(numeric(), c(I,I))
  for (i in 1:I) {
    for (k in 1:I) {
      temp[i,k] <- Q.tilde[i,k, u[n,H.n-1]+1] * dthmm.emission.prob.func(J, mu.tilde[k], y[n,H.n]) * dthmm.intervention.prob.func(L, eta.tilde[y[n,H.n]+1], u[n,H.n]) * alpha[H.n-1,i]
    }
  }
  for (i in 1:I) {
    for (k in 1:I) {
      nu[H.n-1,i,k] <- temp[i,k] / sum(temp)
    }
  }
  for (t in 1:(H.n-2)) {
    temp <- array(numeric(), c(I,I))
    for (i in 1:I) {
      for (k in 1:I) {
        temp[i,k] <- Q.tilde[i,k, u[n,t]+1] * dthmm.emission.prob.func(J, mu.tilde[k], y[n,t+1]) * dthmm.intervention.prob.func(L, eta.tilde[y[n,t+1]+1], u[n,t+1]) * alpha[t,i] * beta[t+1,k]
      }
    }
    for (i in 1:I) {
      for (k in 1:I) {
        nu[t,i,k] <- temp[i,k] / sum(temp)
      }
    }
  }
  
  ret <- list("alpha" = alpha, "beta" = beta, "gamma" = gamma, "nu" = nu)
  
  return(ret)
  
}

