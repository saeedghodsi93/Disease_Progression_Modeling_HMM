
# transition model
cthmm.transition.prob.func <- function(Q.in, i.in, k.in, l.in, tau.in) {
  
  transition.mat <- expm(tau.in * Q.in[, , l.in])
  
  ret <- transition.mat[i.in, k.in]
  
  return(ret)
  
}


# emission model
cthmm.emission.prob.func <- function(J, mu.in, j.in) {
  
  ret <- choose(J-1, j.in) * (mu.in ^ j.in) * ((1-mu.in) ^ (J-1 - j.in))
  
  return(ret)
  
}


# intervention model
cthmm.intervention.prob.func <- function(L, eta.in, l.in) {
  
  ret <- choose(L-1, l.in) * (eta.in ^ l.in) * ((1-eta.in) ^ (L-1 - l.in))
  
  return(ret)
  
}


# the forward-backward algorithm for calculating the end-state posterior probabilities
cthmm.forward.backward.func <- function(I, J, L, tau.obs, y.obs, u.obs, pi.in, Q.in, mu.in, eta.in, H.in, n) {
  
  # forward probabilities
  alpha <- array(numeric(), c(H.in,I))
  for (i in 1:I) {
    alpha[1,i] <- pi.in[i] * cthmm.emission.prob.func(J, mu.in[i], y.obs[n,1]) * cthmm.intervention.prob.func(L, eta.in[y.obs[n,1]+1], u.obs[n,1])
  }
  for (t in 1:(H.in-1)) {
    for (i in 1:I) {
      alpha[t+1,i] <- 0
      for (k in 1:I) {
        alpha[t+1,i] <- alpha[t+1,i] + cthmm.transition.prob.func(Q.in, k, i, u.obs[n,t]+1, tau.obs[n,t+1]-tau.obs[n,t]) * cthmm.emission.prob.func(J, mu.in[i], y.obs[n,t+1]) * cthmm.intervention.prob.func(L, eta.in[y.obs[n,t+1]+1], u.obs[n,t+1]) * alpha[t,k]
      }
    }
  }
  
  # backward probabilities
  beta <- array(numeric(), c(H.in-1,I))
  for (i in 1:I) {
    beta[H.in-1,i] <- 0
    for (k in 1:I) {
      beta[H.in-1,i] <- beta[H.in-1,i] + cthmm.transition.prob.func(Q.in, i, k, u.obs[n,H.in-1]+1, tau.obs[n,H.in]-tau.obs[n,H.in-1]) * cthmm.emission.prob.func(J, mu.in[k], y.obs[n,H.in]) * cthmm.intervention.prob.func(L, eta.in[y.obs[n,H.in]+1], u.obs[n,H.in])
    }
  }
  for (t in (H.in-1):2) {
    for (i in 1:I) {
      beta[t-1,i] <- 0
      for (k in 1:I) {
        beta[t-1,i] <- beta[t-1,i] + cthmm.transition.prob.func(Q.in, i, k, u.obs[n,t-1]+1, tau.obs[n,t]-tau.obs[n,t-1]) * cthmm.emission.prob.func(J, mu.in[k], y.obs[n,t]) * cthmm.intervention.prob.func(L, eta.in[y.obs[n,t]+1], u.obs[n,t]) * beta[t,k]
      } 
    }
  }
  
  # posterior probability of z_t
  gamma <- array(numeric(), c(H.in,I))
  for (i in 1:I) {
    gamma[H.in,i] <- alpha[H.in,i] / sum(alpha[H.in,])
  }
  for (t in 1:(H.in-1)) {
    temp <- array(numeric(), c(I))
    for (i in 1:I) {
      temp[i] <- alpha[t,i] * beta[t,i]
    }
    for (i in 1:I) {
      gamma[t,i] <- temp[i] / sum(temp)
    }
  }
  
  # posterior probability of (z_t, z_{t+1})
  nu <- array(numeric(), c(H.in-1,I,I))
  temp <- array(numeric(), c(I,I))
  for (i in 1:I) {
    for (k in 1:I) {
      temp[i,k] <- cthmm.transition.prob.func(Q.in, i, k, u.obs[n,H.in-1]+1, tau.obs[n,H.in]-tau.obs[n,H.in-1]) * cthmm.emission.prob.func(J, mu.in[k], y.obs[n,H.in]) * cthmm.intervention.prob.func(L, eta.in[y.obs[n,H.in]+1], u.obs[n,H.in]) * alpha[H.in-1,i]
    }
  }
  for (i in 1:I) {
    for (k in 1:I) {
      nu[H.in-1,i,k] <- temp[i,k] / sum(temp)
    }
  }
  for (t in 1:(H.in-2)) {
    temp <- array(numeric(), c(I,I))
    for (i in 1:I) {
      for (k in 1:I) {
        temp[i,k] <- cthmm.transition.prob.func(Q.in, i, k, u.obs[n,t]+1, tau.obs[n,t+1]-tau.obs[n,t]) * cthmm.emission.prob.func(J, mu.in[k], y.obs[n,t+1]) * cthmm.intervention.prob.func(L, eta.in[y.obs[n,t+1]+1], u.obs[n,t+1]) * alpha[t,i] * beta[t+1,k]
      }
    }
    for (i in 1:I) {
      for (k in 1:I) {
        nu[t,i,k] <- temp[i,k] / sum(temp)
      }
    }
  }

  # for debugging  
  # if (n == 100) {
  #   print(alpha)
  #   print(beta)
  #   print(gamma)
  #   print(aperm(nu, c(2,3,1)))
  #   print(z.obs[n,])
  #   print(y.obs[n,])
  # }
  
  ret <- list("alpha" = alpha, "beta" = beta, "gamma" = gamma, "nu" = nu)
  
  return(ret)
  
}

