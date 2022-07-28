
# transition model
cthmm.ext2.transition.prob.func <- function(Q.in, i.in, k.in, l.in, tau.in) {
  
  transition.mat <- expm(tau.in * Q.in[, , l.in])
  
  ret <- transition.mat[i.in, k.in]
  
  return(ret)
  
}


# emission model
cthmm.ext2.emission.prob.func <- function(J, mu.in, j.in) {
  
  ret <- choose(J-1, j.in) * (mu.in ^ j.in) * ((1-mu.in) ^ (J-1 - j.in))
  
  return(ret)
  
}


# intervention model
cthmm.ext2.intervention.prob.func <- function(L, eta.in, l.in) {
  
  ret <- choose(L-1, l.in) * (eta.in ^ l.in) * ((1-eta.in) ^ (L-1 - l.in))
  
  return(ret)
  
}


# the forward-backward algorithm for calculating the end-state posterior probabilities
cthmm.ext2.forward.backward.func <- function(I, J, L, tau.obs, z.acc, y.obs, u.obs, O, pi.in, Q.in, mu.in, eta.in, eta.prime.in, H.in, n) {
  
  # forward probabilities
  alpha <- array(numeric(), c(H.in,I))
  if (O[n,1]==0) {
    for (i in 1:I) {
      alpha[1,i] <- pi.in[i] * cthmm.ext2.emission.prob.func(J, mu.in[i,1], y.obs[n,1]) * cthmm.ext2.intervention.prob.func(L, eta.in[y.obs[n,1]+1], u.obs[n,1])
    }
  } else{
    alpha[1, z.acc[n,1]] <- pi.in[z.acc[n,1]] * cthmm.ext2.intervention.prob.func(L, eta.in[z.acc[n,1]], u.obs[n,1])
  }
  for (t in 1:(H.in-1)) {
    if ((O[n,t]==0) & (O[n,t+1]==0)) {
      for (i in 1:I) {
        alpha[t+1,i] <- 0
        for (k in 1:I) {
          alpha[t+1,i] <- alpha[t+1,i] + cthmm.ext2.transition.prob.func(Q.in, k, i, u.obs[n,t]+1, tau.obs[n,t+1]-tau.obs[n,t]) * cthmm.ext2.emission.prob.func(J, mu.in[i, u.obs[n,t]+1], y.obs[n,t+1]) * cthmm.ext2.intervention.prob.func(L, eta.in[y.obs[n,t+1]+1], u.obs[n,t+1]) * alpha[t,k]
        }
      }
    } else if ((O[n,t]==1) & (O[n,t+1]==0)) {
      for (i in 1:I) {
        alpha[t+1,i] <- cthmm.ext2.transition.prob.func(Q.in, z.acc[n,t], i, u.obs[n,t]+1, tau.obs[n,t+1]-tau.obs[n,t]) * cthmm.ext2.emission.prob.func(J, mu.in[i, u.obs[n,t]+1], y.obs[n,t+1]) * cthmm.ext2.intervention.prob.func(L, eta.in[y.obs[n,t+1]+1], u.obs[n,t+1]) * alpha[t, z.acc[n,t]]
      }
    } else if ((O[n,t]==0) & (O[n,t+1]==1)) {
      alpha[t+1, z.acc[n,t+1]] <- 0
      for (k in 1:I) {
        alpha[t+1, z.acc[n,t+1]] <- alpha[t+1, z.acc[n,t+1]] + cthmm.ext2.transition.prob.func(Q.in, k, z.acc[n,t+1], u.obs[n,t]+1, tau.obs[n,t+1]-tau.obs[n,t]) * cthmm.ext2.intervention.prob.func(L, eta.prime.in[z.acc[n,t+1]], u.obs[n,t+1]) * alpha[t,k]
      }
    } else {
      alpha[t+1, z.acc[n,t+1]] <- cthmm.ext2.transition.prob.func(Q.in, z.acc[n,t], z.acc[n,t+1], u.obs[n,t]+1, tau.obs[n,t+1]-tau.obs[n,t]) * cthmm.ext2.intervention.prob.func(L, eta.prime.in[z.acc[n,t+1]], u.obs[n,t+1]) * alpha[t, z.acc[n,t]]
    }
  }
  
  # backward probabilities
  beta <- array(numeric(), c(H.in-1,I))
  if ((O[n, H.in-1]==0) & (O[n,H.in]==0)) {
    for (i in 1:I) {
      beta[H.in-1,i] <- 0
      for (k in 1:I) {
        beta[H.in-1,i] <- beta[H.in-1,i] + cthmm.ext2.transition.prob.func(Q.in, i, k, u.obs[n,H.in-1]+1, tau.obs[n,H.in]-tau.obs[n,H.in-1]) * cthmm.ext2.emission.prob.func(J, mu.in[k, u.obs[n,H.in-1]+1], y.obs[n,H.in]) * cthmm.ext2.intervention.prob.func(L, eta.in[y.obs[n,H.in]+1], u.obs[n,H.in])
      }
    }
  } else if ((O[n, H.in-1]==0) & (O[n,H.in]==1)) {
    for (i in 1:I) {
      beta[H.in-1,i] <- cthmm.ext2.transition.prob.func(Q.in, i, z.acc[n,H.in], u.obs[n,H.in-1]+1, tau.obs[n,H.in]-tau.obs[n,H.in-1]) * cthmm.ext2.intervention.prob.func(L, eta.prime.in[z.acc[n,H.in]], u.obs[n,H.in])
    }
  } else if ((O[n, H.in-1]==1) & (O[n,H.in]==0)) {
    beta[H.in-1, z.acc[n,H.in-1]] <- 0
    for (k in 1:I) {
      beta[H.in-1, z.acc[n,H.in-1]] <- beta[H.in-1, z.acc[n,H.in-1]] + cthmm.ext2.transition.prob.func(Q.in, z.acc[n,H.in-1], k, u.obs[n,H.in-1]+1, tau.obs[n,H.in]-tau.obs[n,H.in-1]) * cthmm.ext2.emission.prob.func(J, mu.in[k, u.obs[n,H.in-1]+1], y.obs[n,H.in]) * cthmm.ext2.intervention.prob.func(L, eta.in[y.obs[n,H.in]+1], u.obs[n,H.in])
    }
  } else {
    beta[H.in-1, z.acc[n,H.in-1]] <- cthmm.ext2.transition.prob.func(Q.in, z.acc[n,H.in-1], z.acc[n,H.in], u.obs[n,H.in-1]+1, tau.obs[n,H.in]-tau.obs[n,H.in-1]) * cthmm.ext2.intervention.prob.func(L, eta.prime.in[z.acc[n,H.in]], u.obs[n,H.in])
  }
  for (t in (H.in-1):2) {
    if ((O[n,t-1]==0) & (O[n,t]==0)) {
      for (i in 1:I) {
        beta[t-1,i] <- 0
        for (k in 1:I) {
          beta[t-1,i] <- beta[t-1,i] + cthmm.ext2.transition.prob.func(Q.in, i, k, u.obs[n,t-1]+1, tau.obs[n,t]-tau.obs[n,t-1]) * cthmm.ext2.emission.prob.func(J, mu.in[k, u.obs[n,t-1]+1], y.obs[n,t]) * cthmm.ext2.intervention.prob.func(L, eta.in[y.obs[n,t]+1], u.obs[n,t]) * beta[t,k]
        } 
      }
    } else if ((O[n,t-1]==0) & (O[n,t]==1)) {
      for (i in 1:I) {
        beta[t-1,i] <- cthmm.ext2.transition.prob.func(Q.in, i, z.acc[n,t], u.obs[n,t-1]+1, tau.obs[n,t]-tau.obs[n,t-1]) * cthmm.ext2.intervention.prob.func(L, eta.prime.in[z.acc[n,t]], u.obs[n,t]) * beta[t, z.acc[n,t]]
      }
    } else if ((O[n,t-1]==1) & (O[n,t]==0)) {
      beta[t-1, z.acc[n,t-1]] <- 0
      for (k in 1:I) {
        beta[t-1, z.acc[n,t-1]] <- beta[t-1, z.acc[n,t-1]] + cthmm.ext2.transition.prob.func(Q.in, z.acc[n,t-1], k, u.obs[n,t-1]+1, tau.obs[n,t]-tau.obs[n,t-1]) * cthmm.ext2.emission.prob.func(J, mu.in[k, u.obs[n,t-1]+1], y.obs[n,t]) * cthmm.ext2.intervention.prob.func(L, eta.in[y.obs[n,t]+1], u.obs[n,t]) * beta[t,k]
      }
    } else {
      beta[t-1, z.acc[n,t-1]] <- cthmm.ext2.transition.prob.func(Q.in, z.acc[n,t-1], z.acc[n,t], u.obs[n,t-1]+1, tau.obs[n,t]-tau.obs[n,t-1]) * cthmm.ext2.intervention.prob.func(L, eta.prime.in[z.acc[n,t]], u.obs[n,t]) * beta[t, z.acc[n,t]]
    }
  }
  
  # posterior probability of z_t
  gamma <- array(numeric(), c(H.in,I))
  for (i in 1:I) {
    if (O[n,H.in]==0) {
      gamma[H.in,i] <- alpha[H.in,i] / sum(alpha[H.in,])
    } else {
      gamma[H.in,i] <- sum(i == z.acc[n,H.in])
    }
  }
  for (t in 1:(H.in-1)) {
    temp <- array(numeric(), c(I))
    for (i in 1:I) {
      if (O[n,t]==0) {
        temp[i] <- alpha[t,i] * beta[t,i]
      } else {
        temp[i] <- sum(i == z.acc[n,t])
      }
    }
    for (i in 1:I) {
      gamma[t,i] <- temp[i] / sum(temp)
    }
  }
  
  # posterior probability of (z_t, z_{t+1})
  nu <- array(numeric(), c(H.in-1,I,I))
  if ((O[n,H.in-1]==0) & (O[n,H.in]==0)) {
    temp <- array(numeric(), c(I,I))
    for (i in 1:I) {
      for (k in 1:I) {
        temp[i,k] <- cthmm.ext2.transition.prob.func(Q.in, i, k, u.obs[n,H.in-1]+1, tau.obs[n,H.in]-tau.obs[n,H.in-1]) * cthmm.ext2.emission.prob.func(J, mu.in[k, u.obs[n,H.in-1]+1], y.obs[n,H.in]) * cthmm.ext2.intervention.prob.func(L, eta.in[y.obs[n,H.in]+1], u.obs[n,H.in]) * alpha[H.in-1,i]
      }
    }
    for (i in 1:I) {
      for (k in 1:I) {
        nu[H.in-1,i,k] <- temp[i,k] / sum(temp)
      }
    }
  } else if ((O[n,H.in-1]==0) & (O[n,H.in]==1)) {
    for (i in 1:I) {
      for (k in 1:I) {
        nu[H.in-1,i,k] <- gamma[H.in-1, i] * sum(k == z.acc[n,H.in])
      }
    }
  } else if ((O[n,H.in-1]==1) & (O[n,H.in]==0)) {
    for (i in 1:I) {
      for (k in 1:I) {
        nu[H.in-1,i,k] <- gamma[H.in, k] * sum(i == z.acc[n,H.in-1])
      }
    }
  } else {
    for (i in 1:I) {
      for (k in 1:I) {
        nu[H.in-1,i,k] <- sum(i == z.acc[n,H.in-1]) * sum(k == z.acc[n,H.in])
      }
    }
  }
  for (t in 1:(H.in-2)) {
    if ((O[n,t]==0) & (O[n,t+1]==0)) {
      temp <- array(numeric(), c(I,I))
      for (i in 1:I) {
        for (k in 1:I) {
          temp[i,k] <- cthmm.ext2.transition.prob.func(Q.in, i, k, u.obs[n,t]+1, tau.obs[n,t+1]-tau.obs[n,t]) * cthmm.ext2.emission.prob.func(J, mu.in[k, u.obs[n,t]+1], y.obs[n,t+1]) * cthmm.ext2.intervention.prob.func(L, eta.in[y.obs[n,t+1]+1], u.obs[n,t+1]) * alpha[t,i] * beta[t+1,k]
        }
      }
      for (i in 1:I) {
        for (k in 1:I) {
          nu[t,i,k] <- temp[i,k] / sum(temp)
        }
      }
    } else if ((O[n,t]==0) & (O[n,t+1]==1)) {
      for (i in 1:I) {
        for (k in 1:I) {
          nu[t,i,k] <- gamma[t,i] * sum(k == z.acc[n,t+1])
        }
      }
    } else if ((O[n,t]==1) & (O[n,t+1]==0)) {
      for (i in 1:I) {
        for (k in 1:I) {
          nu[t,i,k] <- gamma[t+1,k] * sum(i == z.acc[n,t])
        }
      }
    } else {
      for (i in 1:I) {
        for (k in 1:I) {
          nu[t,i,k] <- sum(i == z.acc[n,t]) * sum(k == z.acc[n,t+1])
        }
      }
    }
  }

  # for debugging  
  # if (n == 6) {
  #   print(alpha)
  #   print(beta)
  #   print(gamma)
  #   print(aperm(nu, c(2,3,1)))
  #   print(z.obs[n,])
  #   print(z.acc[n,])
  #   print(y.obs[n,])
  #   print(u.obs[n,])
  # }
  
  ret <- list("alpha" = alpha, "beta" = beta, "gamma" = gamma, "nu" = nu)
  
  return(ret)
  
}

