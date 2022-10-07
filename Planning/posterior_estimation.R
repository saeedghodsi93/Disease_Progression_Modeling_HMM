
# the forward algorithm for calculating the end-state posterior probabilities
planning.forward.func <- function(I, J, L, tau.obs, y.obs, u.obs, pi.in, Q.in, mu.in, eta.in, H.in, n, start.period) {
  
  # forward probabilities
  alpha <- array(numeric(), c(H.in,I))
  for (i in 1:I) {
    alpha[start.period,i] <- pi.in[i] * cthmm.emission.prob.func(J, mu.in[i], y.obs[n,start.period]) * cthmm.intervention.prob.func(L, eta.in[y.obs[n,start.period]+1], u.obs[n,start.period])
  }
  for (t in start.period:(H.in-1)) {
    for (i in 1:I) {
      alpha[t+1,i] <- 0
      for (k in 1:I) {
        alpha[t+1,i] <- alpha[t+1,i] + cthmm.transition.prob.func(Q.in, k, i, u.obs[n,t]+1, tau.obs[n,t+1]-tau.obs[n,t]) * cthmm.emission.prob.func(J, mu.in[i], y.obs[n,t+1]) * cthmm.intervention.prob.func(L, eta.in[y.obs[n,t+1]+1], u.obs[n,t+1]) * alpha[t,k]
      }
    }
  }
  
  ret <- list("alpha" = alpha)
  
  return(ret)
  
}

# estimate the posterior distribution of the health state in the last time period using only the last physician observation
planning.posterior.estimation.lastobs.func <- function(N.test, I, J, L, H.test, a.test, y.obs.test, pi.hat, mu.hat) {
  
  # estimate the joint probability
  prob.joint <- array(0, c(N.test,I))
  for (n in 1:N.test) {
    for (i in 1:I) {
      
      # calculate the joint probabilities
      H.n <- H.test[n]
      prob.joint[n,i] <- cthmm.emission.prob.func(J, mu.hat[i], y.obs.test[n,H.n]) * pi.hat[i]
      
    }
  }
  
  # estimate the posterior probability
  prob.posterior <- array(0, c(N.test,I))
  for (n in 1:N.test) {
    for (i in 1:I) {
      prob.posterior[n,i] <- prob.joint[n,i] / sum(prob.joint[n,])
    }
  }
  
  ret <- list("prob.joint" = prob.joint, "prob.posterior" = prob.posterior) 
  
  return(ret)  
  
}


# estimate the posterior distribution of the health state in the last time period using all the historical data
planning.posterior.estimation.history.func <- function(N.test, I, J, L, H.test, a.test, tau.obs.test, y.obs.test, u.obs.test, pi.hat, Q.hat, mu.hat, eta.hat) {
  
  # estimate the joint probability
  prob.joint <- array(0, c(N.test,I))
  for (n in 1:N.test) {
    
    # calculate the forward probabilities
    H.n <- H.test[n]
    ret.forward.backward <- planning.forward.func(I, J, L, tau.obs.test, y.obs.test, u.obs.test, pi.hat, Q.hat, mu.hat, eta.hat, H.n, n, 1)
    alpha <- ret.forward.backward$alpha
    
    for (i in 1:I) {
      
      # add the joint probabilities
      for (k in 1:I) {
        prob.joint[n,i] <- prob.joint[n,i] + cthmm.emission.prob.func(J, mu.hat[i], y.obs.test[n,H.n]) * cthmm.transition.prob.func(Q.hat, k, i, u.obs.test[n,H.n-1]+1, tau.obs.test[n,H.n]-tau.obs.test[n,H.n-1]) * alpha[H.n-1,k] 
      }
      
    }
  }
  
  # estimate the posterior probability
  prob.posterior <- array(0, c(N.test,I))
  for (n in 1:N.test) {
    for (i in 1:I) {
      prob.posterior[n,i] <- prob.joint[n,i] / sum(prob.joint[n,])
    }
  }
  
  
  ret <- list("prob.joint" = prob.joint, "prob.posterior" = prob.posterior) 
  
  return(ret)  
  
}

