
# estimate the posterior distribution of the health state in the last time period using all the historical data
planning.posterior.estimation.history.func <- function(N.test, I, J, L, H.test, a.test, tau.obs.test, y.obs.test, u.obs.test, pi.hat, Q.hat, mu.hat, eta.hat) {
  
  # estimate the joint probability
  prob.joint <- array(0, c(N.test,I))
  for (n in 1:N.test) {
    
    # calculate the forward probabilities
    H.n <- H.test[n]
    ret.forward.backward <- cthmm.forward.backward.func(I, J, L, tau.obs.test, y.obs.test, u.obs.test, pi.hat, Q.hat, mu.hat, eta.hat, H.n, n)
    alpha <- ret.forward.backward$alpha
    
    for (i in 1:I) {
      
      # add the joint probabilities
      for (k in 1:I) {
        prob.joint[n,i] <- prob.joint[n,i] + cthmm.transition.prob.func(Q.hat, k, i, u.obs.test[n,H.n-1]+1, tau.obs.test[n,H.n]-tau.obs.test[n,H.n-1]) * cthmm.emission.prob.func(J, mu.hat[i], y.obs.test[n,H.n]) * alpha[H.n-1,k] 
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
