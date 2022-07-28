
# calculate the expected disutility for the test patients
planning.expected.disutility.func <- function(N.test, I, J, L, H.test, a.test, z.obs.test, prob.posterior.lastobs, prob.posterior.history, Q.in, delta.tau.planning, beta.planning) {
  
  # calculate the expected disutility, assuming z is known at the planning time
  expected.disutility.optimal <- array(0, c(N.test,L))
  best.choice.optimal <- array(0, c(N.test))
  actual.disutility.optimal <- array(0, c(N.test))
  for (n in 1:N.test) {
    H.n <- H.test[n]
    
    for (l in 0:(L-1)) {
      for (k in 1:I) {
        expected.disutility.optimal[n,l+1] <- expected.disutility.optimal[n,l+1] + beta.planning[k] * cthmm.transition.prob.func(Q.in, z.obs.test[n,H.n], k, l+1, delta.tau.planning)
      }
    }
    
    best.choice.optimal[n] <- which.min(expected.disutility.optimal[n,]) - 1
    
    actual.disutility.optimal[n] <- expected.disutility.optimal[n, best.choice.optimal[n]+1]
      
  }
  
  # calculate the expected disutility, using only the last observation for estimating the distribution of z at the planning time
  expected.disutility.lastobs <- array(0, c(N.test,L))
  best.choice.lastobs <- array(0, c(N.test))
  actual.disutility.lastobs <- array(0, c(N.test))
  for (n in 1:N.test) {
    H.n <- H.test[n]
    
    for (l in 0:(L-1)) {
      for (i in 1:I) {
        for (k in 1:I) {
          expected.disutility.lastobs[n,l+1] <- expected.disutility.lastobs[n,l+1] + beta.planning[k] * cthmm.transition.prob.func(Q.in, i, k, l+1, delta.tau.planning) * prob.posterior.lastobs[n,i]
        }
      }
    }
    
    best.choice.lastobs[n] <- which.min(expected.disutility.lastobs[n,]) - 1
    
    actual.disutility.lastobs[n] <- expected.disutility.optimal[n, best.choice.lastobs[n]+1]
    
  }
  
  # calculate the expected disutility, using the entire history for estimating the distribution of z at the planning time
  expected.disutility.history <- array(0, c(N.test,L))
  best.choice.history <- array(0, c(N.test))
  actual.disutility.history <- array(0, c(N.test))
  for (n in 1:N.test) {
    H.n <- H.test[n]
    
    for (l in 0:(L-1)) {
      for (i in 1:I) {
        for (k in 1:I) {
          expected.disutility.history[n,l+1] <- expected.disutility.history[n,l+1] + beta.planning[k] * cthmm.transition.prob.func(Q.in, i, k, l+1, delta.tau.planning) * prob.posterior.history[n,i]
        }
      }
    }
    
    best.choice.history[n] <- which.min(expected.disutility.history[n,]) - 1
    
    actual.disutility.history[n] <- expected.disutility.optimal[n, best.choice.history[n]+1]
    
  }
  
  ret <- list("expected.disutility.optimal" = expected.disutility.optimal, "expected.disutility.lastobs" = expected.disutility.lastobs, "expected.disutility.history" = expected.disutility.history, "best.choice.optimal" = best.choice.optimal, "best.choice.lastobs" = best.choice.lastobs, "best.choice.history" = best.choice.history, "actual.disutility.optimal" = actual.disutility.optimal, "actual.disutility.lastobs" = actual.disutility.lastobs, "actual.disutility.history" = actual.disutility.history) 
  
  return(ret)  
  
}

