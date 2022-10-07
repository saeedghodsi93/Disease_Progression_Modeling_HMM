
# calculate the expected disutility for the test patients
planning.expected.disutility.func <- function(N.test, I, J, L, H.test, a.test, z.obs.test, prob.posterior.lastobs, prob.posterior.history, lambda.in, R.in, Q.in, delta.tau.planning, beta.z.planning, beta.u.planning, between.visit.time.params) {
  
  # calculate the actual disutility for each patient and each intervention choice (does not contain the beta.u.planning; add later)
  actual.disutility <- array(0, c(N.test,L))
  generate_next_sample <- FALSE
  if (generate_next_sample == TRUE) {
    for (n in 1:N.test) {
      H.n <- H.test[n]
      
      for (l in 0:(L-1)) {
        
        # true health state change times, and true health state values
        temp.max.num.true.changes <- 5 * (between.visit.time.params[1] / min(1/lambda.in))
        temp.tau.true <- array(numeric(), c(temp.max.num.true.changes))
        temp.z.true <- array(numeric(), c(temp.max.num.true.changes))
        
        # in each interval, iterate over the changes in the true underlying health state variable
        counter <- 0
        temp.tau.true[counter+1] <- rexp(1, rate=lambda.in[l+1, z.obs.test[n, H.n]])
        temp.z.true[counter+1] <- sample(1:I, 1, replace=TRUE, prob=R.in[z.obs.test[n, H.n], , l+1])
        while (temp.tau.true[counter+1] < delta.tau.planning) {
          counter <- counter + 1
          temp.tau.true[counter+1] <- temp.tau.true[counter] + rexp(1, rate=lambda.in[l+1, temp.z.true[counter]])
          temp.z.true[counter+1] <- sample(1:I, 1, replace=TRUE, prob=R.in[temp.z.true[counter], , l+1])
        }
        if (counter > 0) {
          z.end <- temp.z.true[counter]
        } else {
          z.end <- z.obs.test[n, H.n]
        }
        
        actual.disutility[n,l+1] <- beta.z.planning[z.end]
          
      }
    }
  }
  
  # calculate the expected disutility, assuming z is known at the planning time
  expected.disutility.optimal <- array(0, c(N.test,L))
  best.choice.optimal <- array(0, c(N.test))
  actual.disutility.optimal <- array(0, c(N.test))
  for (n in 1:N.test) {
    H.n <- H.test[n]
    
    for (l in 0:(L-1)) {
      for (k in 1:I) {
        expected.disutility.optimal[n,l+1] <- expected.disutility.optimal[n,l+1] + beta.z.planning[k] * cthmm.transition.prob.func(Q.in, z.obs.test[n,H.n], k, l+1, delta.tau.planning)
        # if ((n < 10) & (z.obs.test[n,H.test[n]]==3)) {
        #   print(n)
        #   print(l)
        #   print(k)
        #   print(beta.z.planning[k])
        #   print(cthmm.transition.prob.func(Q.in, z.obs.test[n,H.n], k, l+1, delta.tau.planning))
        # }
      }
      expected.disutility.optimal[n,l+1] <- expected.disutility.optimal[n,l+1] + beta.u.planning[l+1]
    }
    
    best.choice.optimal[n] <- which.min(expected.disutility.optimal[n,]) - 1
    
    actual.disutility.optimal[n] <- expected.disutility.optimal[n, best.choice.optimal[n]+1]
    # actual.disutility.optimal[n] <- actual.disutility[n, best.choice.optimal[n]+1]
    
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
          expected.disutility.lastobs[n,l+1] <- expected.disutility.lastobs[n,l+1] + beta.z.planning[k] * cthmm.transition.prob.func(Q.in, i, k, l+1, delta.tau.planning) * prob.posterior.lastobs[n,i]
        }
      }
      expected.disutility.lastobs[n,l+1] <- expected.disutility.lastobs[n,l+1] + beta.u.planning[l+1]
    }
    
    best.choice.lastobs[n] <- which.min(expected.disutility.lastobs[n,]) - 1
    
    actual.disutility.lastobs[n] <- expected.disutility.optimal[n, best.choice.lastobs[n]+1]
    # actual.disutility.lastobs[n] <- actual.disutility[n, best.choice.lastobs[n]+1]
    
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
          expected.disutility.history[n,l+1] <- expected.disutility.history[n,l+1] + beta.z.planning[k] * cthmm.transition.prob.func(Q.in, i, k, l+1, delta.tau.planning) * prob.posterior.history[n,i]
        }
      }
      expected.disutility.history[n,l+1] <- expected.disutility.history[n,l+1] + beta.u.planning[l+1]
    }
    
    best.choice.history[n] <- which.min(expected.disutility.history[n,]) - 1
    
    actual.disutility.history[n] <- expected.disutility.optimal[n, best.choice.history[n]+1]
    # actual.disutility.history[n] <- actual.disutility[n, best.choice.history[n]+1]
    
  }
  
  ret <- list("expected.disutility.optimal" = expected.disutility.optimal, "expected.disutility.lastobs" = expected.disutility.lastobs, "expected.disutility.history" = expected.disutility.history, "best.choice.optimal" = best.choice.optimal, "best.choice.lastobs" = best.choice.lastobs, "best.choice.history" = best.choice.history, "actual.disutility" = actual.disutility, "actual.disutility.optimal" = actual.disutility.optimal, "actual.disutility.lastobs" = actual.disutility.lastobs, "actual.disutility.history" = actual.disutility.history) 
  
  return(ret)  
  
}

