
# generate synthetic data
cthmm.ext2.data.generation <- function(I, J, L, pi, lambda, R, Q, mu, eta, eta.prime, N, H.params, between.visit.time.params, accurate.examination.params) {
  
  # generate the number of observations for each patient
  H <- round(rnorm(N, mean=H.params[1], sd=H.params[2]))
  H[ H < 3 ] <- 3
  
  # age, accurate examination indicators, visit times, health state at the visit time, observation, and intervention
  a <- array(numeric(), c(N))
  O <- array(numeric(), c(N, max(H)))
  tau.obs <- array(numeric(), c(N, max(H)))
  z.obs <- array(numeric(), c(N, max(H)))
  z.acc <- array(numeric(), c(N, max(H)))
  y.obs <- array(numeric(), c(N, max(H)))
  u.obs <- array(numeric(), c(N, max(H)))
  
  # true health state change times, and true health state values
  max.num.true.changes <- 4 * max(H) * (between.visit.time.params[1] / min(1/lambda))
  tau.true <- array(numeric(), c(N, max.num.true.changes))
  z.true <- array(numeric(), c(N, max.num.true.changes))
  
  # the visit times
  for (n in 1:N) {
    tau.obs[n,1] <- 0
    for (itr in 1:(H[n]-1)) {
      between.visit.time <- round(max(3, rnorm(1, between.visit.time.params[1], between.visit.time.params[2])))
      # between.visit.time <- between.visit.time.params[1]
      tau.obs[n,itr+1] <- tau.obs[n,itr] + between.visit.time
    }
  }
  
  # generate the data for each patient
  for (n in 1:N) {

    # generate the age data
    age.group <- sum(runif(1)>=0.25)
    a[n] <- abs((1-age.group) * rnorm(1, 25, 5) + age.group * rnorm(1, 50, 10))
    
    # the initial value of the accurate examination indicator
    O[n,1] <- rbinom(1, 1, accurate.examination.params[1])
    
    # generate the data at time zero
    tau.true[n, 1] <- 0
    z.true[n, 1] <- sample(1:I, 1, replace=TRUE, prob=pi)
    z.obs[n, 1] <- z.true[n, 1]
    if (O[n,1] == 0) {
      y.obs[n, 1] <- rbinom(1, J-1, prob=mu[z.obs[n,1], 1])
      u.obs[n, 1] <- rbinom(1, L-1, prob=eta[y.obs[n,1]+1])
    } else {
      z.acc[n, 1] <- z.obs[n, 1]
      u.obs[n, 1] <- rbinom(1, L-1, prob=eta.prime[z.obs[n,1]])
    }
    
    # iterate over the inter-visit intervals
    counter <- 0
    for (itr in 1:(H[n]-1)) {
      
      # the value of the next accurate examination indicator
      if (O[n,itr] == 0) {
        Sigmoid.input <- (y.obs[n,itr] - accurate.examination.params[2]) / accurate.examination.params[3] 
      } else {
        Sigmoid.input <- (z.obs[n,itr] - accurate.examination.params[4]) / accurate.examination.params[5]
      }
      O[n, itr+1] <- rbinom(1, 1, exp(Sigmoid.input)/(1+exp(Sigmoid.input)))
      
      # in each interval, iterate over the changes in the true underlying health state variable
      counter <- counter + 1
      tau.true[n, counter+1] <- tau.true[n, counter] + rexp(1, rate=lambda[u.obs[n, itr]+1, z.true[n, counter]])
      z.true[n, counter+1] <- sample(1:I, 1, replace=TRUE, prob=R[z.true[n, counter], , u.obs[n, itr]+1])
      while (tau.true[n, counter+1] < tau.obs[n,itr+1]) {
        counter <- counter + 1
        tau.true[n, counter+1] <- tau.true[n, counter] + rexp(1, rate=lambda[u.obs[n, itr]+1, z.true[n, counter]])
        z.true[n, counter+1] <- sample(1:I, 1, replace=TRUE, prob=R[z.true[n, counter], , u.obs[n, itr]+1])
      }
      tau.true[n, counter+1] <- tau.obs[n,itr+1]
      z.true[n, counter+1] <- z.true[n, counter]
      
      # set the variables at the visit time based on the true underlying health state at the visit time
      z.obs[n, itr+1] <- z.true[n, counter+1]
      if (O[n,itr+1] == 0) {
        y.obs[n, itr+1] <- rbinom(1, J-1, prob=mu[z.obs[n,itr+1], u.obs[n,itr]+1])
        u.obs[n, itr+1] <- rbinom(1, L-1, prob=eta[y.obs[n,itr+1]+1])
      } else{
        z.acc[n, itr+1] <- z.obs[n, itr+1]
        u.obs[n, itr+1] <- rbinom(1, L-1, prob=eta.prime[z.obs[n,itr+1]])
      }
      
    }
    
  }
  
  ret <- list("H"=H, "a"=a, "O"=O, "tau.true"=tau.true, "tau.obs"=tau.obs,
              "z.true"=z.true, "z.obs"=z.obs, "z.acc"=z.acc, "y.obs"=y.obs, "u.obs"=u.obs)
  
  return(ret) 
  
}

