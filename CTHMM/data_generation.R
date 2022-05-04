
# generate synthetic data
data.generation <- function(I, J, L, pi, lambda, R, Q, mu, eta, N, H.params, between.visit.time.params) {
  
  # generate the number of observations for each patient
  H <- round(rnorm(N, mean=H.params[1], sd=H.params[2]))
  H[ H < 3 ] <- 3
  
  # age, visit times, health state at the visit time, observation, and intervention
  a <- array(numeric(), c(N))
  tau.obs <- array(numeric(), c(N, max(H)))
  z.obs <- array(numeric(), c(N, max(H)))
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
    
    # generate the data at time zero
    tau.true[n, 1] <- 0
    z.true[n, 1] <- sample(1:I, 1, replace=TRUE, prob=pi)
    z.obs[n, 1] <- z.true[n, 1]
    y.obs[n, 1] <- rbinom(1, J-1, prob=mu[z.obs[n,1]])
    u.obs[n, 1] <- rbinom(1, L-1, prob=eta[y.obs[n,1]+1])
    # u.obs[n, 1] <- 0

    # iterate over the inter-visit intervals
    counter <- 0
    for (itr in 1:(H[n]-1)) {
      
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
      y.obs[n, itr+1] <- rbinom(1, J-1, prob=mu[z.obs[n,itr+1]])
      u.obs[n, itr+1] <- rbinom(1, L-1, prob=eta[y.obs[n,itr+1]+1])
      # u.obs[n, itr+1] <- 0
      
    }
    
  }
  
  ret <- list("H"=H, "a"=a, "tau.true"=tau.true, "tau.obs"=tau.obs,
              "z.true"=z.true, "z.obs"=z.obs, "y.obs"=y.obs, "u.obs"=u.obs)
  
  return(ret) 
  
}

