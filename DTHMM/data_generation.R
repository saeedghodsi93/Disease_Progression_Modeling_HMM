
# generate synthetic data
dthmm.data.generation <- function(I, J, L, pi, Q, mu, eta, N, H.params) {
  
  # generate the number of observations for each patient
  H <- round(rnorm(N, mean=H.params[1], sd=H.params[2]))
  H[ H < 3 ] <- 3
  
  # health state, observation, and intervention
  z <- array(numeric(), c(N,max(H)))
  y <- array(numeric(), c(N,max(H)))
  u <- array(numeric(), c(N,max(H)))
  
  # generate the data for each patient
  for (n in 1:N) {

    H.n <- H[n]
    z[n,1] <- sample(1:I, 1, replace=TRUE, prob=pi)
    y[n,1] <- rbinom(1, J-1, prob=mu[z[n,1]])
    u[n,1] <- rbinom(1, L-1, prob=eta[y[n,1]+1])
    for (t in 2:H.n) {
      z[n,t] <- sample(1:I, 1, replace=TRUE, prob=Q[z[n,t-1],,u[n,t-1]+1])
      y[n,t] <- rbinom(1, J-1, prob=mu[z[n,t]])
      u[n,t] <- rbinom(1, L-1, prob=eta[y[n,t]+1])
    }
    
  }
  
  ret <- list("H"=H, "z"=z, "y"=y, "u"=u)
  
  return(ret) 
  
}

