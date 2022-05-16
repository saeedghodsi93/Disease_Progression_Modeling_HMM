

# the Gibbs sampling approach for calculating the end-state posterior probabilities
dthmm.gibbs.sampling.func <- function(I, J, L, y, u, pi.tilde, Q.tilde, mu.tilde, eta.tilde, H.n, n, M, B) {
  
  # initialize the sample
  z.current <- sample(1:I, H.n, replace=TRUE, prob=pi.tilde)
  
  # throw away the first B samples and keep the next M samples
  z.samples <- array(numeric(), c(M,H.n))
  counter.MC <- 1
  while (counter.MC <= B + M) {
    
    for (t in 1:H.n) {
      
      # calculate the conditional probabilities (intervention terms cancel out in normalization)
      z.prob <- array(numeric(), c(I))
      temp <- array(numeric(), c(I))
      for (i in 1:I) {
        if (t == 1) {
          temp[i] <- Q.tilde[i, z.current[2], u[n,1]+1] * dthmm.emission.prob.func(J, mu.tilde[i], y[n,1]) * pi.tilde[i] 
        }
        else if (t == H.n) {
          temp[i] <- Q.tilde[z.current[H.n-1], i, u[n,H.n-1]+1] * dthmm.emission.prob.func(J, mu.tilde[i], y[n,H.n])
        } else {
          temp[i] <- Q.tilde[z.current[t-1], i, u[n,t-1]+1] * Q.tilde[i, z.current[t+1], u[n,t]+1] * dthmm.emission.prob.func(J, mu.tilde[i], y[n,t])
        }
      }
      
      # normalize the probabilities 
      if(is.na(sum(temp))) {
        for (i in 1:I) {
          z.prob[i] <- 1 / I
        }  
      } else if (sum(temp) == 0) {
        for (i in 1:I) {
          z.prob[i] <- 1 / I
        }
      } else {
        for (i in 1:I) {
          z.prob[i] <- temp[i] / sum(temp)
        }
      }
        
      # draw one sample
      z.current[t] <- sample(1:I, 1, replace=TRUE, prob=z.prob)
      
    }
    
    # keep the sample if the burn-in period has finished
    if (counter.MC > B) {
      z.samples[counter.MC-B,] <- z.current
    }
    
    # increase the counter
    counter.MC  <- counter.MC + 1
  }
  
  # approximate the posterior probabilities of z using the generated samples
  gamma <- array(0, c(H.n,I))
  for (t in 1:H.n) {
    for (i in 1:I) {
      for (m in 1:M) {
        gamma[t,i] <- gamma[t,i] + sum(z.samples[m,t]==i) 
      }
      gamma[t,i] <- gamma[t,i] / M
    }
  }
  nu <- array(0, c(H.n,I,I))
  for (t in 1:(H.n-1)) {
    for (i in 1:I) {
      for (k in 1:I) {
        for (m in 1:M) {
          nu[t,i,k] <- nu[t,i,k] + sum(z.samples[m,t]==i) * sum(z.samples[m,t+1]==k) 
        }
        nu[t,i,k] <- nu[t,i,k] / M
      }
    }
  }
  
  ret <- list("gamma" = gamma, "nu" = nu)
  
  return(ret)
  
}

