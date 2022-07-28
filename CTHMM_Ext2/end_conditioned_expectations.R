
# calculate the integral using the expm method
cthmm.ext2.calc.integral.func <- function(I, Q.in, l.in, tau.in) {
  
  Q.temp <- Q.in[,, l.in]
  integral <- array(numeric(), dim=c(I,I,I,I))
  for (i in 1:I) {
    for (k in 1:I) {
      B.temp <- array(0, dim=c(I,I))
      B.temp[i,k] <- 1
      integral[i,k,,] <- expm(tau.in * rbind(cbind(Q.temp, B.temp), cbind(array(0,dim=c(I,I)), Q.temp)))[1:I, (I+1):(2*I)]
    }
  }
  
  ret <- integral
  
  return(ret)
  
}


# calculate the expected number of transitions and expected sojourn times between to consecutive visits
cthmm.ext2.end.conditioned.expectations.func <- function(I, Q.in, l.in, tau.in) {
  
  # evaluate the integral
  integral.mat <- cthmm.ext2.calc.integral.func(I, Q.in, l.in, tau.in)
  
  # calculate the end state probabilities
  end.state.probs <- expm(tau.in * Q.in[,,l.in])
  
  # calculate the expected number of transitions
  transition.mat <- array(0, c(I,I,I,I))
  for (i.tilde in 1:I) {
    for (k.tilde in 1:I) {
      for (i in 1:I) {
        for (k in 1:I) {
          if (k != i) {
            transition.mat[i, k, i.tilde, k.tilde] <- (Q.in[i, k, l.in] * integral.mat[i, k, i.tilde, k.tilde]) / end.state.probs[i.tilde, k.tilde]
          }
        }
      }
    }
  }
  
  # calculate the expected sojourn times
  sojourn.time <- array(0, c(I,I,I))
  for (i.tilde in 1:I) {
    for (k.tilde in 1:I) {
      end.state.probs <- expm(tau.in * Q.in[,,l.in])
      for (i in 1:I) {
        for (k in 1:I) {
          sojourn.time[i, i.tilde, k.tilde] <- integral.mat[i, i, i.tilde, k.tilde] / end.state.probs[i.tilde, k.tilde]
        }
      }
    }
  }
  
  ret <- list("transition.mat" = transition.mat, "sojourn.time" = sojourn.time)
  
  return(ret)
  
}

