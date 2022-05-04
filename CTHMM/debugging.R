
# for debugging
debugging.func <- function(H, a, tau.true, tau.obs, z.true, z.obs, y.obs, u.obs) {
  
  # 1) test the transition probability function formula
  transition.mat <- array(0, c(I,I))
  for (n in 1:N) {
    # tau.true.n <- tau.true[n,]
    # z.begin <- z.true[n,1]
    # z.end <- z.true[n,max(which(tau.true.n<1000))]
    transition.mat[z.begin, z.end] <- transition.mat[z.begin, z.end] + 1
    for (s in 2:H[n]) {
      transition.mat[z.obs[n,s-1], z.obs[n,s]] <- transition.mat[z.obs[n,s-1], z.obs[n,s]] + 1
    }
  }
  for (i in 1:I) {
    temp <- transition.mat[i,] / sum(transition.mat[i,])
    transition.mat[i,] <- temp
  }
  # print(transition.mat)
  # print(expm(10*Q[,,1]))
  
  # 2) test the exmp method for calculating end state conditioned expected number of transitions and sojourn times
  transition.mat <- array(0, c(I,I,I,I))
  counter <- array(0, c(I,I))
  for (n in 1:N) {
    for (s in 2:H[n]) {
      counter[z.obs[n,s-1], z.obs[n,s]] <- counter[z.obs[n,s-1], z.obs[n,s]] + 1
      true.idx.begin <- min(which(tau.true[n,]>=tau.obs[n,s-1]))
      true.idx.end <- max(which(tau.true[n,]<tau.obs[n,s]))
      if (true.idx.begin < true.idx.end) {
        for (true.idx in c(true.idx.begin:(true.idx.end-1))) {
          transition.mat[z.true[n,true.idx], z.true[n,true.idx+1], z.obs[n,s-1], z.obs[n,s]] <- transition.mat[z.true[n,true.idx], z.true[n,true.idx+1], z.obs[n,s-1], z.obs[n,s]] + 1
        }
      }
    }
  }
  for (i.tilde in 1:I) {
    for (k.tilde in 1:I) {
      transition.mat[, , i.tilde, k.tilde] <- transition.mat[, , i.tilde, k.tilde] / counter[i.tilde, k.tilde]
    }
  }
  
  transition.mat.analytical <- array(0, c(I,I,I,I))
  integral.temp <- calc.integral.func(Q, 1, 10)
  for (i.tilde in 1:I) {
    for (k.tilde in 1:I) {
      end.probs <- expm(10*Q[,,1])
      for (i in 1:I) {
        for (k in 1:I) {
          if (k != i) {
            transition.mat.analytical[i, k, i.tilde, k.tilde] <- (Q[i, k, 1] * integral.temp[i, k, i.tilde, k.tilde]) / end.probs[i.tilde, k.tilde]
          }
        }
      }
    }
  }
  
  sojourn.time.analytical <- array(0, c(I,I,I))
  integral.temp <- calc.integral.func(Q, 1, 10)
  for (i.tilde in 1:I) {
    for (k.tilde in 1:I) {
      end.probs <- expm(10*Q[,,1])
      for (i in 1:I) {
        for (k in 1:I) {
          sojourn.time.analytical[i, i.tilde, k.tilde] <- integral.temp[i, i, i.tilde, k.tilde] / end.probs[i.tilde, k.tilde]
        }
      }
    }
  }
  
}

