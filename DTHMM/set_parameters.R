
# set the model parameters
dthmm.set.model.parameters.func <- function(latex.table) {
  
  # the number of possible hidden states, number of possible observations, and number of possible interventions
  I <- 3
  J <- 10
  L <- 3
  
  ##### set the true model parameters
  
  # the initial probability of the hidden state
  pi <- c(0.25, 0.45, 0.30)
  
  # the parameters of the transition probability matrix
  Q <- array(c(0.1, 0.5, 0.4,  0.1, 0.3, 0.6,  0.1, 0.1, 0.8,
               0.2, 0.6, 0.2,  0.2, 0.3, 0.5,  0.1, 0.3, 0.6,
               0.4, 0.5, 0.1,  0.4, 0.3, 0.3,  0.2, 0.4, 0.4), dim=c(I, I, L))
  Q[,,1] <- t(Q[,,1])
  Q[,,2] <- t(Q[,,2])
  Q[,,3] <- t(Q[,,3])
  
  # the parameters of the emission and intervention distributions
  mu <- c(0.1, 0.5, 0.9)
  eta <- c(0.04, 0.15, 0.26, 0.32, 0.43, 0.51, 0.62, 0.77, 0.81, 0.90)
  
  
  ##### set the initial model parameters  
  
  # the initial probability of the hidden state
  pi.init <- c(1/3, 1/3, 1/3)
  
  # the parameters of the transition probability matrix
  Q.init <- array(numeric(), c(I,I,L))
  Q.init[,,1] <- exp(matrix(c(log(0.35), log(0.45), log(0.2),  log(0.3), log(0.3), log(0.4),  log(0.2), log(0.2), log(0.6)), nrow=I, ncol=I, byrow=TRUE))
  Q.init[,,2] <- exp(matrix(c(log(0.45), log(0.4), log(0.15),  log(0.3), log(0.5), log(0.2),  log(0.1), log(0.4), log(0.5)), nrow=I, ncol=I, byrow=TRUE))
  Q.init[,,3] <- exp(matrix(c(log(0.55), log(0.35), log(0.1),  log(0.5), log(0.3), log(0.2),  log(0.2), log(0.5), log(0.3)), nrow=I, ncol=I, byrow=TRUE))
  
  # the parameters of the emission and intervention distributions
  mu.init <- c(0.45, 0.65, 0.7)
  eta.init <- c(0.01, 0.03, 0.06, 0.08, 0.45, 0.48, 0.51, 0.55, 0.91, 0.97)
  
  # print the true and initial parameters in Latex tables
  if (latex.table == TRUE) {
    dthmm.latex.calcs.func(Q, Q.init)
  }
  
  ret <- list("I"=I, "J"=J, "L"=L,
              "pi" = pi, "Q" = Q, "mu" = mu, "eta" = eta,
              "pi.init" = pi.init, "Q.init" = Q.init, "mu.init" = mu.init, "eta.init" = eta.init)
  
  return(ret)
  
}

