
# set the model parameters
planning.set.model.parameters.func <- function(latex.table) {
  
  # the number of possible hidden states, number of possible observations, and number of possible interventions
  I <- 3
  J <- 10
  L <- 3
  
  ##### set the true model parameters
  
  # the parameters of the generator matrix
  lambda <- rbind(c(1/18, 1/5, 1/25),
                  c(1/9, 1/11, 1/10),
                  c(1/4, 1/8, 1/6))
  R <- array(c(0, 0.9, 0.1,  0.2, 0, 0.8,  0.1, 0.9, 0,
               0, 0.7, 0.3,  0.65, 0, 0.35,  0.15, 0.85, 0,
               0, 0.55, 0.45,  0.5, 0, 0.5,  0.6, 0.4, 0), dim=c(I, I, L))
  R[,,1] <- t(R[,,1])
  R[,,2] <- t(R[,,2])
  R[,,3] <- t(R[,,3])
  Q <- array(c(-lambda[1,1], lambda[1,1]*R[1,2,1], lambda[1,1]*R[1,3,1],
               lambda[1,2]*R[2,1,1], -lambda[1,2], lambda[1,2]*R[2,3,1],
               lambda[1,3]*R[3,1,1], lambda[1,3]*R[3,2,1], -lambda[1,3],
               -lambda[2,1], lambda[2,1]*R[1,2,2], lambda[2,1]*R[1,3,2],
               lambda[2,2]*R[2,1,2], -lambda[2,2], lambda[2,2]*R[2,3,2],
               lambda[2,3]*R[3,1,2], lambda[2,3]*R[3,2,2], -lambda[2,3],
               -lambda[3,1], lambda[3,1]*R[1,2,3], lambda[3,1]*R[1,3,3],
               lambda[3,2]*R[2,1,3], -lambda[3,2], lambda[3,2]*R[2,3,3],
               lambda[3,3]*R[3,1,3], lambda[3,3]*R[3,2,3], -lambda[3,3]), dim=c(I, I, L)) 
  Q[,,1] <- t(Q[,,1])
  Q[,,2] <- t(Q[,,2])
  Q[,,3] <- t(Q[,,3])
  
  # the covariate coefficients
  rho <- array(c(0, 0.02, 0.02,  0.015, 0, 0.015,  0.01, 0.01, 0,
                 0, 0.02, 0.02,  0.015, 0, 0.015,  0.01, 0.01, 0,
                 0, 0.02, 0.02,  0.015, 0, 0.015,  0.01, 0.01, 0), dim=c(I, I, L))
  rho[,,1] <- t(rho[,,1])
  rho[,,2] <- t(rho[,,2])
  rho[,,3] <- t(rho[,,3])
  
  # the initial probability of the hidden state (steady-state probability associated with u=0 for a long time)
  pi <- expm(100000*Q[,,1])[1,]
  
  # the parameters of the emission and intervention distributions
  mu <- c(0.45, 0.55, 0.75)
  eta <- c(0.04, 0.11, 0.19, 0.32, 0.41, 0.47, 0.63, 0.77, 0.81, 0.90)
  
  
  ##### set the initial model parameters  
  
  # the initial probability of the hidden state
  pi.init <- c(1/3, 1/3, 1/3)
  
  # the parameters of the generator matrix
  lambda.init <- rbind(c(1/5, 1/8, 1/12),
                        c(1/6, 1/10, 1/10),
                        c(1/9, 1/8, 1/4))
  R.init <- array(c(0, 0.5, 0.5,  0.3, 0, 0.7,  0.2, 0.8, 0,
                     0, 0.6, 0.4,  0.4, 0, 0.6,  0.4, 0.6, 0,
                     0, 0.7, 0.3,  0.6, 0, 0.4,  0.4, 0.6, 0), dim=c(I, I, L))
  R.init[,,1] <- t(R.init[,,1])
  R.init[,,2] <- t(R.init[,,2])
  R.init[,,3] <- t(R.init[,,3])
  Q.init <- array(c(-lambda.init[1,1], lambda.init[1,1]*R.init[1,2,1], lambda.init[1,1]*R.init[1,3,1],
                     lambda.init[1,2]*R.init[2,1,1], -lambda.init[1,2], lambda.init[1,2]*R.init[2,3,1],
                     lambda.init[1,3]*R.init[3,1,1], lambda.init[1,3]*R.init[3,2,1], -lambda.init[1,3],
                     -lambda.init[2,1], lambda.init[2,1]*R.init[1,2,2], lambda.init[2,1]*R.init[1,3,2],
                     lambda.init[2,2]*R.init[2,1,2], -lambda.init[2,2], lambda.init[2,2]*R.init[2,3,2],
                     lambda.init[2,3]*R.init[3,1,2], lambda.init[2,3]*R.init[3,2,2], -lambda.init[2,3],
                     -lambda.init[3,1], lambda.init[3,1]*R.init[1,2,3], lambda.init[3,1]*R.init[1,3,3],
                     lambda.init[3,2]*R.init[2,1,3], -lambda.init[3,2], lambda.init[3,2]*R.init[2,3,3],
                     lambda.init[3,3]*R.init[3,1,3], lambda.init[3,3]*R.init[3,2,3], -lambda.init[3,3]), dim=c(I, I, L)) 
  Q.init[,,1] <- t(Q.init[,,1])
  Q.init[,,2] <- t(Q.init[,,2])
  Q.init[,,3] <- t(Q.init[,,3])
  
  # the parameters of the emission and intervention distributions
  mu.init <- c(0.45, 0.65, 0.7)
  eta.init <- c(0.01, 0.03, 0.06, 0.08, 0.45, 0.48, 0.51, 0.55, 0.91, 0.97)
  
  # print the true and initial parameters in Latex tables
  if (latex.table == TRUE) {
    cthmm.latex.calcs.func(I, L, lambda, R, rho, Q, R.init, Q.init)
  }
  
  ret <- list("I"=I, "J"=J, "L"=L,
              "pi" = pi, "lambda"=lambda, "R"=R, "Q" = Q, "mu" = mu, "eta" = eta,
              "pi.init" = pi.init, "Q.init" = Q.init, "mu.init" = mu.init, "eta.init" = eta.init)
  
  return(ret)
  
}

