
# set the model parameters
set.model.parameters.func <- function() {
  
  # the number of possible hidden states, number of possible observations, and number of possible interventions
  I <- 3
  J <- 10
  L <- 3
  
  ##### set the true model parameters
  
  # the initial probability of the hidden state
  pi <- c(0.25, 0.45, 0.3)
  
  # the parameters of the generator matrix
  lambda <- rbind(c(1/4, 1/6, 1/20),
                  c(1/9, 1/12, 1/10),
                  c(1/11, 1/10, 1/7))
  R <- array(c(0, 0.6, 0.4,  0.2, 0, 0.8,  0.1, 0.9, 0,
               0, 0.8, 0.2,  0.6, 0, 0.4,  0.3, 0.7, 0,
               0, 0.9, 0.1,  0.8, 0, 0.2,  0.5, 0.5, 0), dim=c(I, I, L))
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
  
  # the parameters of the emission and intervention distributions
  mu <- matrix(c(0.1, 0.5, 0.9), nrow=I, ncol=1, byrow=TRUE)
  eta <- matrix(c(0.04, 0.15, 0.26, 0.32, 0.43, 0.51, 0.62, 0.77, 0.81, 0.90), nrow=J, ncol=1, byrow=TRUE)
  
  
  ##### set the initial model parameters  
  
  # the initial probability of the hidden state
  pi.tilde <- c(1/3, 1/3, 1/3)
  
  # the parameters of the generator matrix
  lambda.tilde <- rbind(c(1/5, 1/8, 1/12),
                        c(1/6, 1/10, 1/10),
                        c(1/9, 1/8, 1/4))
  R.tilde <- array(c(0, 0.5, 0.5,  0.3, 0, 0.7,  0.2, 0.8, 0,
                     0, 0.6, 0.4,  0.4, 0, 0.6,  0.4, 0.6, 0,
                     0, 0.7, 0.3,  0.6, 0, 0.4,  0.4, 0.6, 0), dim=c(I, I, L))
  R.tilde[,,1] <- t(R.tilde[,,1])
  R.tilde[,,2] <- t(R.tilde[,,2])
  R.tilde[,,3] <- t(R.tilde[,,3])
  Q.tilde <- array(c(-lambda.tilde[1,1], lambda.tilde[1,1]*R.tilde[1,2,1], lambda.tilde[1,1]*R.tilde[1,3,1],
                     lambda.tilde[1,2]*R.tilde[2,1,1], -lambda.tilde[1,2], lambda.tilde[1,2]*R.tilde[2,3,1],
                     lambda.tilde[1,3]*R.tilde[3,1,1], lambda.tilde[1,3]*R.tilde[3,2,1], -lambda.tilde[1,3],
                     -lambda.tilde[2,1], lambda.tilde[2,1]*R.tilde[1,2,2], lambda.tilde[2,1]*R.tilde[1,3,2],
                     lambda.tilde[2,2]*R.tilde[2,1,2], -lambda.tilde[2,2], lambda.tilde[2,2]*R.tilde[2,3,2],
                     lambda.tilde[2,3]*R.tilde[3,1,2], lambda.tilde[2,3]*R.tilde[3,2,2], -lambda.tilde[2,3],
                     -lambda.tilde[3,1], lambda.tilde[3,1]*R.tilde[1,2,3], lambda.tilde[3,1]*R.tilde[1,3,3],
                     lambda.tilde[3,2]*R.tilde[2,1,3], -lambda.tilde[3,2], lambda.tilde[3,2]*R.tilde[2,3,3],
                     lambda.tilde[3,3]*R.tilde[3,1,3], lambda.tilde[3,3]*R.tilde[3,2,3], -lambda.tilde[3,3]), dim=c(I, I, L)) 
  Q.tilde[,,1] <- t(Q.tilde[,,1])
  Q.tilde[,,2] <- t(Q.tilde[,,2])
  Q.tilde[,,3] <- t(Q.tilde[,,3])
  
  # the parameters of the emission and intervention distributions
  mu.tilde <- matrix(c(0.45, 0.65, 0.7), nrow=I, ncol=1, byrow=TRUE)
  eta.tilde <- matrix(c(0.01, 0.03, 0.06, 0.08, 0.45, 0.48, 0.51, 0.55, 0.91, 0.97), nrow=J, ncol=1, byrow=TRUE)
  
  ret <- list("I"=I, "J"=J, "L"=L,
              "pi" = pi, "lambda"=lambda, "R"=R, "Q" = Q, "mu" = mu, "eta" = eta,
              "pi.tilde" = pi.tilde, "Q.tilde" = Q.tilde, "mu.tilde" = mu.tilde, "eta.tilde" = eta.tilde)
  
  return(ret)
  
}

