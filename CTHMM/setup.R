
# initial setup
setup.func <- function(constant.seed) {
  
  # set the seed
  if (constant.seed) {
    set.seed(1834772)
  }
  
  # set options
  options(width = 80)
  
  # set the conda environment for reticulate
  use_condaenv("base")
  
}


# progress bar
progress.bar <- function(counter.in, counter.max.in, flag, debugging.mode) {
  
  # in the debugging mode, do not run the progress bar
  if (debugging.mode == FALSE) {
    
    if (flag==TRUE) {
      extra <- nchar('||100%')
      width <- options()$width
      step <- round(counter.in / counter.max.in * (width - extra))
      text <- sprintf('|%s%s|% 3s%%', strrep('=', step), strrep(' ', width - step - extra), round(counter.in / counter.max.in * 100))
      cat(text) 
    }
    else{
      cat(if (counter.in == counter.max.in) '\n' else '\014')
    }
    
  }
  
}


# print the input matrix as Latex table
latex.table.func <- function(mat.in) {
  
  x = xtable(mat.in, align=rep("",ncol(mat.in)+1))
  print(x, floating=FALSE, tabular.environment="bmatrix", hline.after=NULL, include.rownames=FALSE, include.colnames=FALSE)
  
}


# prepare the matrices for usage in Latex
dthmm.latex.calcs.func <- function(Q, Q.init) {
  
  # print the matrices into Latex tables
  print("Q: ")
  latex.table.func(Q[,,1])
  latex.table.func(Q[,,2])
  latex.table.func(Q[,,3])
  print("Q.init: ")
  latex.table.func(Q.init[,,1])
  latex.table.func(Q.init[,,2])
  latex.table.func(Q.init[,,3])
  
}


# prepare the matrices for usage in Latex
cthmm.latex.calcs.func <- function(I, L, lambda, R, rho, Q, R.init, Q.init) {

  # print the matrices into Latex tables
  print("R: ")
  latex.table.func(R[,,1])
  latex.table.func(R[,,2])
  latex.table.func(R[,,3])
  print("Q: ")
  latex.table.func(Q[,,1])
  latex.table.func(Q[,,2])
  latex.table.func(Q[,,3])
  print("R.init: ")
  latex.table.func(R.init[,,1])
  latex.table.func(R.init[,,2])
  latex.table.func(R.init[,,3])
  print("Q.init: ")
  latex.table.func(Q.init[,,1])
  latex.table.func(Q.init[,,2])
  latex.table.func(Q.init[,,3])
  
  # calculate Q for the example case
  a <- 50
  Q.example <- array(0, dim=c(I,I,L))
  for (l in 0:(L-1)) {
    for (i in 1:I){
      for (k in 1:I){
        if (k != i){
          Q.example[i,k,l+1] <- lambda[l+1,i] * R[i,k,l+1] * exp(a * rho[i,k,l+1])
        }
      }
      temp <- -sum(Q.example[i,,l+1])
      Q.example[i,i,l+1] <- temp
    }
  }
  print("Q.example: ")
  latex.table.func(Q.example[,,1])
  latex.table.func(Q.example[,,2])
  latex.table.func(Q.example[,,3])
  
}
