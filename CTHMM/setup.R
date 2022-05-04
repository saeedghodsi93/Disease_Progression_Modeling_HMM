
# initial setup
setup.func <- function() {
  
  # add some libraries
  library(MASS)
  library(tidyverse)
  library(tinytex)
  library(ggplot2)
  library(plotly) 
  library(reshape2)
  library(graphics)
  library(AER)
  library(expm)
  library(msm)
  library(cthmm)
  library(hmm.discnp)
  
  # set the seed
  set.seed(1834772)
  
  # set options
  options(width = 80)
  
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

