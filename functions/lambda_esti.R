lambda_esti <- function(rets, alfa, VaRm, ESm, FZmat) {
  #################################################################################################
  # Estimates the tuning parameter lambda of the RSC combining method by Taylor (2020, IJF)
  # Output: 
  # list: *$esti= estimated lambda; *$AL: minimized AL score
  # Inputs:
  # rets: training data (returns)
  # alfa: nominal VaR coverage
  # VaRm: matrix of training VaR forecasts - T x nmod (nmod is the size of the model universe)
  # ESm: matrix of training ES forecasts - T x nmod (nmod is the size of the model universe)
  # FZmat: matrix of joint losses in the training period - T x nmod 
  #################################################################################################
  
  # Initialization and settings
  N <- length(rets)
  OUT <- 1

  set.seed(123)  # Set random seed for replicability of results
  ninit <- 1000  # Number of initial samples
  nelect <- 5    # Set of best solutions
  nrep <- 5      # Number of iterations in estimation loops
  
  # Compute AL Loss over a random set of initial values
  fstore <- matrix(0, ninit, 2)
  out <- 1

  fstore[,1] <- runif(ninit) * (0.4)
  
  temp <- as.matrix(fstore[,1], ncol = 1)
  
  fstore[,2] <- sapply(temp[,1], function(x) AL_lambda(x, rets, alfa, VaRm, ESm, FZmat, out))
  
  # Take best nelect solutions
  ranksol <- fstore[order(fstore[, 2]), ]
  startmat <- ranksol[1:nelect, 1]
  
  # Use each "best" solution as initial point for iterating between simplex
  # and quasi newton optimizers until convergence
  ffstore <- matrix(0, nelect, 2)
  exitflag <- integer(nelect)
  
  for (jj in 1:nelect) {
    fun <- function(x) AL_lambda(x, rets, alfa, VaRm, ESm, FZmat, out)
    
    # First optimization using (simplex method)
    result <- optim(startmat[jj], fun, method = "Nelder-Mead", 
                    control = list(warn.1d.NelderMead = FALSE))
    
    ffstore[jj, 1] <- result$par
    ffstore[jj, 2] <- result$value
    
    # Iteration loop for quasi newton and simplex optimization
    for (ii in 1:nrep) {  
      result <- optim(ffstore[jj, 1], fun, method = "BFGS")
      ffstore[jj, 1] <- result$par
      ffstore[jj, 2] <- result$value
      
      result2 <- optim(ffstore[jj, 1], fun, method = "Nelder-Mead", 
                       control = list(warn.1d.NelderMead = FALSE))
      ffstore[jj, 1] <- result2$par
      ffstore[jj, 2] <- result2$value
      
      if (result$convergence == 0) {  # Check if convergence is achieved
        exitflag[jj] <- 1
        break
      }
    }
  }
  
  # Select best estimate and report output
  index <- which(ffstore[,2]<1e+100)
  ffstore <- ffstore[index,]
  ffstore <- ffstore[order(ffstore[, 2]), ]
  esti <- ffstore[1, 1]
  AL <- ffstore[1, 2]
  
  return(list(esti = esti, AL = AL))
}