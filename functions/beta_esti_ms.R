beta_esti_ms <- function(rets, alfa,VaRm,ESm) {
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
  nelect <-5     # Set of best solutions
  nrep <- 5      # Number of iterations in estimation loops
  
  nmod=dim(VaRm)[2]
# Compute Tick Loss over a random set of initial values
  fstore <- matrix(0, ninit,((2*nmod)+1))
  out=1
  temp <- runif(ninit*2*nmod)
  fstore[,1:((2*nmod))]=matrix(temp,nrow=ninit,ncol=2*nmod)    

temp=as.matrix(fstore[,1:(2*nmod)],ncol=nmod)
last=2*nmod+1 
fstore[,last] <- apply(temp,1,'AL_msc',rets, alfa, VaRm, ESm,out)


# Take best nselect solutions
  ranksol <- fstore[order(fstore[,last]), ]
  startmat <- ranksol[1:nelect, -last]
  
# Use each "best" solution as initial point for iterating between simplex
# and quasi newton optimizers until convergence
  ffstore <- matrix(0, nelect, 2*nmod+1)
  exitflag <- integer(nelect)
  
  for (jj in 1:nelect) {
    fun <- function(x) AL_msc(x, rets, alfa, VaRm, ESm,out)
    

# First optimization using (simplex method)
result <- optim(startmat[jj,], fun, method = "Nelder-Mead",control=list(warn.1d.NelderMead="FALSE"))#control = list(maxit = 2000), silent = TRUE)
    
      ffstore[jj, -last] <- result$par
      ffstore[jj, last] <- result$value
    
# Iteration loop for quasi newton and simplex optimization
    for (jj in 1:nrep) {
      result <- optim(ffstore[jj, 1:2*nmod], fun, method = "BFGS")
      ffstore[jj, -last] <- result$par
      ffstore[jj, last] <- result$value
      
      result2 <- optim(ffstore[jj, 1:2*nmod], fun, method = "Nelder-Mead",control=list(warn.1d.NelderMead="FALSE"))
      ffstore[jj, -last] <- result2$par
      ffstore[jj, last] <- result2$value
      
      if (result$convergence == 0) {  # Check if convergence is achieved
        exitflag[jj] <- 1
        break
      }
    }
  }
  

  
# Select best estimate and report output
  index=which(ffstore[,last]<1e+100)
  ffstore=ffstore[index,]
  ffstore <- ffstore[order(ffstore[,last]), ]
  esti <- ffstore[1, -last]
  AL <- ffstore[1, last]
  return(list(esti = esti, AL = AL)) 
}
