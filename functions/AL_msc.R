AL_msc <- function(beta, rets, alfa, VaRm, ESm,OUT) {

############################################################################################################################
  # Computes the AL loss or the (VaR, ES) forecast for the relative score combining method by Taylor (2020, IJF)
  # Output: 
  # if OUT=1, returns AL loss value; if OUT=2, returns (VaR, ES) forecast series
  # Inputs:
  # lambda: tuning parameter
  # rets: training data (returns) Tx1 matrix or vector
  # alfa: nominal VaR coverage
  # VaRm: matrix of training VaR forecasts - T x nmod (nmod is the size of the model universe)
  # ESm: matrix of training ES forecasts - T x nmod (nmod is the size of the model universe)
  # FZmat: matrix of joint losses in the training period - T x nmod 
  # OUT: selects the output of the function (see above) 
  ###########################################################################################################################

# Compute the model weights as a function of lambda
T=length(rets)
betaq=beta[1:ncol(VaRm)]
betae=beta[(ncol(VaRm)+1):(2*ncol(VaRm))]
betaq=exp(betaq)/sum(exp(betaq))
betae=exp(betae)/sum(exp(betae))

# Compute the VaR and ES series
vt=rowSums(betaq*VaRm)  
est=vt+rowSums(betae*(ESm-VaRm))  

#Compute the AL loss
 AD <- (rets - vt) * (alfa - (rets <= vt))
 iest <- est^(-1)
 AL <- -log((alfa - 1) * iest) - ((AD / alfa) * iest)
 AL <- mean(AL)

  if (!is.finite(AL)) {
    AL <- 1e+100
  }

  
 # Select the output profile
  if (OUT == 1) {
    output <- AL
  } else if (OUT == 2) {
    output <- cbind(vt, est)
  } else {
    stop("Wrong output selected. Choose OUT = 1 for RQ, or OUT = 2 for [VaR, Hit].")
  }
  
#Return the selected output
  return(AL)
}
