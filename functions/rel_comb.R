rel_comb<- function(fzmat,lambda){
##################################################################################
#This function computes the model weights for the RSC method by Taylor (2020, IJF)
# Input:
# fzmat: nmod x 1 vector of scoring function values
# lambda: tuning parameter (to be estimated)
# Output:
# w: vector of model weights
###################################################################################

fzmat <- as.numeric(fzmat)
  n <- length(fzmat)

  # Numerical stability trick (log-sum-exp):
  # subtracting the minimum value does not affect the relative weights,
  # but prevents underflow in the exponential function
  fzmat_shift <- fzmat - min(fzmat)

  w <- exp(-lambda * fzmat_shift)
  s <- sum(w)

  if (!is.finite(s) || s == 0) {
    warning("rel_comb: weight sum is zero or non-finite; using uniform weights.")
    return(rep(1/n, n))
  }

  w <- w / s
  return(w)
}



