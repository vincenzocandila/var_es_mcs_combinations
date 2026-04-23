###########################################################################################################################
# transbeta: Computes normalized combination weights for the Minimum Score Combining (MSC) method
#            as proposed by Taylor (2020, IJF).
#
# Input:
#   beta : numeric vector (2*nmod x 1) of combination parameters,
#          where nmod is the number of models being combined.
#          The first nmod elements correspond to VaR weights;
#          the remaining nmod elements correspond to ES weights.
#
# Output:
#   A (nmod x 2) matrix of normalized weights, where:
#          column 1 (betaq) = normalized weights for VaR combination
#          column 2 (betae) = normalized weights for ES combination
###########################################################################################################################


transbeta <- function(beta) {
  
  nmod  <- length(beta) / 2          # Number of models (beta contains 2 parameters per model)
  betaq <- beta[1:nmod]              #  VaR combination parameters
  betae <- beta[(nmod + 1):(2*nmod)] #  ES combination parameters
  
  # Transform betaq and betae
  # This ensures that: (i) all weights are strictly positive, and
  #                    (ii) weights within each group sum to one.
  betaq <- exp(betaq) / sum(exp(betaq))
  betae <- exp(betae) / sum(exp(betae))
  
  return(cbind(betaq, betae))
}
