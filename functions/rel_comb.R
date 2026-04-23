rel_comb<- function(fzmat,lambda){
##################################################################################
#This function computes the model weights for the RSC method by Taylor (2020, IJF)
# Output:
# w: vector of model weights
# Input:
# fzmat: nmod x 1 vector of scoring function values
# lambda: tuning parameter (to be estimated)
###################################################################################

n=length(fzmat)
w=exp(-lambda*(fzmat))
w=w/sum(w)
return(w)
}



