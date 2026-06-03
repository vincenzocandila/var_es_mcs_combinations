############################################################################################################################
### This code generates the six forecast combinations proposed in the paper:                                             ###
### 'Combining Value-at-Risk and Expected Shortfall Forecasts via the Model Confidence Set'                              ###
############################################################################################################################

#### Analysis performed using R version 4.4.3

#####################################################################################
# 1. Load libraries and functions
#####################################################################################

library(xts)			# Version used: 0.14.1
library(rugarch)		# Version used: 1.5.3
library(np)				# Version used: 0.60.18
library(esback)			# Version used: 0.3.1
library(GAS)			# Version used: 0.3.4.1 
						# Currently, the GAS package has been removed from the CRAN. But previous versions
						# are available at https://cran.r-project.org/src/contrib/Archive/GAS/
library(zoo)			# Version used: 1.8.13

### To install GAS:

#install.packages("https://cran.r-project.org/src/contrib/Archive/GAS/GAS_0.3.4.1.tar.gz",
#                 repos = NULL,
#                 type = "source")

### Load user-defined functions

source('functions/AL_msc.R')
source('functions/transbeta.R')   
source('functions/beta_esti_ms.R')

source('functions/rel_comb.R')
source('functions/lambda_esti.R')
source('functions/AL_lambda.R')

source('functions/main_functions.R')


#####################################################################################
# 2. Load data
#####################################################################################

# Load the dataset containing all inputs required for the six proposed forecast combinations
# The user can choose sp500_tau_0.025_files or shanghai_comp_tau_0.025_files 

load("data/intermediary/sp500_tau_0.025_files.RData")  

#####################################################################################
#### Variables stored in sp500_tau_0.025_files and shanghai_comp_tau_0.025_files:
##   Tin					 : length of the training period
##   nstep				 : number of forecasting steps (model re-estimations) and out-of-sample period length
##   N_model				 : number of competing models (nmod)
##   tau					 : coverage level used	
##   VaR_training_data_mod : array of VaR forecasts over the training period  - Tin x nmod x nstep
##   ES_training_data_mod  : array of ES forecasts over the training period   - Tin x nmod x nstep
##   r_t_in_s_matrix       : matrix of returns over the training period       - Tin x nstep
##   VaR_oos               : matrix of out-of-sample VaR forecasts            - nstep x nmod
##   ES_oos                : matrix of out-of-sample ES forecasts             - nstep x nmod
##   r_t_oos_full          : matrix of out-of-sample returns (xts/zoo format) - nstep x 1
##   r_t_oos               : matrix of out-of-sample returns                  - nstep x 1
##   list_of_models        : character vector of labels for the nmod competing models
#####################################################################################

# Output file name for saving results
filename <- "data/results/sp500_partial_results.RData" # or "data/results/shanghai_comp_partial_results.RData"

#####################################################################################
# 4. Set parameters for MCS and forecast combination methods
#####################################################################################

### Significance level for MCS tests
alpha <- 0.25

### Number of bootstrap replicates for the MCS procedure
B <- 1000

### Decay parameter for the EWMA weighting scheme
theta <- 0.06

### Total number of competing models
M <- N_model 


#####################################################################################
# 5. Initialize output objects
#####################################################################################

### Vectors of length nstep to store VaR and ES forecasts from the six combination methods
### and two benchmarks (mean andmedian combinations)

Mean_Comb_VaR   <- Mean_Comb_ES   <-
  Median_Comb_VaR <- Median_Comb_ES <-
  MCS_Comb_VaR    <- MCS_Comb_ES    <-
  WL_MCS_Comb_VaR <- WL_MCS_Comb_ES <-
  MCS_RSC_Comb_VaR    <- MCS_RSC_Comb_ES    <-
  WL_MCS_RSC_Comb_VaR <- WL_MCS_RSC_Comb_ES <-
  MCS_MSC_Comb_VaR    <- MCS_MSC_Comb_ES    <-
  WL_MCS_MSC_Comb_VaR <- WL_MCS_MSC_Comb_ES <-
  rep(NA, nstep)

### Extend model labels to include the eight combination method names
lab_full <- c(list_of_models,
              "EW-Comb",
              "Median-Comb",
              "MCS-Comb",
              "WL-MCS-Comb",
              "MCS-RS-Comb",
              "WL-MCS-RS-Comb",
              "MCS-MS-Comb",
              "WL-MCS-MS-Comb")

### Initialize lists for storing intermediate MCS results and combination diagnostics
MCS_est_included   <- MCS_est_included_w <- list()

RSC_lambda_store    <- RSC_AL_store    <-
WL_RSC_lambda_store <- WL_RSC_AL_store <-
MSC_AL_store        <- WL_MSC_AL_store <- list()

### Inizialize array for the backtesting

Backtesting_pvalues<- array(NA,dim=c(N_model,6,nstep))


#####################################################################################
# 6. Main loop: compute backtesting pvalues, MCS and forecast combinations at each step 
#####################################################################################

for (tt in 1:nstep) {

  # Compute the FZ loss for each model over the training period
  db_loss_training <- matrix(NA, nrow = Tin, ncol = M)
  
  for (i in 1:M) {
    db_loss_training[, i] <- FZLoss(r_t_in_s_matrix[, tt],
                                    VaR_training_data_mod[, i, tt],
                                    ES_training_data_mod[, i, tt],
                                    tau)
  }
  
  colnames(db_loss_training) <- list_of_models
  
  # Replace any NA values in the loss matrix with column means
  db_loss_training <- apply(db_loss_training, 2, function(col) {
    ifelse(is.na(col), mean(col, na.rm = TRUE), col)
  })

  # store daily-log returns, VaR and ES for the analysis
  r_s   <- r_t_in_s_matrix[, tt]
  VaR_s <- VaR_training_data_mod[, , tt]   # VaR_training_data_mod  #VaR_in_s_array
  ES_s  <- ES_training_data_mod[, , tt]    # ES_training_data_mod #ES_in_s_array

  # ----------------------------------------------------------------------------------
  # Backtesting 
  # ----------------------------------------------------------------------------------

  for (i in seq_len(N_model)) {
    v <- VaR_s[, i]
    e <- ES_s[, i]

    bt1 <- BacktestVaR(r_s, v, tau)
    bt3 <- BacktestVaR(r_s, v, tau)

    Backtesting_pvalues[i, , tt] <- round(c(
      as.numeric(bt1$LRuc[2]),
      as.numeric(bt1$LRcc[2]),
      as.numeric(bt3$DQ[2]),
      esr_backtest(r_s, v, e, tau, version = 1)$pvalue_twosided_asymptotic,
      esr_backtest(r_s, v, e, tau, version = 2)$pvalue_twosided_asymptotic,
      esr_backtest(r_s, v, e, tau, version = 3)$pvalue_twosided_asymptotic
    ), 3)
  }
  

  # ----------------------------------------------------------------------------------
  # Benchmark 1: EW-Comb — Equal-weighted combination across all competing models
  # ----------------------------------------------------------------------------------
  
  Mean_Comb_VaR[tt] <- mean(VaR_oos[tt, ])
  Mean_Comb_ES[tt]  <- mean(ES_oos[tt, ])
  
  
  # ----------------------------------------------------------------------------------
  # Benchmark 2: Median-Comb — Median combination across all competing models
  # ----------------------------------------------------------------------------------
  
  Median_Comb_VaR[tt] <- median(VaR_oos[tt, ])
  Median_Comb_ES[tt]  <- median(ES_oos[tt, ])

  
  # ----------------------------------------------------------------------------------
  # Method 1: MCS-Comb — Equal-weighted combination of MCS-selected models (unweighted FZ loss)
  # ----------------------------------------------------------------------------------
  
  MCS_est <- MCS_f(db_loss_training, B = B, alpha = alpha)
  MCS_est_included[[tt]] <- MCS_est$includedR
  
  # Average VaR and ES over MCS-selected models
  MCS_Comb_VaR[tt] <- mean(VaR_oos[tt, MCS_est$includedR])
  MCS_Comb_ES[tt]  <- mean(ES_oos[tt,  MCS_est$includedR])
  
  
  # ----------------------------------------------------------------------------------
  # Method 2: WL-MCS-Comb — Equal-weighted combination of MCS-selected models (EWMA-weighted FZ loss)
  # ----------------------------------------------------------------------------------
  
  # Apply EWMA weighting to the loss matrix before running MCS
  db_loss_training_w <- ewma_v(theta, db_loss_training)
  MCS_est_w <- MCS_f(db_loss_training_w, B = B, alpha = alpha)
  MCS_est_included_w[[tt]] <- MCS_est_w$includedR
  
  # Average VaR and ES over MCS-selected models (based on EWMA-weighted losses)
  WL_MCS_Comb_VaR[tt] <- mean(VaR_oos[tt, MCS_est_w$includedR])
  WL_MCS_Comb_ES[tt]  <- mean(ES_oos[tt,  MCS_est_w$includedR])
  
  
  # ----------------------------------------------------------------------------------
  # Method 3: MCS-RS-Comb — Relative Score Combination (RSC) restricted to MCS-selected models
  #           (unweighted FZ loss; RSC method as in Taylor, 2020)
  # ----------------------------------------------------------------------------------
  
  MCS_models   <- MCS_est$includedR
  n_mcs_models <- length(MCS_models)
  
  if (n_mcs_models == 0) {
    # No models survive MCS: assign NA
    MCS_RSC_Comb_VaR[tt] <- NA
    MCS_RSC_Comb_ES[tt]  <- NA
    RSC_lambda_store[tt] <- NA
    RSC_AL_store[tt]     <- NA
    
  } else if (n_mcs_models == 1) {
    # Only one model survives MCS: use it directly
    MCS_RSC_Comb_VaR[tt] <- VaR_oos[tt, MCS_models]
    MCS_RSC_Comb_ES[tt]  <- ES_oos[tt,  MCS_models]
    RSC_lambda_store[tt] <- NA
    RSC_AL_store[tt]     <- NA
    
  } else {
    # Extract training-period forecasts and losses for MCS-selected models
    VaRm_mcs  <- as.matrix(VaR_training_data_mod[, MCS_models, tt])
    ESm_mcs   <- as.matrix(ES_training_data_mod[,  MCS_models, tt])
    FZmat_mcs <- as.matrix(db_loss_training[,       MCS_models])
    
    VaRm <- VaRm_mcs
    ESm  <- ESm_mcs
    rets <- r_t_in_s_matrix[, tt]
    nmod <- n_mcs_models
    
    # Build cumulative FZ loss matrix (row-constant, equal to column sums)
    FZmat <- matrix(colSums(FZmat_mcs),
                    nrow = nrow(FZmat_mcs),
                    ncol = ncol(FZmat_mcs),
                    byrow = TRUE)
    
    # Estimate the RSC lambda parameter
    temp <- lambda_esti(rets, tau, VaRm, ESm, FZmat)
    RSC_lambda_store[[tt]] <- lambda <- temp$esti
    RSC_AL_store[[tt]]     <- AL     <- temp$AL
    
    # Compute total FZ loss over the training period for each model
    FZmat_last <- matrix(rep(NA, nmod), ncol = nmod)
    for (j in 1:nmod) {
      FZmat_last[, j] <- sum(FZLoss(rets, VaRm[, j], ESm[, j], tau))
    }
    
    # Derive model weights and compute combined VaR and ES forecasts
    varf  <- VaR_oos[tt, MCS_models]
    esf   <- ES_oos[tt,  MCS_models]
    w     <- rel_comb(FZmat_last, lambda)  # RSC model weights
    cvarf <- sum(w * varf)
    cesf  <- sum(w * esf)
    
    MCS_RSC_Comb_VaR[tt] <- cvarf
    MCS_RSC_Comb_ES[tt]  <- cesf
  }
  
  
  # ----------------------------------------------------------------------------------
  # Method 4: WL-MCS-RS-Comb — Relative Score Combination (RSC) restricted to MCS-selected models
  #           (EWMA-weighted FZ loss; RSC method as in Taylor, 2020)
  # ----------------------------------------------------------------------------------
  
  MCS_models   <- MCS_est_w$includedR
  n_mcs_models <- length(MCS_models)
  
  if (n_mcs_models == 0) {
    # No models survive MCS: assign NA
    WL_MCS_RSC_Comb_VaR[tt] <- NA
    WL_MCS_RSC_Comb_ES[tt]  <- NA
    WL_RSC_lambda_store[tt] <- NA
    WL_RSC_AL_store[tt]     <- NA
    
  } else if (n_mcs_models == 1) {
    # Only one model survives MCS: use it directly
    WL_MCS_RSC_Comb_VaR[tt] <- VaR_oos[tt, MCS_models]
    WL_MCS_RSC_Comb_ES[tt]  <- ES_oos[tt,  MCS_models]
    WL_RSC_lambda_store[tt] <- NA
    WL_RSC_AL_store[tt]     <- NA
    
  } else {
    # Extract training-period forecasts and losses for MCS-selected models
    VaRm_mcs  <- as.matrix(VaR_training_data_mod[, MCS_models, tt])
    ESm_mcs   <- as.matrix(ES_training_data_mod[,  MCS_models, tt])
    FZmat_mcs <- as.matrix(db_loss_training[,       MCS_models])
    
    VaRm <- VaRm_mcs
    ESm  <- ESm_mcs
    rets <- r_t_in_s_matrix[, tt]
    nmod <- n_mcs_models
    
    # Build cumulative FZ loss matrix (row-constant, equal to column sums)
    FZmat <- matrix(colSums(FZmat_mcs),
                    nrow = nrow(FZmat_mcs),
                    ncol = ncol(FZmat_mcs),
                    byrow = TRUE)
    
    # Estimate the RSC lambda parameter
    temp <- lambda_esti(rets, tau, VaRm, ESm, FZmat)
    WL_RSC_lambda_store[[tt]] <- lambda <- temp$esti
    WL_RSC_AL_store[[tt]]     <- AL     <- temp$AL
    
    # Compute total FZ loss over the training period for each model
    FZmat_last <- matrix(rep(NA, nmod), ncol = nmod)
    for (j in 1:nmod) {
      FZmat_last[, j] <- sum(FZLoss(rets, VaRm[, j], ESm[, j], tau))
    }
    
    # Derive model weights and compute combined VaR and ES forecasts
    varf  <- VaR_oos[tt, MCS_models]
    esf   <- ES_oos[tt,  MCS_models]
    w     <- rel_comb(FZmat_last, lambda)  # RSC model weights
    cvarf <- sum(w * varf)
    cesf  <- sum(w * esf)
    
    WL_MCS_RSC_Comb_VaR[tt] <- cvarf
    WL_MCS_RSC_Comb_ES[tt]  <- cesf
  }
  
  
  # ----------------------------------------------------------------------------------
  # Method 5: MCS-MS-Comb — Minimum Score Combination (MSC) restricted to MCS-selected models
  #           (unweighted FZ loss; MSC method as in Taylor, 2020)
  # ----------------------------------------------------------------------------------
  
  MCS_models   <- MCS_est$includedR
  n_mcs_models <- length(MCS_models)
  
  if (n_mcs_models == 0) {
    # No models survive MCS: assign NA
    MCS_MSC_Comb_VaR[tt] <- NA
    MCS_MSC_Comb_ES[tt]  <- NA
    MSC_AL_store[tt]     <- NA
    
  } else if (n_mcs_models == 1) {
    # Only one model survives MCS: use it directly
    MCS_MSC_Comb_VaR[tt] <- VaR_oos[tt, MCS_models]
    MCS_MSC_Comb_ES[tt]  <- ES_oos[tt,  MCS_models]
    MSC_AL_store[tt]     <- NA
    
  } else {
    # Extract training-period forecasts for MCS-selected models
    VaRm_mcs <- as.matrix(VaR_training_data_mod[, MCS_models, tt])
    ESm_mcs  <- as.matrix(ES_training_data_mod[,  MCS_models, tt])
    rets     <- r_t_in_s_matrix[, tt]
    nmod     <- n_mcs_models
    
    # Estimate MSC beta coefficients
    fit  <- beta_esti_ms(rets, tau, VaRm_mcs, ESm_mcs)
    beta <- fit$esti
    AL   <- fit$AL
    MSC_AL_store[tt] <- AL
    
    # Transform beta to obtain quantile (betaq) and ES (betae) weights
    betamat <- transbeta(beta)
    betaq   <- betamat[, 1]
    betae   <- betamat[, 2]
    
    # Compute combined VaR and ES forecasts using MSC weights
    varf_mcs <- as.numeric(VaR_oos[tt, MCS_models])
    esf_mcs  <- as.numeric(ES_oos[tt,  MCS_models])
    
    cvarf <- sum(betaq * varf_mcs)
    cesf  <- cvarf + sum(betae * (esf_mcs - varf_mcs))
    
    MCS_MSC_Comb_VaR[tt] <- cvarf
    MCS_MSC_Comb_ES[tt]  <- cesf
  }
  
  
  # ----------------------------------------------------------------------------------
  # Method 6: WL-MCS-MS-Comb — Minimum Score Combination (MSC) restricted to MCS-selected models
  #           (EWMA-weighted FZ loss; MSC method as in Taylor, 2020)
  # ----------------------------------------------------------------------------------
  
  MCS_models   <- MCS_est_w$includedR
  n_mcs_models <- length(MCS_models)
  
  if (n_mcs_models == 0) {
    # No models survive MCS: assign NA
    WL_MCS_MSC_Comb_VaR[tt] <- NA
    WL_MCS_MSC_Comb_ES[tt]  <- NA
    WL_MSC_AL_store[tt]     <- NA
    
  } else if (n_mcs_models == 1) {
    # Only one model survives MCS: use it directly
    WL_MCS_MSC_Comb_VaR[tt] <- VaR_oos[tt, MCS_models]
    WL_MCS_MSC_Comb_ES[tt]  <- ES_oos[tt,  MCS_models]
    WL_MSC_AL_store[tt]     <- NA
    
  } else {
    # Extract training-period forecasts for MCS-selected models
    VaRm_mcs <- as.matrix(VaR_training_data_mod[, MCS_models, tt])
    ESm_mcs  <- as.matrix(ES_training_data_mod[,  MCS_models, tt])
    rets     <- r_t_in_s_matrix[, tt]
    nmod     <- n_mcs_models
    
    # Estimate MSC beta coefficients
    fit  <- beta_esti_ms(rets, tau, VaRm_mcs, ESm_mcs)
    beta <- fit$esti
    AL   <- fit$AL
    WL_MSC_AL_store[tt] <- AL
    
    # Transform beta to obtain quantile (betaq) and ES (betae) weights
    betamat <- transbeta(beta)
    betaq   <- betamat[, 1]
    betae   <- betamat[, 2]
    
    # Compute combined VaR and ES forecasts using MSC weights
    varf  <- VaR_oos[tt, MCS_models]
    esf   <- ES_oos[tt,  MCS_models]
    cvarf <- sum(betaq * varf)
    cesf  <- cvarf + sum(betae * (esf - varf))
    
    WL_MCS_MSC_Comb_VaR[tt] <- cvarf
    WL_MCS_MSC_Comb_ES[tt]  <- cesf
  }
  
  
   print(paste('estimation step',tt,'out of',nstep))
  
  
  #####################################################################################
  # 7. Save all results to file at each step (allows recovery if the loop is interrupted)
  #####################################################################################
  
  save(
    Backtesting_pvalues,
    Mean_Comb_VaR,   Mean_Comb_ES,
    Median_Comb_VaR, Median_Comb_ES,
    MCS_Comb_VaR,    MCS_Comb_ES,
    WL_MCS_Comb_VaR, WL_MCS_Comb_ES,
    MCS_RSC_Comb_VaR,    MCS_RSC_Comb_ES,
    WL_MCS_RSC_Comb_VaR, WL_MCS_RSC_Comb_ES,
    MCS_MSC_Comb_VaR,    MCS_MSC_Comb_ES,
    WL_MCS_MSC_Comb_VaR, WL_MCS_MSC_Comb_ES,
    MCS_est_included, MCS_est_included_w,
    tau,
    nstep,
    theta,
    alpha,
    lab_full,
    file = filename)
  
}



#####################################################################################
# 8. Relative Score Combining (Taylor, 2020) - Benchmark 3
#####################################################################################

source('functions/rsc_main.R')

RS_Comb_VaR<-cvar_rc_store
RS_Comb_ES<-ces_rc_store


#####################################################################################
# 9. Minimum Score Combining (Taylor, 2020) - Benchmark 4
#####################################################################################

source('functions/msc_main.R')

MS_Comb_VaR<-cvar_mc_store
MS_Comb_ES<-ces_mc_store


#####################################################################################
# 10. Store all the VaR and ES into VaR_oos_ev and ES_oos_ev objects
#####################################################################################

VaR_oos_ev<-cbind(
VaR_oos,
Mean_Comb_VaR, 		# Equally weighted combination
Median_Comb_VaR,		# Median combination
RS_Comb_VaR,			# RSC
MS_Comb_VaR,			# MSC
MCS_Comb_VaR,			# MCS based (loss unweighted)	
WL_MCS_Comb_VaR,		# MCS based (loss weighted)
MCS_RSC_Comb_VaR, 
WL_MCS_RSC_Comb_VaR, 
MCS_MSC_Comb_VaR, 
WL_MCS_MSC_Comb_VaR)

dim(VaR_oos_ev)

ES_oos_ev<-cbind(
ES_oos,
Mean_Comb_ES, 		# Equally weighted combination
Median_Comb_ES,		# Median combination
RS_Comb_ES,			# RSC
MS_Comb_ES,			# MSC
MCS_Comb_ES,			# MCS based (loss unweighted)	
WL_MCS_Comb_ES,		# MCS based (loss weighted)
MCS_RSC_Comb_ES, 
WL_MCS_RSC_Comb_ES, 
MCS_MSC_Comb_ES, 
WL_MCS_MSC_Comb_ES)

dim(ES_oos_ev)

################ full labels

lab_full<-c(list_of_models,
"EW-Comb","Median-Comb",
"RS-Comb","MS-Comb",
"MCS-Comb","WL-MCS-Comb",
"MCS-RS-Comb","WL-MCS-RS-Comb",
"MCS-MS-Comb","WL-MCS-MS-Comb")

colnames(VaR_oos_ev)<-colnames(ES_oos_ev)<-lab_full

############### Compute the FZLoss for the evaluation MCS

db_loss_oos<-matrix(NA,nrow=nstep,ncol=ncol(VaR_oos_ev))

for(i in 1:ncol(db_loss_oos)){
db_loss_oos[,i]<-FZLoss(coredata(r_t_oos_full),
coredata(VaR_oos_ev[,i]),coredata(ES_oos_ev[,i]),tau)
}

#####################################################################################
# 11. Save all the results
#####################################################################################

r_t_oos_full_plot<-r_t_oos_full
r_t_oos_full<- r_t_oos_full_plot[501:nstep]		# remove the first 500 observations as in the paper
db_loss_oos<-db_loss_oos[501:nstep,]		 	# remove the first 500 observations as in the paper
VaR_oos_ev<-VaR_oos_ev[501:nstep,]				# remove the first 500 observations as in the paper
ES_oos_ev<-ES_oos_ev[501:nstep,]				# remove the first 500 observations as in the paper

save(r_t_oos_full_plot,r_t_oos_full,nstep,
lab_full,tau,
Backtesting_pvalues,
db_loss_oos, VaR_oos_ev, ES_oos_ev, 
MCS_est_included,MCS_est_included_w,
file="data/results/sp500_final_results_generated.RData") # or "shanghai_comp_final_results_generated.RData"

