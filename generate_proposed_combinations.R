############################################################################################################################
### This code generates the six forecast combinations proposed in the paper:                                             ###
### Combining Value-at-Risk and Expected Shortfall Forecasts via the Model Confidence Set                                ###
### For ease of implementation and faster execution, we rely on a simplified (reduced) version of the original exercise. ###
############################################################################################################################

#### Analysis performed using R version 4.4.3

#####################################################################################
# 1. Load libraries and functions
#####################################################################################

### Load required libraries

library(np)
library(esback)
library(Rsolnp)
library(GAS)
library(xts)

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

# Load the example dataset containing all inputs required for the six forecast combinations
load("data/shanghai_comp_tau_0.025_comb_example.RData") 

#####################################################################################
#### Variables stored in "sp500_tau_0.025_comb_example.RData":
##   VaR_training_data_mod : array of VaR forecasts over the training period  - Tin x nmod x nstep
##   ES_training_data_mod  : array of ES forecasts over the training period   - Tin x nmod x nstep
##   r_t_in_s_matrix       : matrix of returns over the training period       - Tin x nstep
##   VaR_oos               : matrix of out-of-sample VaR forecasts            - nstep x nmod
##   ES_oos                : matrix of out-of-sample ES forecasts             - nstep x nmod
##   r_t_oos_full          : matrix of out-of-sample returns (xts/zoo format) - nstep x 1
##   r_t_oos               : matrix of out-of-sample returns                  - nstep x 1
##   list_of_models        : character vector of labels for the nmod competing models
##   N_model               : number of competing models (nmod)
#####################################################################################

# Output file name for saving results
filename <- "combinations_shanghai_comp_tau_0.025_theta_0.06_MCS_0.25_example.RData"


#####################################################################################
# 3. Set time span and risk level parameters
#####################################################################################

nstep <- 10    # Number of forecasting steps (model re-estimations) and out-of-sample period length
Tin   <- 1000  # Length of the training (in-sample) period
tau   <- 0.025 # Nominal coverage level for VaR and ES


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

MCS_Comb_VaR    <- MCS_Comb_ES    <-
  WL_MCS_Comb_VaR <- WL_MCS_Comb_ES <-
  MCS_RSC_Comb_VaR    <- MCS_RSC_Comb_ES    <-
  WL_MCS_RSC_Comb_VaR <- WL_MCS_RSC_Comb_ES <-
  MCS_MSC_Comb_VaR    <- MCS_MSC_Comb_ES    <-
  WL_MCS_MSC_Comb_VaR <- WL_MCS_MSC_Comb_ES <-
  rep(NA, nstep)

### Extend model labels to include the six combination method names
lab_full <- c(list_of_models,
              "MCS-Comb",
              "WL-MCS-Comb",
              "MCS-RS-Comb",
              "WL-MCS-RS-Comb",
              "MCS-MS-Comb",
              "WL-MCS-MS-Comb")

### Initialize lists for storing intermediate MCS results and combination diagnostics
MCS_est_included   <- MCS_est_included_w <- list()

RSC_lambda_store   <- RSC_AL_store   <-
  WL_RSC_lambda_store <- WL_RSC_AL_store <-
  MSC_AL_store       <- WL_MSC_AL_store <- list()


#####################################################################################
# 6. Main loop: compute forecast combinations at each step using MCS-based methods
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
  
  
  message(tt)  # Print progress: current forecasting step
  
  
  #####################################################################################
  # 7. Save all results to file at each step (allows recovery if the loop is interrupted)
  #####################################################################################
  
  save(
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





