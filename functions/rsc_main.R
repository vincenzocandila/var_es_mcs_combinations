# =============================================================================
# RSC FORECASTING — MAIN CODE
# =============================================================================
# Implements the Relative Score Combining (RSC) method for VaR and ES forecast
# combination as proposed in:
#   Taylor, J.W. (2020). Forecast combinations for Value at Risk and
#   Expected Shortfall. International Journal of Forecasting, 36(2), 428-441.
# =============================================================================


# -----------------------------------------------------------------------------
# INITIALIZATION
# -----------------------------------------------------------------------------

alfa <- tau          # quantile level for VaR/ES (e.g., 0.01 or 0.025)
Tin  <- dim(VaR_training_data_mod)[1]   # number of in-sample observations
nmod <- N_model                         # number of individual models


# --- Storage matrices --------------------------------------------------------

FZmat            <- matrix(0, nrow = 1000, ncol = N_model)   # cumulative FZ scores over training window [Tin x nmod]
FZmat_for        <- matrix(0, nrow = nstep, ncol = N_model)  # FZ scores at each out-of-sample step
wstore           <- matrix(nrow = nstep, ncol = N_model)     # RSC combination weights [nstep x nmod]
cvar_rc_store    <- matrix(nrow = nstep, ncol = 1)           # combined VaR forecasts
ces_rc_store     <- matrix(nrow = nstep, ncol = 1)           # combined ES forecasts
ms_cvar_rc_store <- matrix(nrow = nstep, ncol = 1)           # model-selected VaR forecasts (best single model)
ms_ces_rc_store  <- matrix(nrow = nstep, ncol = 1)           # model-selected ES forecasts (best single model)
lambda_store     <- matrix(0, nrow = nstep, ncol = 1)        # estimated temperature parameter lambda
AL_store         <- matrix(0, nrow = nstep, ncol = 1)        # minimized AL score at each step


# =============================================================================
# FORECASTING LOOP
# =============================================================================

for (i in 1:nstep) {

  print(paste('estimation step', i, 'out of', nstep))

  # Training data for step i
  VaRm <- VaR_training_data_mod[, , i]   # in-sample VaR forecasts [Tin x nmod]
  ESm  <- ES_training_data_mod[, , i]    # in-sample ES forecasts  [Tin x nmod]
  rets <- r_t_in_s_matrix[, i]           # realized returns        [Tin x 1]


  # Compute cumulative FZ (Fissler-Ziegel) loss score for each model
  FZmat_last <- matrix(rep(NA, nmod), ncol = nmod)   # [1 x nmod]

  for (j in 1:nmod) {
    FZmat_last[, j] <- sum(FZLoss(rets, VaRm[, j], ESm[, j], alfa))
  }

  # Replicate the score row Tin times for use in lambda_esti()
  FZmat <- matrix(rep(FZmat_last, Tin), byrow = TRUE, ncol = nmod)


  # Estimate lambda: temperature parameter controlling weight concentration.
  # lambda -> 0 gives equal weights; lambda -> inf selects the best model only.
  temp   <- lambda_esti(rets, alfa, VaRm, ESm, FZmat)

  lambda <- temp$esti   # estimated temperature parameter
  AL     <- temp$AL     # minimized AL score

  lambda_store[i] <- lambda
  AL_store[i]     <- AL


  # Out-of-sample forecasts
  varf <- VaR_oos[i, ]   # individual VaR forecasts [1 x nmod]
  esf  <- ES_oos[i, ]    # individual ES forecasts  [1 x nmod]

  # rel_comb() applies a softmax over FZ scores scaled by lambda -> weight vector [1 x nmod]
  w <- rel_comb(FZmat_last, lambda)

  # Combined VaR and ES: weighted sums using RSC weights
  cvarf <- sum(w * varf)
  cesf  <- sum(w * esf)

  cvar_rc_store[i] <- cvarf   # combined VaR
  ces_rc_store[i]  <- cesf    # combined ES
  wstore[i, ]      <- w       # RSC weights

} # end forecasting loop