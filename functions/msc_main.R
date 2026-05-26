# =============================================================================
# MSC FORECASTING - MAIN CODE
# =============================================================================
# Implements the Minimum Score Combining (MSC) method for VaR and ES forecast
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

wstore        <- matrix(nrow = nstep, ncol = nmod)   # combination weights [nstep x nmod]
cvar_mc_store <- matrix(nrow = nstep, ncol = 1)      # combined VaR forecasts
ces_mc_store  <- matrix(nrow = nstep, ncol = 1)      # combined ES forecasts
betaq_store   <- matrix(0, nrow = nstep, ncol = nmod) # VaR combination coefficients (beta_q)
betae_store   <- matrix(0, nrow = nstep, ncol = nmod) # ES combination coefficients (beta_e)
AL_store      <- matrix(0, nrow = nstep, ncol = 1)   # minimized AL score at each step
check         <- rep(1, nstep)   # convergence flag: 0 = success, 1 = fallback


# --- Warm-start values -------------------------------------------------------
# Equal weights across models; carried forward if the optimizer fails

betaq_old <- c(0, rep(1 / nmod, nmod))
betae_old <- c(0, rep(1 / nmod, nmod))


# =============================================================================
# FORECASTING LOOP
# =============================================================================

for (i in 1:nstep) {

  print(paste('estimation step', i, 'out of', nstep))

  # Training data for step i
  VaRm <- VaR_training_data_mod[, , i]   # in-sample VaR forecasts [Tin x nmod]
  ESm  <- ES_training_data_mod[, , i]    # in-sample ES forecasts  [Tin x nmod]
  rets <- r_t_in_s_matrix[, i]           # realized returns        [Tin x 1]


  # Estimate MSC weights by minimizing the AL (joint VaR-ES) loss.
  # tryCatch returns NA if the optimizer encounters a numerical error.
  temp <- tryCatch(
    { fit <- beta_esti_ms(rets, alfa, VaRm, ESm) },
    error = function(msg) { return(NA) }
  )


  if (class(temp) != "list") {
    # Optimizer failed: reuse last valid coefficients
    betae <- betae_old
    betaq <- betaq_old

  } else {
    beta <- fit$esti   # raw parameter vector from the optimizer
    AL   <- fit$AL     # minimized AL score

    # transbeta() maps beta to constrained (betaq, betae) weight matrices
    betamat <- transbeta(beta)
    betaq   <- betamat[, 1]   # VaR weights
    betae   <- betamat[, 2]   # ES weights

    check[i]    <- 0
    AL_store[i] <- AL
  }

  # Update warm-start and store coefficients
  betaq_old        <- betaq
  betae_old        <- betae
  betaq_store[i, ] <- betaq
  betae_store[i, ] <- betae


  # Out-of-sample forecasts
  varf <- VaR_oos[i, ]   # individual VaR forecasts [1 x nmod]
  esf  <- ES_oos[i, ]    # individual ES forecasts  [1 x nmod]

  # Combined VaR: weighted sum of individual VaR forecasts
  cvarf <- sum(betaq * varf)

  # Combined ES: VaR combination + tail correction term
  # ES_combined = VaR_combined + sum(betae * (ES_i - VaR_i))
  cesf <- cvarf + sum(betae * (esf - varf))

  cvar_mc_store[i] <- cvarf   # combined VaR
  ces_mc_store[i]  <- cesf    # combined ES

} # end forecasting loop