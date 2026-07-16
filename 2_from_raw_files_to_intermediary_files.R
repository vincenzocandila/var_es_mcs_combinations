############################################################################################################################
### This code generates the VaR and ES forecasts of the 32 models used in the paper:                                     ###
### 'Combining Value-at-Risk and Expected Shortfall Forecasts via the Model Confidence Set'                              ###
############################################################################################################################

#### Analysis performed using R version 4.4.3

###################################
#### Load the libraries
###################################

library(xts)				# Version used: 0.14.1
library(zoo)				# Version used: 1.8.13
library(fGarch)				# Version used: 4033.92
library(rumidas)			# Version used: 0.1.3
library(rugarch)			# Version used: 1.5.3
library(quantreg)			# Version used: 6.1

###################################
#### Load the functions
###################################

source("functions/main_functions.R")

###################################
#### Set the period of interest
###################################

period_of_interest<-"2013-01/2022-06"

###################################
#### Set seed
###################################

set.seed(123)

#####################################################################
#### Load S&P 500 or Shanghai Composite or
#### one of the other seven indices used in the paper:
#### BOVESPA, BSESN, EUROSTOXX50, HSI, IXIC, MXX, and NIKKEI
#####################################################################

load("data/raw/sp500_data.RData")  

r_t_est	<-	r_t[period_of_interest]
rvol5	<-	rvol5[period_of_interest]
rb_ss	<-	rb_ss[period_of_interest]
rk		<-	rk[period_of_interest]

TT<-length(r_t_est)


#####################################################################
#### sp500_data (as well as the datasets for the other indices) includes:
#### r_t    : daily close-to-close log returns;
#### rb_ss  : realized volatility based on realized bipower variation
####          with subsampling;
#### rk     : realized kernel volatility;
#### rvol5  : realized volatility computed from 5-minute returns.
#####################################################################

###### Economic Policy Uncertainty data collected from:
###### https://www.policyuncertainty.com/global_monthly.html
###### We use the column labeled "GEPU_current"
###### The EPU indices were downloaded on May 26, 2023. 
###### Note that the series may be subject to revisions over time

epu_full<-read.csv(file="data/raw/Global_Policy_Uncertainty_Data.csv",sep=";",dec=",")
epu_full_i <- make_xts(epu_full, value_col = names(epu_full)[3])

MV<-diff(epu_full_i)/lag(epu_full_i)
MV[1]<-0

##################################################
################################################## VaR, ES, and MIDAS settings
##################################################

tau	<-	0.025	# Coverage level

K	<-	36 		# Number of lagged realizations entering the long-run equation

##################################################
################################################## OOS settings
##################################################

Tin		<- 1000 	# length of rolling period
lstep	<- 1 	# frequency of re-estimation

nstep<-(TT-Tin)

nstep


################################################## Historical Simulation

ST<-"2005"
EN<-last(index(r_t_est))

################################################## HS-25

w<-25 # length of rolling window

HS_25<-hs_f(r_t,w,tau,ST,EN)

VaR_HS_25<-HS_25$VaR_HS
ES_HS_25<-HS_25$ES_HS


################################################## HS-50

w<-50 # length of rolling window

HS_50<-hs_f(r_t,w,tau,ST,EN)

VaR_HS_50<-HS_50$VaR_HS

ES_HS_50<-HS_50$ES_HS

################################################## HS-100

w<-100 # length of rolling window

HS_100<-hs_f(r_t,w,tau,ST,EN)

VaR_HS_100<-HS_100$VaR_HS
ES_HS_100<-HS_100$ES_HS

################################################## HS-250

w<-250 # length of rolling window

HS_250<-hs_f(r_t,w,tau,ST,EN)

VaR_HS_250<-HS_250$VaR_HS
ES_HS_250<-HS_250$ES_HS

################################################## HS-500

w<-500 # length of rolling window

HS_500<-hs_f(r_t,w,tau,ST,EN)

VaR_HS_500<-HS_500$VaR_HS
ES_HS_500<-HS_500$ES_HS

##################################################
################################################## list and array settings
################################################## 

list_of_models<-lab<-c(
"RM-N","RM-N-CF","RM-t",
"GARCH-N","GARCH-N-CF","GARCH-t",
"GJR-N","GJR-N-CF","GJR-t",
"RGARCH-RVOL5-N","RGARCH-RVOL5-N-CF","RGARCH-RVOL5-t",
"RGARCH-RB-SS-N","RGARCH-RB-SS-N-CF","RGARCH-RB-SS-t",
"RGARCH-RK-N","RGARCH-RK-N-CF","RGARCH-RK-t",
"HS-25","HS-50","HS-100","HS-250","HS-500",
"SAV","AS","IG",
"CAViaR-X-RVOL5","CAViaR-X-RB-SS","CAViaR-X-RK",
"MF-X-RVOL5","MF-X-RB-SS","MF-X-RK")

M<-N_model<-length(list_of_models)

db_VaR_combined_full<-db_ES_combined_full<-db_VaR_combined<-db_ES_combined<-MCS_est_included<-list()
db_VaR_combined_w<-db_ES_combined_w<-MCS_est_included_w<-list()

db_loss_in_s_array<-VaR_in_s_array<-ES_in_s_array<-array(NA,dim=c(Tin,M,nstep))
dim(db_loss_in_s_array)

db_loss_oos_array<-VaR_oos_array<-ES_oos_array<-array(NA,dim=c(lstep,M,nstep))
dim(db_loss_in_s_array)

VaR_training_data_mod<-ES_training_data_mod<-array(NA,dim=c(Tin,M,nstep))

VaR_oos<-ES_oos<-matrix(NA,ncol=M,nrow=nstep)
colnames(VaR_oos)<-colnames(ES_oos)<-list_of_models

r_t_in_s_matrix<-matrix(NA,nrow=Tin,ncol=nstep)
dim(r_t_in_s_matrix)

coef_sav_previous <- NULL
VaR_SAV_in_s_previous <- NULL
ES_SAV_in_s_previous <- NULL

failed_SAV <- integer(0)

coef_as_previous <- NULL
VaR_AS_in_s_previous <- NULL
ES_AS_in_s_previous <- NULL

failed_AS <- integer(0)

coef_ig_previous <- NULL
VaR_IG_in_s_previous <- NULL
ES_IG_in_s_previous <- NULL

failed_IG <- integer(0)

coef_sav_rvol_5_previous <- NULL
VaR_SAVX_rvol_5_in_s_previous <- NULL
ES_SAVX_rvol_5_in_s_previous <- NULL

failed_SAVX_rvol_5 <- integer(0)

coef_sav_rb_ss_previous <- NULL
VaR_SAVX_rb_ss_in_s_previous <- NULL
ES_SAVX_rb_ss_in_s_previous <- NULL

failed_SAVX_rb_ss <- integer(0)

coef_sav_rk_previous <- NULL
VaR_SAVX_rk_in_s_previous <- NULL
ES_SAVX_rk_in_s_previous <- NULL

failed_SAVX_rk<- integer(0)

##################################################
################################################## MODEL specs
################################################## 

spec_garch_n <- ugarchspec(variance.model=list(model="sGARCH", garchOrder=c(1,1)), 
		mean.model=list(armaOrder=c(0,0), include.mean=FALSE),  
		distribution.model="norm")

spec_garch_t <- ugarchspec(variance.model=list(model="sGARCH", garchOrder=c(1,1)), 
		mean.model=list(armaOrder=c(0,0), include.mean=FALSE),  
		distribution.model="std")

spec_gjr_n <- ugarchspec(variance.model=list(model="gjrGARCH", garchOrder=c(1,1)), 
		mean.model=list(armaOrder=c(0,0), include.mean=FALSE),  
		distribution.model="norm")

spec_gjr_t <- ugarchspec(variance.model=list(model="gjrGARCH", garchOrder=c(1,1)), 
		mean.model=list(armaOrder=c(0,0), include.mean=FALSE),  
		distribution.model="std")

spec_re_garch_n <- ugarchspec(variance.model=list(model="realGARCH", garchOrder=c(1,1)), 
		mean.model=list(armaOrder=c(0,0), include.mean=FALSE),
		distribution.model="norm")

spec_re_garch_t <- ugarchspec(variance.model=list(model="realGARCH", garchOrder=c(1,1)), 
		mean.model=list(armaOrder=c(0,0), include.mean=FALSE),  
		distribution.model="std")

lambda_rm<-0.94 # RiskMetrics decay factor

df_l<-list()

z_a<-qnorm(tau)

# Numerical stability thresholds used to discard unreliable model estimates

# Maximum acceptable Nyblom stability statistic
NYBLOM_MAX <- 10

# Maximum acceptable Hessian condition number
KAPPA_MAX <- 1e6

# Acceptable range for the Student-t degrees-of-freedom parameter
SHAPE_MIN <- 2.05
SHAPE_MAX <- 100

# Number of random restarts used by the gosolnp optimizer
GOSOLNP_RESTARTS <- 15



filename<-"data/intermediary/sp500_tau_0.025_gen_files.RData" 

##################################################
##################################################
################################################## FULL Cycle
################################################## 
##################################################


for(tt in 1:nstep){  

day_end_est<-tt+Tin

r_t_est_cycle<-r_t_est[tt:day_end_est]

r_t_in_s<-r_t_est_cycle
r_t_in_s_tin<-coredata(r_t_est_cycle)[-length(r_t_est_cycle)]
r_t_in_s_tin_xts<-(r_t_est_cycle)[-length(r_t_est_cycle)]

r_t_oos<-r_t_est[tt+Tin]

r_t_oos_full<-r_t_est[(Tin+1):length(r_t_est)]

r_t_in_s_matrix[,tt]<-coredata(r_t_in_s_tin)

################################################## X variables

rvol_5_est_cycle<-rvol5[tt:day_end_est]
rvol_5_est_cycle<-fix_zeros(rvol_5_est_cycle)

rb_ss_est_cycle<-rb_ss[tt:day_end_est]
rb_ss_est_cycle<-fix_zeros(rb_ss_est_cycle)

rk_est_cycle<-rk[tt:day_end_est]
rk_est_cycle<-fix_zeros(rk_est_cycle)

################################################## MIDAS variables

epu_mv_cycle<-mv_into_mat(r_t_est_cycle,MV,K=K,"monthly")

##################################################
################################################## Parametric methods
################################################## (from M1 to M18)

############################################ 
############################################ RISKMETRICS (NORMAL)
############################################ M1

m<-1

sigma_rm<-(riskmetrics_f(r_t_est_cycle,lambda_rm))^0.5
sigma_rm_in_s<-sigma_rm[1:Tin]
sigma_rm_oos<-sigma_rm[(Tin+1)]

VaR_rm_n_in_s<-z_a*sigma_rm_in_s
ES_rm_n_in_s<- -sigma_rm_in_s*dnorm(z_a)/tau

VaR_oos[tt,m]<-qnorm(tau)*sigma_rm_oos
ES_oos[tt,m]<- -sigma_rm_oos*dnorm(qnorm(tau))/tau

if (tt == 1) {
      # first step, only in-sample
      VaR_training_data_mod[, m, tt] <- VaR_rm_n_in_s
      ES_training_data_mod[, m, tt] <- ES_rm_n_in_s
      
    } else if (tt <= Tin) {
      # Mixed window: first in-sample, last oos
      n_in <- Tin - tt + 1
      n_oos <- (n_in+1):Tin
      
      VaR_training_data_mod[1:n_in, m, tt] <- VaR_rm_n_in_s[1:n_in]
      VaR_training_data_mod[n_oos, m, tt] <- VaR_oos[1:(tt-1), m]
      
      ES_training_data_mod[1:n_in, m, tt] <- ES_rm_n_in_s[1:n_in]
      ES_training_data_mod[n_oos, m, tt] <- ES_oos[1:(tt-1), m]
      
    } else {
      # Only oos 
      VaR_training_data_mod[, m, tt] <- VaR_oos[(tt - Tin):(tt-1), m]
      ES_training_data_mod[, m, tt] <- ES_oos[(tt - Tin):(tt-1), m]
    }
  


############################################ 
############################################ RISKMETRICS (NORMAL + cf)
############################################ M2

m<-2

z_cf_VaR<-Quantile_CF_f(r_t_in_s_tin,sigma_rm_in_s,tau)[,1]
z_cf_ES<-Quantile_CF_f(r_t_in_s_tin,sigma_rm_in_s,tau)[,2]

VaR_rm_n_cf_in_s<- sigma_rm_in_s*z_cf_VaR
ES_rm_n_cf_in_s<- sigma_rm_in_s*z_cf_ES

VaR_oos[tt,m]<- sigma_rm_oos*z_cf_VaR
ES_oos[tt,m]<- sigma_rm_oos*z_cf_ES

if (tt == 1) {
      # first step, only in-sample
      VaR_training_data_mod[, m, tt] <- VaR_rm_n_cf_in_s
      ES_training_data_mod[, m, tt] <- ES_rm_n_cf_in_s
      
    } else if (tt <= Tin) {
      # Mixed window: first in-sample, last oos
      n_in <- Tin - tt + 1
      n_oos <- (n_in+1):Tin
      
      VaR_training_data_mod[1:n_in, m, tt] <- VaR_rm_n_cf_in_s[1:n_in]
      VaR_training_data_mod[n_oos, m, tt] <- VaR_oos[1:(tt-1), m]
      
      ES_training_data_mod[1:n_in, m, tt] <- ES_rm_n_cf_in_s[1:n_in]
      ES_training_data_mod[n_oos, m, tt] <- ES_oos[1:(tt-1), m]
      
    } else {
      # Only oos 
      VaR_training_data_mod[, m, tt] <- VaR_oos[(tt - Tin):(tt-1), m]
      ES_training_data_mod[, m, tt] <- ES_oos[(tt - Tin):(tt-1), m]
    }
  
############################################ 
############################################ RISKMETRICS (Student-t)
############################################ M3

m<-3

df_l[[tt]]<-df_hat<-as.numeric(stdFit(r_t_in_s_tin/sigma_rm_in_s)$par[3])

if(df_l[[tt]]<2.01){

df_hat<-df_l[[tt-1]] 

} 

correc<-((df_hat-2)/df_hat)^0.5

VaR_oos[tt,m]<-qt(tau,df_hat)*sigma_rm_oos*correc

## first term:
first_term<-(dt(qt(tau,df_hat),df_hat)/tau)

## second term:
second_term<- (df_hat+(qt(tau,df_hat))^2)/(df_hat-1)

ES_oos[tt,m]<- -sigma_rm_oos*first_term*second_term*correc

VaR_rm_t_in_s<-qt(tau,df_hat)*sigma_rm_in_s*correc
ES_rm_t_in_s<- -sigma_rm_in_s*first_term*second_term*correc


if (tt == 1) {
      # first step, only in-sample
      VaR_training_data_mod[, m, tt] <- VaR_rm_t_in_s
      ES_training_data_mod[, m, tt] <- ES_rm_t_in_s
      
    } else if (tt <= Tin) {
      # Mixed window: first in-sample, last oos
      n_in <- Tin - tt + 1
      n_oos <- (n_in+1):Tin
      
      VaR_training_data_mod[1:n_in, m, tt] <- VaR_rm_t_in_s[1:n_in]
      VaR_training_data_mod[n_oos, m, tt] <- VaR_oos[1:(tt-1), m]
      
      ES_training_data_mod[1:n_in, m, tt] <- ES_rm_t_in_s[1:n_in]
      ES_training_data_mod[n_oos, m, tt] <- ES_oos[1:(tt-1), m]
      
    } else {
      # Only oos 
      VaR_training_data_mod[, m, tt] <- VaR_oos[(tt - Tin):(tt-1), m]
      ES_training_data_mod[, m, tt] <- ES_oos[(tt - Tin):(tt-1), m]
    }
  



############################################ 
############################################ GARCH(1,1) with a normal distribution
############################################ M4

m<-4

fit_garch_n<-ugarchfit(spec=spec_garch_n, data=r_t_est_cycle,out.sample=lstep,solver="hybrid")

sigma_garch_n_in_s<-fit_garch_n@fit$sigma
sigma_garch_n_oos<-as.numeric(sigma(ugarchforecast(fit_garch_n,  n.ahead = 1)))

VaR_garch_n_in_s<-z_a*sigma_garch_n_in_s
ES_garch_n_in_s<- -sigma_garch_n_in_s*dnorm(z_a)/tau

VaR_oos[tt,m]<- z_a*sigma_garch_n_oos
ES_oos[tt,m]<- -sigma_garch_n_oos*dnorm(z_a)/tau



if (tt == 1) {
      # first step, only in-sample
      VaR_training_data_mod[, m, tt] <- VaR_garch_n_in_s
      ES_training_data_mod[, m, tt] <- ES_garch_n_in_s
      
    } else if (tt <= Tin) {
      # Mixed window: first in-sample, last oos
      n_in <- Tin - tt + 1
      n_oos <- (n_in+1):Tin
      
      VaR_training_data_mod[1:n_in, m, tt] <- VaR_garch_n_in_s[1:n_in]
      VaR_training_data_mod[n_oos, m, tt] <- VaR_oos[1:(tt-1), m]
      
      ES_training_data_mod[1:n_in, m, tt] <- ES_garch_n_in_s[1:n_in]
      ES_training_data_mod[n_oos, m, tt] <- ES_oos[1:(tt-1), m]
      
    } else {
      # Only oos 
      VaR_training_data_mod[, m, tt] <- VaR_oos[(tt - Tin):(tt-1), m]
      ES_training_data_mod[, m, tt] <- ES_oos[(tt - Tin):(tt-1), m]
    }
  


############################################ 
############################################ GARCH(1,1) with a normal distribution and CF approx
############################################ M5

m<-5

z_cf_VaR<-Quantile_CF_f(r_t_in_s_tin,sigma_garch_n_in_s,tau)[,1]
z_cf_ES<-Quantile_CF_f(r_t_in_s_tin,sigma_garch_n_in_s,tau)[,2]

VaR_garch_n_cf_in_s<- sigma_garch_n_in_s*z_cf_VaR
ES_garch_n_cf_in_s<- sigma_garch_n_in_s*z_cf_ES

VaR_oos[tt,m]<- sigma_garch_n_oos*z_cf_VaR
ES_oos[tt,m]<- sigma_garch_n_oos*z_cf_ES


if (tt == 1) {
      # first step, only in-sample
      VaR_training_data_mod[, m, tt] <- VaR_garch_n_cf_in_s
      ES_training_data_mod[, m, tt] <- ES_garch_n_cf_in_s
      
    } else if (tt <= Tin) {
      # Mixed window: first in-sample, last oos
      n_in <- Tin - tt + 1
      n_oos <- (n_in+1):Tin
      
      VaR_training_data_mod[1:n_in, m, tt] <- VaR_garch_n_cf_in_s[1:n_in]
      VaR_training_data_mod[n_oos, m, tt] <- VaR_oos[1:(tt-1), m]
      
      ES_training_data_mod[1:n_in, m, tt] <- ES_garch_n_cf_in_s[1:n_in]
      ES_training_data_mod[n_oos, m, tt] <- ES_oos[1:(tt-1), m]
      
    } else {
      # Only oos 
      VaR_training_data_mod[, m, tt] <- VaR_oos[(tt - Tin):(tt-1), m]
      ES_training_data_mod[, m, tt] <- ES_oos[(tt - Tin):(tt-1), m]
    }
  


############################################ 
############################################ GARCH(1,1) with a student-t distribution
############################################ M6

m<-6

fit_garch_t<-ugarchfit(spec=spec_garch_t, data=r_t_est_cycle,out.sample=lstep,solver="hybrid")

sigma_garch_t_oos<-as.numeric(sigma(ugarchforecast(fit_garch_t,  n.ahead = 1)))

df_hat<-as.numeric(fit_garch_t@fit$coef['shape'])

correc<-((df_hat-2)/df_hat)^0.5

VaR_oos[tt,m]<-qt(tau,df_hat)*sigma_garch_t_oos*correc

## first term:
first_term<-(dt(qt(tau,df_hat),df_hat)/tau)

## second term:
second_term<- (df_hat+(qt(tau,df_hat))^2)/(df_hat-1)

ES_oos[tt,m]<- -sigma_garch_t_oos*first_term*second_term*correc

VaR_garch_t_in_s<-qt(tau,df_hat)*fit_garch_t@fit$sigma*correc
ES_garch_t_in_s<- -fit_garch_t@fit$sigma*first_term*second_term*correc

if (tt == 1) {
      # first step, only in-sample
      VaR_training_data_mod[, m, tt] <- VaR_garch_t_in_s
      ES_training_data_mod[, m, tt] <- ES_garch_t_in_s
      
    } else if (tt <= Tin) {
      # Mixed window: first in-sample, last oos
      n_in <- Tin - tt + 1
      n_oos <- (n_in+1):Tin
      
      VaR_training_data_mod[1:n_in, m, tt] <- VaR_garch_t_in_s[1:n_in]
      VaR_training_data_mod[n_oos, m, tt] <- VaR_oos[1:(tt-1), m]
      
      ES_training_data_mod[1:n_in, m, tt] <- ES_garch_t_in_s[1:n_in]
      ES_training_data_mod[n_oos, m, tt] <- ES_oos[1:(tt-1), m]
      
    } else {
      # Only oos 
      VaR_training_data_mod[, m, tt] <- VaR_oos[(tt - Tin):(tt-1), m]
      ES_training_data_mod[, m, tt] <- ES_oos[(tt - Tin):(tt-1), m]
    }
  



############################################ 
############################################ GJR(1,1) with a normal distribution
############################################ M7

m<-7

fit_gjr_n<-ugarchfit(spec=spec_gjr_n, data=r_t_est_cycle,out.sample=lstep,solver="hybrid")

sigma_gjr_n_in_s<-fit_gjr_n@fit$sigma
sigma_gjr_n_oos<-as.numeric(sigma(ugarchforecast(fit_gjr_n,  n.ahead = 1)))

VaR_oos[tt,m]<-z_a*sigma_gjr_n_oos
ES_oos[tt,m]<- -sigma_gjr_n_oos*dnorm(z_a)/tau

VaR_gjr_n_in_s<-z_a*fit_gjr_n@fit$sigma
ES_gjr_n_in_s<- -fit_gjr_n@fit$sigma*dnorm(z_a)/tau

if (tt == 1) {
      # first step, only in-sample
      VaR_training_data_mod[, m, tt] <- VaR_gjr_n_in_s
      ES_training_data_mod[, m, tt] <- ES_gjr_n_in_s
      
    } else if (tt <= Tin) {
      # Mixed window: first in-sample, last oos
      n_in <- Tin - tt + 1
      n_oos <- (n_in+1):Tin
      
      VaR_training_data_mod[1:n_in, m, tt] <- VaR_gjr_n_in_s[1:n_in]
      VaR_training_data_mod[n_oos, m, tt] <- VaR_oos[1:(tt-1), m]
      
      ES_training_data_mod[1:n_in, m, tt] <- ES_gjr_n_in_s[1:n_in]
      ES_training_data_mod[n_oos, m, tt] <- ES_oos[1:(tt-1), m]
      
    } else {
      # Only oos 
      VaR_training_data_mod[, m, tt] <- VaR_oos[(tt - Tin):(tt-1), m]
      ES_training_data_mod[, m, tt] <- ES_oos[(tt - Tin):(tt-1), m]
    }
  



############################################ 
############################################ GJR(1,1) with a normal distribution and CF approx
############################################ M8

m<-8

z_cf_VaR<-Quantile_CF_f(r_t_in_s_tin,sigma_gjr_n_in_s,tau)[,1]
z_cf_ES<-Quantile_CF_f(r_t_in_s_tin,sigma_gjr_n_in_s,tau)[,2]

VaR_gjr_n_cf_in_s<- sigma_gjr_n_in_s*z_cf_VaR
ES_gjr_n_cf_in_s<- sigma_gjr_n_in_s*z_cf_ES

VaR_oos[tt,m]<- sigma_gjr_n_oos*z_cf_VaR
ES_oos[tt,m]<- sigma_gjr_n_oos*z_cf_ES

if (tt == 1) {
      # first step, only in-sample
      VaR_training_data_mod[, m, tt] <- VaR_gjr_n_cf_in_s
      ES_training_data_mod[, m, tt] <- ES_gjr_n_cf_in_s
      
    } else if (tt <= Tin) {
      # Mixed window: first in-sample, last oos
      n_in <- Tin - tt + 1
      n_oos <- (n_in+1):Tin
      
      VaR_training_data_mod[1:n_in, m, tt] <- VaR_gjr_n_cf_in_s[1:n_in]
      VaR_training_data_mod[n_oos, m, tt] <- VaR_oos[1:(tt-1), m]
      
      ES_training_data_mod[1:n_in, m, tt] <- ES_gjr_n_cf_in_s[1:n_in]
      ES_training_data_mod[n_oos, m, tt] <- ES_oos[1:(tt-1), m]
      
    } else {
      # Only oos 
      VaR_training_data_mod[, m, tt] <- VaR_oos[(tt - Tin):(tt-1), m]
      ES_training_data_mod[, m, tt] <- ES_oos[(tt - Tin):(tt-1), m]
    }
  


############################################ 
############################################ GJR(1,1) with a student-t distribution
############################################ M9

m<-9

fit_gjr_t<-ugarchfit(spec=spec_gjr_t, data=r_t_est_cycle,out.sample=lstep,solver="hybrid")

sigma_gjr_t_oos<-as.numeric(sigma(ugarchforecast(fit_gjr_t,  n.ahead = 1)))

df_hat<-as.numeric(fit_gjr_t@fit$coef['shape'])

correc<-((df_hat-2)/df_hat)^0.5

VaR_oos[tt,m]<-qt(tau,df_hat)*sigma_gjr_t_oos*correc

## first term:
first_term<-(dt(qt(tau,df_hat),df_hat)/tau)

## second term:
second_term<- (df_hat+(qt(tau,df_hat))^2)/(df_hat-1)

ES_oos[tt,m]<- -sigma_gjr_t_oos*first_term*second_term*correc

VaR_gjr_t_in_s<-qt(tau,df_hat)*fit_gjr_t@fit$sigma*correc
ES_gjr_t_in_s<- -fit_gjr_t@fit$sigma*first_term*second_term*correc

if (tt == 1) {
      # first step, only in-sample
      VaR_training_data_mod[, m, tt] <- VaR_gjr_t_in_s
      ES_training_data_mod[, m, tt] <- ES_gjr_t_in_s
      
    } else if (tt <= Tin) {
      # Mixed window: first in-sample, last oos
      n_in <- Tin - tt + 1
      n_oos <- (n_in+1):Tin
      
      VaR_training_data_mod[1:n_in, m, tt] <- VaR_gjr_t_in_s[1:n_in]
      VaR_training_data_mod[n_oos, m, tt] <- VaR_oos[1:(tt-1), m]
      
      ES_training_data_mod[1:n_in, m, tt] <- ES_gjr_t_in_s[1:n_in]
      ES_training_data_mod[n_oos, m, tt] <- ES_oos[1:(tt-1), m]
      
    } else {
      # Only oos 
      VaR_training_data_mod[, m, tt] <- VaR_oos[(tt - Tin):(tt-1), m]
      ES_training_data_mod[, m, tt] <- ES_oos[(tt - Tin):(tt-1), m]
    }
  



############################################ 
############################################ RGARCH with a normal distribution (with rvol_5)
############################################ M10

m<-10

rg_est_rvol_5_n<-fit_robust(r_t_est_cycle, rvol_5_est_cycle, 
spec_re_garch_n, is_t = FALSE)$fit

sigma_rg_n_rvol_5_in_s<-rg_est_rvol_5_n@fit$sigma/100
sigma_rg_n_rvol_5_oos<-as.numeric(sigma(ugarchforecast(rg_est_rvol_5_n,  n.ahead = 1)))/100

VaR_oos[tt,m]<-z_a*sigma_rg_n_rvol_5_oos
ES_oos[tt,m]<- -sigma_rg_n_rvol_5_oos*dnorm(z_a)/tau

VaR_rg_n_rvol_5_in_s<-z_a*rg_est_rvol_5_n@fit$sigma
ES_rg_n_rvol_5_in_s<- -rg_est_rvol_5_n@fit$sigma*dnorm(z_a)/tau

if (tt == 1) {
      # first step, only in-sample
      VaR_training_data_mod[, m, tt] <- VaR_rg_n_rvol_5_in_s
      ES_training_data_mod[, m, tt] <- ES_rg_n_rvol_5_in_s
      
    } else if (tt <= Tin) {
      # Mixed window: first in-sample, last oos
      n_in <- Tin - tt + 1
      n_oos <- (n_in+1):Tin
      
      VaR_training_data_mod[1:n_in, m, tt] <- VaR_rg_n_rvol_5_in_s[1:n_in]
      VaR_training_data_mod[n_oos, m, tt] <- VaR_oos[1:(tt-1), m]
      
      ES_training_data_mod[1:n_in, m, tt] <- ES_rg_n_rvol_5_in_s[1:n_in]
      ES_training_data_mod[n_oos, m, tt] <- ES_oos[1:(tt-1), m]
      
    } else {
      # Only oos 
      VaR_training_data_mod[, m, tt] <- VaR_oos[(tt - Tin):(tt-1), m]
      ES_training_data_mod[, m, tt] <- ES_oos[(tt - Tin):(tt-1), m]
    }
  


############################################ 
############################################ RGARCH with a normal distribution (with rvol_5) and CF approx
############################################ M11

m<-11

z_cf_VaR<-Quantile_CF_f(r_t_in_s_tin,sigma_rg_n_rvol_5_in_s,tau)[,1]
z_cf_ES<-Quantile_CF_f(r_t_in_s_tin,sigma_rg_n_rvol_5_in_s,tau)[,2]

VaR_rg_n_cf_rvol_5_in_s<- sigma_rg_n_rvol_5_in_s*z_cf_VaR
ES_rg_n_cf_rvol_5_in_s<- sigma_rg_n_rvol_5_in_s*z_cf_ES

VaR_oos[tt,m]<- sigma_rg_n_rvol_5_oos*z_cf_VaR
ES_oos[tt,m]<- sigma_rg_n_rvol_5_oos*z_cf_ES

if (tt == 1) {
      # first step, only in-sample
      VaR_training_data_mod[, m, tt] <- VaR_rg_n_cf_rvol_5_in_s
      ES_training_data_mod[, m, tt] <- ES_rg_n_cf_rvol_5_in_s
      
    } else if (tt <= Tin) {
      # Mixed window: first in-sample, last oos
      n_in <- Tin - tt + 1
      n_oos <- (n_in+1):Tin
      
      VaR_training_data_mod[1:n_in, m, tt] <- VaR_rg_n_cf_rvol_5_in_s[1:n_in]
      VaR_training_data_mod[n_oos, m, tt] <- VaR_oos[1:(tt-1), m]
      
      ES_training_data_mod[1:n_in, m, tt] <- ES_rg_n_cf_rvol_5_in_s[1:n_in]
      ES_training_data_mod[n_oos, m, tt] <- ES_oos[1:(tt-1), m]
      
    } else {
      # Only oos 
      VaR_training_data_mod[, m, tt] <- VaR_oos[(tt - Tin):(tt-1), m]
      ES_training_data_mod[, m, tt] <- ES_oos[(tt - Tin):(tt-1), m]
    }
  

############################################ 
############################################ RGARCH with a student-t distribution (with rvol_5)
############################################ M12

m<-12


rg_est_rvol_5_t<-fit_robust(r_t_est_cycle, rvol_5_est_cycle, 
spec_re_garch_t, is_t = TRUE)$fit



sigma_rg_t_rvol_5_in_s <- rg_est_rvol_5_t@fit$sigma / 100
sigma_rg_t_rvol_5_oos  <- as.numeric(sigma(ugarchforecast(rg_est_rvol_5_t, n.ahead = 1))) / 100
df_hat <- as.numeric(rg_est_rvol_5_t@fit$coef['shape'])
correc <- ((df_hat - 2) / df_hat)^0.5



VaR_oos[tt,m]<-qt(tau,df_hat)*sigma_rg_t_rvol_5_oos*correc

## first term:
first_term<-(dt(qt(tau,df_hat),df_hat)/tau)

## second term:
second_term<- (df_hat+(qt(tau,df_hat))^2)/(df_hat-1)

ES_oos[tt,m]<- -sigma_rg_t_rvol_5_oos*first_term*second_term*correc

VaR_rg_t_rvol_5_in_s<-qt(tau,df_hat)*rg_est_rvol_5_t@fit$sigma*correc
ES_rg_t_rvol_5_in_s<- -rg_est_rvol_5_t@fit$sigma*first_term*second_term*correc

if (tt == 1) {
      # first step, only in-sample
      VaR_training_data_mod[, m, tt] <- VaR_rg_t_rvol_5_in_s
      ES_training_data_mod[, m, tt] <- ES_rg_t_rvol_5_in_s
      
    } else if (tt <= Tin) {
      # Mixed window: first in-sample, last oos
      n_in <- Tin - tt + 1
      n_oos <- (n_in+1):Tin
      
      VaR_training_data_mod[1:n_in, m, tt] <- VaR_rg_t_rvol_5_in_s[1:n_in]
      VaR_training_data_mod[n_oos, m, tt] <- VaR_oos[1:(tt-1), m]
      
      ES_training_data_mod[1:n_in, m, tt] <- ES_rg_t_rvol_5_in_s[1:n_in]
      ES_training_data_mod[n_oos, m, tt] <- ES_oos[1:(tt-1), m]
      
    } else {
      # Only oos 
      VaR_training_data_mod[, m, tt] <- VaR_oos[(tt - Tin):(tt-1), m]
      ES_training_data_mod[, m, tt] <- ES_oos[(tt - Tin):(tt-1), m]
    }
  


############################################ 
############################################ RGARCH with a normal distribution (with rb_ss)
############################################ M13

m<-13


rg_est_rb_ss_n<-fit_robust(r_t_est_cycle, rb_ss_est_cycle, spec_re_garch_n, is_t = FALSE)$fit


sigma_rg_n_rb_ss_in_s<-rg_est_rb_ss_n@fit$sigma/100
sigma_rg_n_rb_ss_oos<-as.numeric(sigma(ugarchforecast(rg_est_rb_ss_n,  n.ahead = 1)))/100

VaR_oos[tt,m]<-z_a*sigma_rg_n_rb_ss_oos
ES_oos[tt,m]<- -sigma_rg_n_rb_ss_oos*dnorm(z_a)/tau

VaR_rg_n_rb_ss_in_s<- z_a*sigma_rg_n_rb_ss_in_s
ES_rg_n_rb_ss_in_s<- -sigma_rg_n_rb_ss_in_s*dnorm(z_a)/tau

if (tt == 1) {
      # first step, only in-sample
      VaR_training_data_mod[, m, tt] <- VaR_rg_n_rb_ss_in_s
      ES_training_data_mod[, m, tt] <- ES_rg_n_rb_ss_in_s
      
    } else if (tt <= Tin) {
      # Mixed window: first in-sample, last oos
      n_in <- Tin - tt + 1
      n_oos <- (n_in+1):Tin
      
      VaR_training_data_mod[1:n_in, m, tt] <- VaR_rg_n_rb_ss_in_s[1:n_in]
      VaR_training_data_mod[n_oos, m, tt] <- VaR_oos[1:(tt-1), m]
      
      ES_training_data_mod[1:n_in, m, tt] <- ES_rg_n_rb_ss_in_s[1:n_in]
      ES_training_data_mod[n_oos, m, tt] <- ES_oos[1:(tt-1), m]
      
    } else {
      # Only oos 
      VaR_training_data_mod[, m, tt] <- VaR_oos[(tt - Tin):(tt-1), m]
      ES_training_data_mod[, m, tt] <- ES_oos[(tt - Tin):(tt-1), m]
    }
  


############################################ 
############################################ RGARCH with a normal distribution (with rb_ss) and CF approx
############################################ M14

m<-14

z_cf_VaR<-Quantile_CF_f(r_t_in_s_tin,sigma_rg_n_rb_ss_in_s,tau)[,1]
z_cf_ES<-Quantile_CF_f(r_t_in_s_tin,sigma_rg_n_rb_ss_in_s,tau)[,2]

VaR_rg_n_cf_rb_ss_in_s<- sigma_rg_n_rb_ss_in_s*z_cf_VaR
ES_rg_n_cf_rb_ss_in_s<- sigma_rg_n_rb_ss_in_s*z_cf_ES

VaR_oos[tt,m]<- sigma_rg_n_rb_ss_oos*z_cf_VaR
ES_oos[tt,m]<- sigma_rg_n_rb_ss_oos*z_cf_ES

if (tt == 1) {
      # first step, only in-sample
      VaR_training_data_mod[, m, tt] <- VaR_rg_n_cf_rb_ss_in_s
      ES_training_data_mod[, m, tt] <- ES_rg_n_cf_rb_ss_in_s
      
    } else if (tt <= Tin) {
      # Mixed window: first in-sample, last oos
      n_in <- Tin - tt + 1
      n_oos <- (n_in+1):Tin
      
      VaR_training_data_mod[1:n_in, m, tt] <- VaR_rg_n_cf_rb_ss_in_s[1:n_in]
      VaR_training_data_mod[n_oos, m, tt] <- VaR_oos[1:(tt-1), m]
      
      ES_training_data_mod[1:n_in, m, tt] <- ES_rg_n_cf_rb_ss_in_s[1:n_in]
      ES_training_data_mod[n_oos, m, tt] <- ES_oos[1:(tt-1), m]
      
    } else {
      # Only oos 
      VaR_training_data_mod[, m, tt] <- VaR_oos[(tt - Tin):(tt-1), m]
      ES_training_data_mod[, m, tt] <- ES_oos[(tt - Tin):(tt-1), m]
    }
  



############################################ 
############################################ RGARCH with a student-t distribution (with rb_ss)
############################################ M15

m<-15


rg_est_rb_ss_t<-fit_robust(r_t_est_cycle, rb_ss_est_cycle, spec_re_garch_t, is_t = TRUE)$fit

sigma_rg_t_rb_ss_in_s <- rg_est_rb_ss_t@fit$sigma / 100
sigma_rg_t_rb_ss_oos  <- as.numeric(sigma(ugarchforecast(rg_est_rb_ss_t, n.ahead = 1))) / 100
df_hat <- as.numeric(rg_est_rb_ss_t@fit$coef['shape'])
correc <- ((df_hat - 2) / df_hat)^0.5


VaR_oos[tt,m]<-qt(tau,df_hat)*sigma_rg_t_rb_ss_oos*correc

## first term:
first_term<-(dt(qt(tau,df_hat),df_hat)/tau)

## second term:
second_term<- (df_hat+(qt(tau,df_hat))^2)/(df_hat-1)

ES_oos[tt,m]<- -sigma_rg_t_rb_ss_oos*first_term*second_term*correc

VaR_rg_t_rb_ss_in_s<-qt(tau,df_hat)*rg_est_rb_ss_t@fit$sigma*correc
ES_rg_t_rb_ss_in_s<- -rg_est_rb_ss_t@fit$sigma*first_term*second_term*correc

if (tt == 1) {
      # first step, only in-sample
      VaR_training_data_mod[, m, tt] <- VaR_rg_t_rb_ss_in_s
      ES_training_data_mod[, m, tt] <- ES_rg_t_rb_ss_in_s
      
    } else if (tt <= Tin) {
      # Mixed window: first in-sample, last oos
      n_in <- Tin - tt + 1
      n_oos <- (n_in+1):Tin
      
      VaR_training_data_mod[1:n_in, m, tt] <- VaR_rg_t_rb_ss_in_s[1:n_in]
      VaR_training_data_mod[n_oos, m, tt] <- VaR_oos[1:(tt-1), m]
      
      ES_training_data_mod[1:n_in, m, tt] <- ES_rg_t_rb_ss_in_s[1:n_in]
      ES_training_data_mod[n_oos, m, tt] <- ES_oos[1:(tt-1), m]
      
    } else {
      # Only oos 
      VaR_training_data_mod[, m, tt] <- VaR_oos[(tt - Tin):(tt-1), m]
      ES_training_data_mod[, m, tt] <- ES_oos[(tt - Tin):(tt-1), m]
    }


############################################ 
############################################ RGARCH with a normal distribution (with rk)
############################################ M16

m<-16


rg_est_rk_n<-fit_robust(r_t_est_cycle, rk_est_cycle, spec_re_garch_n, is_t = FALSE)$fit

sigma_rg_n_rk_in_s <- rg_est_rk_n@fit$sigma / 100
sigma_rg_n_rk_oos  <- as.numeric(sigma(ugarchforecast(rg_est_rk_n, n.ahead = 1))) / 100



VaR_oos[tt,m]<-z_a*sigma_rg_n_rk_oos
ES_oos[tt,m]<- -sigma_rg_n_rk_oos*dnorm(z_a)/tau

VaR_rg_n_rk_in_s<- z_a*sigma_rg_n_rk_in_s
ES_rg_n_rk_in_s<- -sigma_rg_n_rk_in_s*dnorm(z_a)/tau


if (tt == 1) {
      # first step, only in-sample
      VaR_training_data_mod[, m, tt] <- VaR_rg_n_rk_in_s
      ES_training_data_mod[, m, tt] <- ES_rg_n_rk_in_s
      
    } else if (tt <= Tin) {
      # Mixed window: first in-sample, last oos
      n_in <- Tin - tt + 1
      n_oos <- (n_in+1):Tin
      
      VaR_training_data_mod[1:n_in, m, tt] <- VaR_rg_n_rk_in_s[1:n_in]
      VaR_training_data_mod[n_oos, m, tt] <- VaR_oos[1:(tt-1), m]
      
      ES_training_data_mod[1:n_in, m, tt] <- ES_rg_n_rk_in_s[1:n_in]
      ES_training_data_mod[n_oos, m, tt] <- ES_oos[1:(tt-1), m]
      
    } else {
      # Only oos 
      VaR_training_data_mod[, m, tt] <- VaR_oos[(tt - Tin):(tt-1), m]
      ES_training_data_mod[, m, tt] <- ES_oos[(tt - Tin):(tt-1), m]
    }
  



############################################ 
############################################ RGARCH with a normal distribution (with rk) and CF approx
############################################ M17

m<-17

z_cf_VaR<-Quantile_CF_f(r_t_in_s_tin,sigma_rg_n_rk_in_s,tau)[,1]
z_cf_ES<-Quantile_CF_f(r_t_in_s_tin,sigma_rg_n_rk_in_s,tau)[,2]

VaR_rg_n_cf_rk_in_s<- sigma_rg_n_rk_in_s*z_cf_VaR
ES_rg_n_cf_rk_in_s<- sigma_rg_n_rk_in_s*z_cf_ES

VaR_oos[tt,m]<- sigma_rg_n_rk_oos*z_cf_VaR
ES_oos[tt,m]<- sigma_rg_n_rk_oos*z_cf_ES


if (tt == 1) {
      # first step, only in-sample
      VaR_training_data_mod[, m, tt] <- VaR_rg_n_cf_rk_in_s
      ES_training_data_mod[, m, tt] <- ES_rg_n_cf_rk_in_s
      
    } else if (tt <= Tin) {
      # Mixed window: first in-sample, last oos
      n_in <- Tin - tt + 1
      n_oos <- (n_in+1):Tin
      
      VaR_training_data_mod[1:n_in, m, tt] <- VaR_rg_n_cf_rk_in_s[1:n_in]
      VaR_training_data_mod[n_oos, m, tt] <- VaR_oos[1:(tt-1), m]
      
      ES_training_data_mod[1:n_in, m, tt] <- ES_rg_n_cf_rk_in_s[1:n_in]
      ES_training_data_mod[n_oos, m, tt] <- ES_oos[1:(tt-1), m]
      
    } else {
      # Only oos 
      VaR_training_data_mod[, m, tt] <- VaR_oos[(tt - Tin):(tt-1), m]
      ES_training_data_mod[, m, tt] <- ES_oos[(tt - Tin):(tt-1), m]
    }
  


############################################ 
############################################ RGARCH with a student-t distribution (with rk)
############################################ M18

m<-18

rg_est_rk_t<-fit_robust(r_t_est_cycle, rk_est_cycle, spec_re_garch_t, is_t = TRUE)$fit

sigma_rg_t_rk_in_s <- rg_est_rk_t@fit$sigma / 100
sigma_rg_t_rk_oos  <- as.numeric(sigma(ugarchforecast(rg_est_rk_t, n.ahead = 1))) / 100
df_hat <- as.numeric(rg_est_rk_t@fit$coef['shape'])
correc <- ((df_hat - 2) / df_hat)^0.5

VaR_oos[tt,m]<-qt(tau,df_hat)*sigma_rg_t_rk_oos*correc

## first term:
first_term<-(dt(qt(tau,df_hat),df_hat)/tau)

## second term:
second_term<- (df_hat+(qt(tau,df_hat))^2)/(df_hat-1)

ES_oos[tt,m]<- -sigma_rg_t_rk_oos*first_term*second_term*correc

VaR_rg_t_rk_in_s<-qt(tau,df_hat)*rg_est_rk_t@fit$sigma*correc
ES_rg_t_rk_in_s<- -rg_est_rk_t@fit$sigma*first_term*second_term*correc

if (tt == 1) {
      # first step, only in-sample
      VaR_training_data_mod[, m, tt] <- VaR_rg_t_rk_in_s
      ES_training_data_mod[, m, tt] <- ES_rg_t_rk_in_s
      
    } else if (tt <= Tin) {
      # Mixed window: first in-sample, last oos
      n_in <- Tin - tt + 1
      n_oos <- (n_in+1):Tin
      
      VaR_training_data_mod[1:n_in, m, tt] <- VaR_rg_t_rk_in_s[1:n_in]
      VaR_training_data_mod[n_oos, m, tt] <- VaR_oos[1:(tt-1), m]
      
      ES_training_data_mod[1:n_in, m, tt] <- ES_rg_t_rk_in_s[1:n_in]
      ES_training_data_mod[n_oos, m, tt] <- ES_oos[1:(tt-1), m]
      
    } else {
      # Only oos 
      VaR_training_data_mod[, m, tt] <- VaR_oos[(tt - Tin):(tt-1), m]
      ES_training_data_mod[, m, tt] <- ES_oos[(tt - Tin):(tt-1), m]
    }
  

##################################################
################################################## Non-Parametric methods
################################################## (from M19 to M23)

############################################### HS 25 (M19)

m<-19

VaR_HS_25_in_s<-VaR_HS_25[index(r_t_in_s)[-(Tin+1)]]
ES_HS_25_in_s<-ES_HS_25[index(r_t_in_s)[-(Tin+1)]]

VaR_oos[tt,m]<-VaR_HS_25[index(r_t_oos)]
ES_oos[tt,m]<-ES_HS_25[index(r_t_oos)]

if (tt == 1) {
      # first step, only in-sample
      VaR_training_data_mod[, m, tt] <- VaR_HS_25_in_s
      ES_training_data_mod[, m, tt] <- ES_HS_25_in_s
      
    } else if (tt <= Tin) {
      # Mixed window: first in-sample, last oos
      n_in <- Tin - tt + 1
      n_oos <- (n_in+1):Tin
      
      VaR_training_data_mod[1:n_in, m, tt] <- VaR_HS_25_in_s[1:n_in]
      VaR_training_data_mod[n_oos, m, tt] <- VaR_oos[1:(tt-1), m]
      
      ES_training_data_mod[1:n_in, m, tt] <- ES_HS_25_in_s[1:n_in]
      ES_training_data_mod[n_oos, m, tt] <- ES_oos[1:(tt-1), m]
      
    } else {
      # Only oos 
      VaR_training_data_mod[, m, tt] <- VaR_oos[(tt - Tin):(tt-1), m]
      ES_training_data_mod[, m, tt] <- ES_oos[(tt - Tin):(tt-1), m]
    }
  

############################################### HS 50 (M20)

m<-20

VaR_HS_50_in_s<-VaR_HS_50[index(r_t_in_s)[-(Tin+1)]]
ES_HS_50_in_s<-ES_HS_50[index(r_t_in_s)[-(Tin+1)]]

VaR_oos[tt,m]<-VaR_HS_50[index(r_t_oos)]
ES_oos[tt,m]<-ES_HS_50[index(r_t_oos)]

if (tt == 1) {
      # first step, only in-sample
      VaR_training_data_mod[, m, tt] <- VaR_HS_50_in_s
      ES_training_data_mod[, m, tt] <- ES_HS_50_in_s
      
    } else if (tt <= Tin) {
      # Mixed window: first in-sample, last oos
      n_in <- Tin - tt + 1
      n_oos <- (n_in+1):Tin
      
      VaR_training_data_mod[1:n_in, m, tt] <- VaR_HS_50_in_s[1:n_in]
      VaR_training_data_mod[n_oos, m, tt] <- VaR_oos[1:(tt-1), m]
      
      ES_training_data_mod[1:n_in, m, tt] <- ES_HS_50_in_s[1:n_in]
      ES_training_data_mod[n_oos, m, tt] <- ES_oos[1:(tt-1), m]
      
    } else {
      # Only oos 
      VaR_training_data_mod[, m, tt] <- VaR_oos[(tt - Tin):(tt-1), m]
      ES_training_data_mod[, m, tt] <- ES_oos[(tt - Tin):(tt-1), m]
    }
  

############################################### HS 100 (M21)

m<-21

VaR_HS_100_in_s<-VaR_HS_100[index(r_t_in_s)[-(Tin+1)]]
ES_HS_100_in_s<-ES_HS_100[index(r_t_in_s)[-(Tin+1)]]

VaR_oos[tt,m]<-VaR_HS_100[index(r_t_oos)]
ES_oos[tt,m]<-ES_HS_100[index(r_t_oos)]

if (tt == 1) {
      # first step, only in-sample
      VaR_training_data_mod[, m, tt] <- VaR_HS_100_in_s
      ES_training_data_mod[, m, tt] <- ES_HS_100_in_s
      
    } else if (tt <= Tin) {
      # Mixed window: first in-sample, last oos
      n_in <- Tin - tt + 1
      n_oos <- (n_in+1):Tin
      
      VaR_training_data_mod[1:n_in, m, tt] <- VaR_HS_100_in_s[1:n_in]
      VaR_training_data_mod[n_oos, m, tt] <- VaR_oos[1:(tt-1), m]
      
      ES_training_data_mod[1:n_in, m, tt] <- ES_HS_100_in_s[1:n_in]
      ES_training_data_mod[n_oos, m, tt] <- ES_oos[1:(tt-1), m]
      
    } else {
      # Only oos 
      VaR_training_data_mod[, m, tt] <- VaR_oos[(tt - Tin):(tt-1), m]
      ES_training_data_mod[, m, tt] <- ES_oos[(tt - Tin):(tt-1), m]
    }
  

############################################### HS 250 (M22)

m<-22

VaR_HS_250_in_s<-VaR_HS_250[index(r_t_in_s)[-(Tin+1)]]
ES_HS_250_in_s<-ES_HS_250[index(r_t_in_s)[-(Tin+1)]]

VaR_oos[tt,m]<-VaR_HS_250[index(r_t_oos)]
ES_oos[tt,m]<-ES_HS_250[index(r_t_oos)]

if (tt == 1) {
      # first step, only in-sample
      VaR_training_data_mod[, m, tt] <- VaR_HS_250_in_s
      ES_training_data_mod[, m, tt] <- ES_HS_250_in_s
      
    } else if (tt <= Tin) {
      # Mixed window: first in-sample, last oos
      n_in <- Tin - tt + 1
      n_oos <- (n_in+1):Tin
      
      VaR_training_data_mod[1:n_in, m, tt] <- VaR_HS_250_in_s[1:n_in]
      VaR_training_data_mod[n_oos, m, tt] <- VaR_oos[1:(tt-1), m]
      
      ES_training_data_mod[1:n_in, m, tt] <- ES_HS_250_in_s[1:n_in]
      ES_training_data_mod[n_oos, m, tt] <- ES_oos[1:(tt-1), m]
      
    } else {
      # Only oos 
      VaR_training_data_mod[, m, tt] <- VaR_oos[(tt - Tin):(tt-1), m]
      ES_training_data_mod[, m, tt] <- ES_oos[(tt - Tin):(tt-1), m]
    }
  

############################################### HS 500 (M23)

m<-23

VaR_HS_500_in_s<-VaR_HS_500[index(r_t_in_s)[-(Tin+1)]]
ES_HS_500_in_s<-ES_HS_500[index(r_t_in_s)[-(Tin+1)]]

VaR_oos[tt,m]<-VaR_HS_500[index(r_t_oos)]
ES_oos[tt,m]<-ES_HS_500[index(r_t_oos)]


if (tt == 1) {
      # first step, only in-sample
      VaR_training_data_mod[, m, tt] <- VaR_HS_500_in_s
      ES_training_data_mod[, m, tt] <- ES_HS_500_in_s
      
    } else if (tt <= Tin) {
      # Mixed window: first in-sample, last oos
      n_in <- Tin - tt + 1
      n_oos <- (n_in+1):Tin
      
      VaR_training_data_mod[1:n_in, m, tt] <- VaR_HS_500_in_s[1:n_in]
      VaR_training_data_mod[n_oos, m, tt] <- VaR_oos[1:(tt-1), m]
      
      ES_training_data_mod[1:n_in, m, tt] <- ES_HS_500_in_s[1:n_in]
      ES_training_data_mod[n_oos, m, tt] <- ES_oos[1:(tt-1), m]
      
    } else {
      # Only oos 
      VaR_training_data_mod[, m, tt] <- VaR_oos[(tt - Tin):(tt-1), m]
      ES_training_data_mod[, m, tt] <- ES_oos[(tt - Tin):(tt-1), m]
    }
  

##################################################
################################################## Semi-Parametric methods
################################################## (from M24 to M32)

############################################
############################################ SAV
############################################ M24

m <- 24

############################################
# Estimate the SAV model
############################################

fit_sav_es <- tryCatch(
  ucaviarfit(
    model = "SAV",
    tau,
    daily_ret = r_t_in_s_tin_xts,
    R = 1000,
    B = 1,
    Exp_Short = "Yes",
    std_err = "not_boot"
  ),
  error = function(e) NULL
)

############################################
# Check whether the estimation is valid
############################################

valid_fit_sav <- (
  !is.null(fit_sav_es) &&
  !is.null(fit_sav_es$coef_mat) &&
  NCOL(fit_sav_es$coef_mat) >= 1 &&
  all(is.finite(fit_sav_es$coef_mat[, 1])) &&
  !is.null(fit_sav_es$VaR) &&
  !is.null(fit_sav_es$ES) &&
  all(is.finite(fit_sav_es$VaR)) &&
  all(is.finite(fit_sav_es$ES))
)

if (valid_fit_sav) {

  ##########################################
  # Successful estimation
  ##########################################

  coef_sav <- fit_sav_es$coef_mat[, 1]

  VaR_SAV_in_s <- fit_sav_es$VaR
  ES_SAV_in_s  <- fit_sav_es$ES

} else {

  ##########################################
  # Failed estimation:
  # use the coefficients from the previous
  # point and reconstruct the VaR and ES
  # paths on the updated rolling window
  ##########################################

  if (
    is.null(coef_sav_previous) ||
    is.null(VaR_SAV_in_s_previous)
  ) {
    stop(
      paste0(
        "SAV estimation failed at tt = ",
        tt,
        " and no previous valid estimates are available."
      )
    )
  }

  failed_SAV <- c(failed_SAV, tt)

  # Use the coefficients estimated at the previous point
  coef_sav <- coef_sav_previous

  # Current updated in-sample return window
  r_sav_current <- as.numeric(r_t_in_s_tin_xts)
  n_sav <- length(r_sav_current)

  # Initialize the new in-sample VaR path
  VaR_SAV_in_s <- rep(NA_real_, n_sav)

  # The first observation of the updated rolling window
  # corresponds to the second observation of the previous window
  VaR_SAV_in_s[1] <- VaR_SAV_in_s_previous[2]

  # Reconstruct the in-sample VaR path
  for (i in 2:n_sav) {

    VaR_SAV_in_s[i] <- as.numeric(
      coef_sav[1] +
      coef_sav[2] * VaR_SAV_in_s[i - 1] +
      coef_sav[3] * abs(r_sav_current[i - 1])
    )
  }

  # Reconstruct the corresponding in-sample ES path
  ES_SAV_in_s <- -abs(
    (1 + exp(coef_sav[4])) * VaR_SAV_in_s
  )

  warning(
    paste0(
      "SAV estimation failed at tt = ",
      tt,
      ". Previous coefficients were used."
    ),
    call. = FALSE
  )
}

############################################
# One-step-ahead VaR and ES forecasts
############################################

fit_sav_es_VaR_oos <- as.numeric(
  coef_sav[1] +
  coef_sav[2] * last(VaR_SAV_in_s) +
  coef_sav[3] * last(abs(r_t_in_s_tin_xts))
)

fit_sav_es_ES_oos <- -abs(
  (1 + exp(coef_sav[4])) * fit_sav_es_VaR_oos
)

VaR_oos[tt, m] <- fit_sav_es_VaR_oos
ES_oos[tt, m]  <- fit_sav_es_ES_oos

############################################
# Construct the rolling training data
############################################

if (tt == 1) {

  # First step: only in-sample observations
  VaR_training_data_mod[, m, tt] <- VaR_SAV_in_s
  ES_training_data_mod[, m, tt]  <- ES_SAV_in_s

} else if (tt <= Tin) {

  # Mixed window: first in-sample, last out-of-sample
  n_in  <- Tin - tt + 1
  n_oos <- (n_in + 1):Tin

  VaR_training_data_mod[1:n_in, m, tt] <-
    VaR_SAV_in_s[1:n_in]

  VaR_training_data_mod[n_oos, m, tt] <-
    VaR_oos[1:(tt - 1), m]

  ES_training_data_mod[1:n_in, m, tt] <-
    ES_SAV_in_s[1:n_in]

  ES_training_data_mod[n_oos, m, tt] <-
    ES_oos[1:(tt - 1), m]

} else {

  # Only out-of-sample observations
  VaR_training_data_mod[, m, tt] <-
    VaR_oos[(tt - Tin):(tt - 1), m]

  ES_training_data_mod[, m, tt] <-
    ES_oos[(tt - Tin):(tt - 1), m]
}

############################################
# Store the current estimates for the
# following rolling estimation point
############################################

coef_sav_previous <- coef_sav
VaR_SAV_in_s_previous <- VaR_SAV_in_s
ES_SAV_in_s_previous <- ES_SAV_in_s


############################################
############################################ AS
############################################ M25

m <- 25

############################################
# Estimate the AS model
############################################

fit_as_es <- tryCatch(
  ucaviarfit(
    model = "AS",
    tau,
    daily_ret = r_t_in_s_tin_xts,
    R = 1000,
    B = 1,
    Exp_Short = "Yes",
    std_err = "not_boot"
  ),
  error = function(e) NULL
)

############################################
# Check whether the estimation is valid
############################################

valid_fit_as <- (
  !is.null(fit_as_es) &&
  !is.null(fit_as_es$coef_mat) &&
  NCOL(fit_as_es$coef_mat) >= 1 &&
  all(is.finite(fit_as_es$coef_mat[, 1])) &&
  !is.null(fit_as_es$VaR) &&
  !is.null(fit_as_es$ES) &&
  all(is.finite(fit_as_es$VaR)) &&
  all(is.finite(fit_as_es$ES))
)

if (valid_fit_as) {

  ##########################################
  # Successful estimation
  ##########################################

  coef_as <- fit_as_es$coef_mat[, 1]

  VaR_AS_in_s <- as.numeric(fit_as_es$VaR)
  ES_AS_in_s  <- as.numeric(fit_as_es$ES)

} else {

  ##########################################
  # Failed estimation:
  # use the coefficients from the previous
  # point and reconstruct the VaR and ES
  # paths on the updated rolling window
  ##########################################

  if (
    is.null(coef_as_previous) ||
    is.null(VaR_AS_in_s_previous)
  ) {
    stop(
      paste0(
        "AS estimation failed at tt = ",
        tt,
        " and no previous valid estimates are available."
      )
    )
  }

  failed_AS <- c(failed_AS, tt)

  # Use the coefficients estimated at the previous point
  coef_as <- coef_as_previous

  # Current updated in-sample return window
  r_as_current <- as.numeric(r_t_in_s_tin_xts)
  n_as <- length(r_as_current)

  # Initialize the new in-sample VaR path
  VaR_AS_in_s <- rep(NA_real_, n_as)

  # The first observation of the updated rolling window
  # corresponds to the second observation of the previous window
  VaR_AS_in_s[1] <- VaR_AS_in_s_previous[2]

  # Reconstruct the in-sample VaR path
  for (i in 2:n_as) {

    ret_previous <- r_as_current[i - 1]

    ind_pos_previous <- ifelse(ret_previous > 0, 1, 0)
    ind_neg_previous <- ifelse(ret_previous < 0, 1, 0)

    VaR_AS_in_s[i] <- as.numeric(
      coef_as[1] +
      coef_as[2] * VaR_AS_in_s[i - 1] +
      (
        coef_as[3] * ind_pos_previous +
        coef_as[4] * ind_neg_previous
      ) * abs(ret_previous)
    )
  }

  # Reconstruct the corresponding in-sample ES path
  ES_AS_in_s <- -abs(
    (1 + exp(coef_as[5])) * VaR_AS_in_s
  )

  warning(
    paste0(
      "AS estimation failed at tt = ",
      tt,
      ". Previous coefficients were used."
    ),
    call. = FALSE
  )
}

############################################
# One-step-ahead VaR and ES forecasts
############################################

ret <- as.numeric(last(r_t_in_s_tin_xts))

ind_pos <- ifelse(ret > 0, 1, 0)
ind_neg <- ifelse(ret < 0, 1, 0)

fit_as_es_VaR_oos <- as.numeric(
  coef_as[1] +
  coef_as[2] * last(VaR_AS_in_s) +
  (
    coef_as[3] * ind_pos +
    coef_as[4] * ind_neg
  ) * abs(ret)
)

fit_as_es_ES_oos <- -abs(
  (1 + exp(coef_as[5])) * fit_as_es_VaR_oos
)

VaR_oos[tt, m] <- fit_as_es_VaR_oos
ES_oos[tt, m]  <- fit_as_es_ES_oos

############################################
# Construct the rolling training data
############################################

if (tt == 1) {

  # First step: only in-sample observations
  VaR_training_data_mod[, m, tt] <- VaR_AS_in_s
  ES_training_data_mod[, m, tt]  <- ES_AS_in_s

} else if (tt <= Tin) {

  # Mixed window: first in-sample, last out-of-sample
  n_in  <- Tin - tt + 1
  n_oos <- (n_in + 1):Tin

  VaR_training_data_mod[1:n_in, m, tt] <-
    VaR_AS_in_s[1:n_in]

  VaR_training_data_mod[n_oos, m, tt] <-
    VaR_oos[1:(tt - 1), m]

  ES_training_data_mod[1:n_in, m, tt] <-
    ES_AS_in_s[1:n_in]

  ES_training_data_mod[n_oos, m, tt] <-
    ES_oos[1:(tt - 1), m]

} else {

  # Only out-of-sample observations
  VaR_training_data_mod[, m, tt] <-
    VaR_oos[(tt - Tin):(tt - 1), m]

  ES_training_data_mod[, m, tt] <-
    ES_oos[(tt - Tin):(tt - 1), m]
}

############################################
# Store the current estimates for the
# following rolling estimation point
############################################

coef_as_previous <- coef_as
VaR_AS_in_s_previous <- VaR_AS_in_s
ES_AS_in_s_previous <- ES_AS_in_s

############################################ 
############################################ Indirect GARCH
############################################ M26

m <- 26

############################################
# Estimate the IG model
############################################

fit_ig_es <- tryCatch(
  ucaviarfit(
    model = "IG",
    tau,
    daily_ret = r_t_in_s_tin_xts,
    R = 1000,
    B = 1,
    Exp_Short = "Yes",
    std_err = "not_boot"
  ),
  error = function(e) NULL
)

############################################
# Check whether the estimation is valid
############################################

valid_fit_ig <- (
  !is.null(fit_ig_es) &&
  !is.null(fit_ig_es$coef_mat) &&
  NCOL(fit_ig_es$coef_mat) >= 1 &&
  all(is.finite(fit_ig_es$coef_mat[, 1])) &&
  !is.null(fit_ig_es$VaR) &&
  !is.null(fit_ig_es$ES) &&
  all(is.finite(fit_ig_es$VaR)) &&
  all(is.finite(fit_ig_es$ES))
)

if (valid_fit_ig) {

  coef_ig_candidate <- fit_ig_es$coef_mat[, 1]

  ##########################################
  # Check the one-step-ahead radicand
  ##########################################

  radicand_ig_oos <- as.numeric(
    coef_ig_candidate[1] +
    coef_ig_candidate[2] * last(fit_ig_es$VaR)^2 +
    coef_ig_candidate[3] * last(r_t_in_s_tin_xts)^2
  )

  valid_fit_ig <- (
    length(radicand_ig_oos) == 1 &&
    is.finite(radicand_ig_oos) &&
    radicand_ig_oos >= 0
  )
}

if (valid_fit_ig) {

  ##########################################
  # Successful estimation
  ##########################################

  coef_ig <- coef_ig_candidate

  VaR_IG_in_s <- as.numeric(fit_ig_es$VaR)
  ES_IG_in_s  <- as.numeric(fit_ig_es$ES)

} else {

  ##########################################
  # Failed estimation or invalid forecast:
  # use the coefficients from the previous
  # point and reconstruct the VaR and ES
  # paths on the updated rolling window
  ##########################################

  if (
    is.null(coef_ig_previous) ||
    is.null(VaR_IG_in_s_previous)
  ) {
    stop(
      paste0(
        "IG estimation failed at tt = ",
        tt,
        " and no previous valid estimates are available."
      )
    )
  }

  failed_IG <- c(failed_IG, tt)

  # Use the coefficients from the previous point
  coef_ig <- coef_ig_previous

  # Current updated in-sample return window
  r_ig_current <- as.numeric(r_t_in_s_tin_xts)
  n_ig <- length(r_ig_current)

  # Initialize the new in-sample VaR path
  VaR_IG_in_s <- rep(NA_real_, n_ig)

  # The first observation of the updated window corresponds
  # to the second observation of the previous window
  VaR_IG_in_s[1] <- VaR_IG_in_s_previous[2]

  ##########################################
  # Reconstruct the in-sample VaR path
  ##########################################

  for (i in 2:n_ig) {

    radicand_ig_in_s <- as.numeric(
      coef_ig[1] +
      coef_ig[2] * VaR_IG_in_s[i - 1]^2 +
      coef_ig[3] * r_ig_current[i - 1]^2
    )

    if (
      length(radicand_ig_in_s) != 1 ||
      !is.finite(radicand_ig_in_s) ||
      radicand_ig_in_s < 0
    ) {
      stop(
        paste0(
          "Invalid IG in-sample radicand at tt = ",
          tt,
          ", position i = ",
          i,
          ", using the previous coefficients."
        )
      )
    }

    VaR_IG_in_s[i] <- -sqrt(radicand_ig_in_s)
  }

  # Reconstruct the corresponding in-sample ES path
  ES_IG_in_s <- -abs(
    (1 + exp(coef_ig[4])) * VaR_IG_in_s
  )

  warning(
    paste0(
      "IG estimation failed or produced an invalid forecast at tt = ",
      tt,
      ". Previous coefficients were used."
    ),
    call. = FALSE
  )
}

############################################
# One-step-ahead VaR and ES forecasts
############################################

radicand_ig_oos <- as.numeric(
  coef_ig[1] +
  coef_ig[2] * last(VaR_IG_in_s)^2 +
  coef_ig[3] * last(r_t_in_s_tin_xts)^2
)

if (
  length(radicand_ig_oos) != 1 ||
  !is.finite(radicand_ig_oos) ||
  radicand_ig_oos < 0
) {
  stop(
    paste0(
      "Invalid IG out-of-sample radicand at tt = ",
      tt,
      "."
    )
  )
}

fit_ig_es_VaR_oos <- as.numeric(
  -sqrt(radicand_ig_oos)
)

fit_ig_es_ES_oos <- -abs(
  (1 + exp(coef_ig[4])) * fit_ig_es_VaR_oos
)

VaR_oos[tt, m] <- fit_ig_es_VaR_oos
ES_oos[tt, m]  <- fit_ig_es_ES_oos

############################################
# Construct the rolling training data
############################################

if (tt == 1) {

  # First step: only in-sample observations
  VaR_training_data_mod[, m, tt] <- VaR_IG_in_s
  ES_training_data_mod[, m, tt]  <- ES_IG_in_s

} else if (tt <= Tin) {

  # Mixed window: first in-sample, last out-of-sample
  n_in  <- Tin - tt + 1
  n_oos <- (n_in + 1):Tin

  VaR_training_data_mod[1:n_in, m, tt] <-
    VaR_IG_in_s[1:n_in]

  VaR_training_data_mod[n_oos, m, tt] <-
    VaR_oos[1:(tt - 1), m]

  ES_training_data_mod[1:n_in, m, tt] <-
    ES_IG_in_s[1:n_in]

  ES_training_data_mod[n_oos, m, tt] <-
    ES_oos[1:(tt - 1), m]

} else {

  # Only out-of-sample observations
  VaR_training_data_mod[, m, tt] <-
    VaR_oos[(tt - Tin):(tt - 1), m]

  ES_training_data_mod[, m, tt] <-
    ES_oos[(tt - Tin):(tt - 1), m]
}

############################################
# Store the current estimates for the
# following rolling estimation point
############################################

coef_ig_previous <- coef_ig
VaR_IG_in_s_previous <- VaR_IG_in_s
ES_IG_in_s_previous <- ES_IG_in_s
 
############################################
############################################ SAV-X: rvol_5
############################################ M27

m <- 27

############################################
# Updated in-sample X window
############################################

X_sav_rvol_5_current <- as.numeric(
  rvol_5_est_cycle[-(Tin + 1)]
)

############################################
# Estimate the SAV-X model
############################################

fit_sav_es_rvol_5 <- tryCatch(
  ucaviarfit(
    model = "SAVX",
    tau,
    daily_ret = r_t_in_s_tin_xts,
    X = X_sav_rvol_5_current,
    R = 1000,
    B = 1,
    Exp_Short = "Yes",
    std_err = "not_boot"
  ),
  error = function(e) NULL
)

############################################
# Check whether the estimation is valid
############################################

valid_fit_sav_rvol_5 <- (
  !is.null(fit_sav_es_rvol_5) &&
  !is.null(fit_sav_es_rvol_5$coef_mat) &&
  NCOL(fit_sav_es_rvol_5$coef_mat) >= 1 &&
  all(is.finite(fit_sav_es_rvol_5$coef_mat[, 1])) &&
  !is.null(fit_sav_es_rvol_5$VaR) &&
  !is.null(fit_sav_es_rvol_5$ES) &&
  all(is.finite(fit_sav_es_rvol_5$VaR)) &&
  all(is.finite(fit_sav_es_rvol_5$ES))
)

if (valid_fit_sav_rvol_5) {

  ##########################################
  # Successful estimation
  ##########################################

  coef_sav_rvol_5 <-
    fit_sav_es_rvol_5$coef_mat[, 1]

  VaR_SAVX_rvol_5_in_s <-
    as.numeric(fit_sav_es_rvol_5$VaR)

  ES_SAVX_rvol_5_in_s <-
    as.numeric(fit_sav_es_rvol_5$ES)

} else {

  ##########################################
  # Failed estimation:
  # use the coefficients from the previous
  # point and reconstruct VaR and ES on the
  # updated X window
  ##########################################

  if (
    is.null(coef_sav_rvol_5_previous) ||
    is.null(VaR_SAVX_rvol_5_in_s_previous)
  ) {
    stop(
      paste0(
        "SAV-X rvol_5 estimation failed at tt = ",
        tt,
        " and no previous valid estimates are available."
      )
    )
  }

  failed_SAVX_rvol_5 <-
    c(failed_SAVX_rvol_5, tt)

  # Use the coefficients estimated at the previous point
  coef_sav_rvol_5 <-
    coef_sav_rvol_5_previous

  n_sav_rvol_5 <- length(X_sav_rvol_5_current)

  if (n_sav_rvol_5 != length(r_t_in_s_tin_xts)) {
    stop(
      paste0(
        "The lengths of X and daily_ret differ at tt = ",
        tt,
        "."
      )
    )
  }

  if (any(!is.finite(X_sav_rvol_5_current))) {
    stop(
      paste0(
        "Non-finite values found in the SAV-X rvol_5 ",
        "in-sample window at tt = ",
        tt,
        "."
      )
    )
  }

  # Initialize the new in-sample VaR path
  VaR_SAVX_rvol_5_in_s <-
    rep(NA_real_, n_sav_rvol_5)

  # The first observation of the updated rolling window
  # corresponds to the second observation of the previous window
  VaR_SAVX_rvol_5_in_s[1] <-
    VaR_SAVX_rvol_5_in_s_previous[2]

  ##########################################
  # Reconstruct the in-sample VaR path
  ##########################################

  for (i in 2:n_sav_rvol_5) {

    VaR_SAVX_rvol_5_in_s[i] <- as.numeric(
      coef_sav_rvol_5[1] +
      coef_sav_rvol_5[2] *
        VaR_SAVX_rvol_5_in_s[i - 1] +
      coef_sav_rvol_5[3] *
        X_sav_rvol_5_current[i - 1]
    )
  }

  ##########################################
  # Reconstruct the in-sample ES path
  ##########################################

  ES_SAVX_rvol_5_in_s <- -abs(
    (1 + exp(coef_sav_rvol_5[4])) *
      VaR_SAVX_rvol_5_in_s
  )

  warning(
    paste0(
      "SAV-X rvol_5 estimation failed at tt = ",
      tt,
      ". Previous coefficients were used."
    ),
    call. = FALSE
  )
}

############################################
# One-step-ahead VaR and ES forecasts
############################################

fit_sav_es_VaR_oos_rvol_5 <- as.numeric(
  coef_sav_rvol_5[1] +
  coef_sav_rvol_5[2] *
    last(VaR_SAVX_rvol_5_in_s) +
  coef_sav_rvol_5[3] *
    last(X_sav_rvol_5_current)
)

fit_sav_es_ES_oos_rvol_5 <- -abs(
  (1 + exp(coef_sav_rvol_5[4])) *
    fit_sav_es_VaR_oos_rvol_5
)

VaR_oos[tt, m] <-
  fit_sav_es_VaR_oos_rvol_5

ES_oos[tt, m] <-
  fit_sav_es_ES_oos_rvol_5

############################################
# Construct the rolling training data
############################################

if (tt == 1) {

  # First step: only in-sample observations
  VaR_training_data_mod[, m, tt] <-
    VaR_SAVX_rvol_5_in_s

  ES_training_data_mod[, m, tt] <-
    ES_SAVX_rvol_5_in_s

} else if (tt <= Tin) {

  # Mixed window: first in-sample, last out-of-sample
  n_in  <- Tin - tt + 1
  n_oos <- (n_in + 1):Tin

  VaR_training_data_mod[1:n_in, m, tt] <-
    VaR_SAVX_rvol_5_in_s[1:n_in]

  VaR_training_data_mod[n_oos, m, tt] <-
    VaR_oos[1:(tt - 1), m]

  ES_training_data_mod[1:n_in, m, tt] <-
    ES_SAVX_rvol_5_in_s[1:n_in]

  ES_training_data_mod[n_oos, m, tt] <-
    ES_oos[1:(tt - 1), m]

} else {

  # Only out-of-sample observations
  VaR_training_data_mod[, m, tt] <-
    VaR_oos[(tt - Tin):(tt - 1), m]

  ES_training_data_mod[, m, tt] <-
    ES_oos[(tt - Tin):(tt - 1), m]
}

############################################
# Store the current estimates for the
# following rolling estimation point
############################################

coef_sav_rvol_5_previous <-
  coef_sav_rvol_5

VaR_SAVX_rvol_5_in_s_previous <-
  VaR_SAVX_rvol_5_in_s

ES_SAVX_rvol_5_in_s_previous <-
  ES_SAVX_rvol_5_in_s

############################################
############################################ SAV-X: rb_ss
############################################ M28

m <- 28

############################################
# Updated in-sample X window
############################################

X_sav_rb_ss_current <- as.numeric(
  rb_ss_est_cycle[-(Tin + 1)]
)

############################################
# Estimate the SAV-X model
############################################

fit_sav_es_rb_ss <- tryCatch(
  ucaviarfit(
    model = "SAVX",
    tau,
    daily_ret = r_t_in_s_tin_xts,
    X = X_sav_rb_ss_current,
    R = 1000,
    B = 1,
    Exp_Short = "Yes",
    std_err = "not_boot"
  ),
  error = function(e) NULL
)

############################################
# Check whether the estimation is valid
############################################

valid_fit_sav_rb_ss <- (
  !is.null(fit_sav_es_rb_ss) &&
  !is.null(fit_sav_es_rb_ss$coef_mat) &&
  NCOL(fit_sav_es_rb_ss$coef_mat) >= 1 &&
  all(is.finite(fit_sav_es_rb_ss$coef_mat[, 1])) &&
  !is.null(fit_sav_es_rb_ss$VaR) &&
  !is.null(fit_sav_es_rb_ss$ES) &&
  all(is.finite(fit_sav_es_rb_ss$VaR)) &&
  all(is.finite(fit_sav_es_rb_ss$ES))
)

if (valid_fit_sav_rb_ss) {

  ##########################################
  # Successful estimation
  ##########################################

  coef_sav_rb_ss <-
    fit_sav_es_rb_ss$coef_mat[, 1]

  VaR_SAVX_rb_ss_in_s <-
    as.numeric(fit_sav_es_rb_ss$VaR)

  ES_SAVX_rb_ss_in_s <-
    as.numeric(fit_sav_es_rb_ss$ES)

} else {

  ##########################################
  # Failed estimation:
  # use the coefficients from the previous
  # point and reconstruct VaR and ES on the
  # updated X window
  ##########################################

  if (
    is.null(coef_sav_rb_ss_previous) ||
    is.null(VaR_SAVX_rb_ss_in_s_previous)
  ) {
    stop(
      paste0(
        "SAV-X rb_ss estimation failed at tt = ",
        tt,
        " and no previous valid estimates are available."
      )
    )
  }

  failed_SAVX_rb_ss <-
    c(failed_SAVX_rb_ss, tt)

  # Use the coefficients estimated at the previous point
  coef_sav_rb_ss <-
    coef_sav_rb_ss_previous

  n_sav_rb_ss <- length(X_sav_rb_ss_current)

  if (n_sav_rb_ss != length(r_t_in_s_tin_xts)) {
    stop(
      paste0(
        "The lengths of X and daily_ret differ at tt = ",
        tt,
        "."
      )
    )
  }

  if (any(!is.finite(X_sav_rb_ss_current))) {
    stop(
      paste0(
        "Non-finite values found in the SAV-X rb_ss ",
        "in-sample window at tt = ",
        tt,
        "."
      )
    )
  }

  # Initialize the new in-sample VaR path
  VaR_SAVX_rb_ss_in_s <-
    rep(NA_real_, n_sav_rb_ss)

  # The first observation of the updated rolling window
  # corresponds to the second observation of the previous window
  VaR_SAVX_rb_ss_in_s[1] <-
    VaR_SAVX_rb_ss_in_s_previous[2]

  ##########################################
  # Reconstruct the in-sample VaR path
  ##########################################

  for (i in 2:n_sav_rb_ss) {

    VaR_SAVX_rb_ss_in_s[i] <- as.numeric(
      coef_sav_rb_ss[1] +
      coef_sav_rb_ss[2] *
        VaR_SAVX_rb_ss_in_s[i - 1] +
      coef_sav_rb_ss[3] *
        X_sav_rb_ss_current[i - 1]
    )
  }

  ##########################################
  # Reconstruct the in-sample ES path
  ##########################################

  ES_SAVX_rb_ss_in_s <- -abs(
    (1 + exp(coef_sav_rb_ss[4])) *
      VaR_SAVX_rb_ss_in_s
  )

  warning(
    paste0(
      "SAV-X rb_ss estimation failed at tt = ",
      tt,
      ". Previous coefficients were used."
    ),
    call. = FALSE
  )
}

############################################
# One-step-ahead VaR and ES forecasts
############################################

fit_sav_es_VaR_oos_rb_ss <- as.numeric(
  coef_sav_rb_ss[1] +
  coef_sav_rb_ss[2] *
    last(VaR_SAVX_rb_ss_in_s) +
  coef_sav_rb_ss[3] *
    last(X_sav_rb_ss_current)
)

fit_sav_es_ES_oos_rb_ss <- -abs(
  (1 + exp(coef_sav_rb_ss[4])) *
    fit_sav_es_VaR_oos_rb_ss
)

VaR_oos[tt, m] <-
  fit_sav_es_VaR_oos_rb_ss

ES_oos[tt, m] <-
  fit_sav_es_ES_oos_rb_ss

############################################
# Construct the rolling training data
############################################

if (tt == 1) {

  # First step: only in-sample observations
  VaR_training_data_mod[, m, tt] <-
    VaR_SAVX_rb_ss_in_s

  ES_training_data_mod[, m, tt] <-
    ES_SAVX_rb_ss_in_s

} else if (tt <= Tin) {

  # Mixed window: first in-sample, last out-of-sample
  n_in  <- Tin - tt + 1
  n_oos <- (n_in + 1):Tin

  VaR_training_data_mod[1:n_in, m, tt] <-
    VaR_SAVX_rb_ss_in_s[1:n_in]

  VaR_training_data_mod[n_oos, m, tt] <-
    VaR_oos[1:(tt - 1), m]

  ES_training_data_mod[1:n_in, m, tt] <-
    ES_SAVX_rb_ss_in_s[1:n_in]

  ES_training_data_mod[n_oos, m, tt] <-
    ES_oos[1:(tt - 1), m]

} else {

  # Only out-of-sample observations
  VaR_training_data_mod[, m, tt] <-
    VaR_oos[(tt - Tin):(tt - 1), m]

  ES_training_data_mod[, m, tt] <-
    ES_oos[(tt - Tin):(tt - 1), m]
}

############################################
# Store the current estimates for the
# following rolling estimation point
############################################

coef_sav_rb_ss_previous <-
  coef_sav_rb_ss

VaR_SAVX_rb_ss_in_s_previous <-
  VaR_SAVX_rb_ss_in_s

ES_SAVX_rb_ss_in_s_previous <-
  ES_SAVX_rb_ss_in_s

############################################ 
############################################ SAV+RK
############################################ M29

m <- 29

############################################
# Updated in-sample X window
############################################

X_sav_rk_current <- as.numeric(
  rk_est_cycle[-(Tin + 1)]
)

############################################
# Estimate the SAV-X model
############################################

fit_sav_es_rk <- tryCatch(
  ucaviarfit(
    model = "SAVX",
    tau,
    daily_ret = r_t_in_s_tin_xts,
    X = X_sav_rk_current,
    R = 1000,
    B = 1,
    Exp_Short = "Yes",
    std_err = "not_boot"
  ),
  error = function(e) NULL
)

############################################
# Check whether the estimation is valid
############################################

valid_fit_sav_rk <- (
  !is.null(fit_sav_es_rk) &&
  !is.null(fit_sav_es_rk$coef_mat) &&
  NCOL(fit_sav_es_rk$coef_mat) >= 1 &&
  all(is.finite(fit_sav_es_rk$coef_mat[, 1])) &&
  !is.null(fit_sav_es_rk$VaR) &&
  !is.null(fit_sav_es_rk$ES) &&
  all(is.finite(fit_sav_es_rk$VaR)) &&
  all(is.finite(fit_sav_es_rk$ES))
)

if (valid_fit_sav_rk) {

  ##########################################
  # Successful estimation
  ##########################################

  coef_sav_rk <-
    fit_sav_es_rk$coef_mat[, 1]

  VaR_SAVX_rk_in_s <-
    as.numeric(fit_sav_es_rk$VaR)

  ES_SAVX_rk_in_s <-
    as.numeric(fit_sav_es_rk$ES)

} else {

  ##########################################
  # Failed estimation:
  # use the coefficients from the previous
  # point and reconstruct VaR and ES on the
  # updated X window
  ##########################################

  if (
    is.null(coef_sav_rk_previous) ||
    is.null(VaR_SAVX_rk_in_s_previous)
  ) {
    stop(
      paste0(
        "SAV-X rk estimation failed at tt = ",
        tt,
        " and no previous valid estimates are available."
      )
    )
  }

  failed_SAVX_rk <-
    c(failed_SAVX_rk, tt)

  # Use the coefficients estimated at the previous point
  coef_sav_rk <-
    coef_sav_rk_previous

  n_sav_rk <- length(X_sav_rk_current)

  if (n_sav_rk != length(r_t_in_s_tin_xts)) {
    stop(
      paste0(
        "The lengths of X and daily_ret differ at tt = ",
        tt,
        "."
      )
    )
  }

  if (any(!is.finite(X_sav_rk_current))) {
    stop(
      paste0(
        "Non-finite values found in the SAV-X rk ",
        "in-sample window at tt = ",
        tt,
        "."
      )
    )
  }

  # Initialize the new in-sample VaR path
  VaR_SAVX_rk_in_s <-
    rep(NA_real_, n_sav_rk)

  # The first observation of the updated rolling window
  # corresponds to the second observation of the previous window
  VaR_SAVX_rk_in_s[1] <-
    VaR_SAVX_rk_in_s_previous[2]

  ##########################################
  # Reconstruct the in-sample VaR path
  ##########################################

  for (i in 2:n_sav_rk) {

    VaR_SAVX_rk_in_s[i] <- as.numeric(
      coef_sav_rk[1] +
      coef_sav_rk[2] *
        VaR_SAVX_rk_in_s[i - 1] +
      coef_sav_rk[3] *
        X_sav_rk_current[i - 1]
    )
  }

  ##########################################
  # Reconstruct the in-sample ES path
  ##########################################

  ES_SAVX_rk_in_s <- -abs(
    (1 + exp(coef_sav_rk[4])) *
      VaR_SAVX_rk_in_s
  )

  warning(
    paste0(
      "SAV-X rk estimation failed at tt = ",
      tt,
      ". Previous coefficients were used."
    ),
    call. = FALSE
  )
}

############################################
# One-step-ahead VaR and ES forecasts
############################################

fit_sav_es_VaR_oos_rk <- as.numeric(
  coef_sav_rk[1] +
  coef_sav_rk[2] *
    last(VaR_SAVX_rk_in_s) +
  coef_sav_rk[3] *
    last(X_sav_rk_current)
)

fit_sav_es_ES_oos_rk <- -abs(
  (1 + exp(coef_sav_rk[4])) *
    fit_sav_es_VaR_oos_rk
)

VaR_oos[tt, m] <-
  fit_sav_es_VaR_oos_rk

ES_oos[tt, m] <-
  fit_sav_es_ES_oos_rk

############################################
# Construct the rolling training data
############################################

if (tt == 1) {

  # First step: only in-sample observations
  VaR_training_data_mod[, m, tt] <-
    VaR_SAVX_rk_in_s

  ES_training_data_mod[, m, tt] <-
    ES_SAVX_rk_in_s

} else if (tt <= Tin) {

  # Mixed window: first in-sample, last out-of-sample
  n_in  <- Tin - tt + 1
  n_oos <- (n_in + 1):Tin

  VaR_training_data_mod[1:n_in, m, tt] <-
    VaR_SAVX_rk_in_s[1:n_in]

  VaR_training_data_mod[n_oos, m, tt] <-
    VaR_oos[1:(tt - 1), m]

  ES_training_data_mod[1:n_in, m, tt] <-
    ES_SAVX_rk_in_s[1:n_in]

  ES_training_data_mod[n_oos, m, tt] <-
    ES_oos[1:(tt - 1), m]

} else {

  # Only out-of-sample observations
  VaR_training_data_mod[, m, tt] <-
    VaR_oos[(tt - Tin):(tt - 1), m]

  ES_training_data_mod[, m, tt] <-
    ES_oos[(tt - Tin):(tt - 1), m]
}

############################################
# Store the current estimates for the
# following rolling estimation point
############################################

coef_sav_rk_previous <-
  coef_sav_rk

VaR_SAVX_rk_in_s_previous <-
  VaR_SAVX_rk_in_s

ES_SAVX_rk_in_s_previous <-
  ES_SAVX_rk_in_s

############################################ 
############################################ MF-QR-X (EPU+RVOL 5 min)
############################################ M30

m<-30

Q_linear_arch<-seq_test(r_t_est_cycle,tau)
Q_linear_arch<-round(Q_linear_arch[[1]],3)

q_hat<-which(Q_linear_arch>0.05)[1]

q_hat<-ifelse(q_hat<4|q_hat>9,4,q_hat)

repeat {
fit_mfx_es_rvol_5<-NULL

fit_mfx_es_rvol_5<-tryCatch(
uqfit(model="lARCHMIDASX", tau, daily_ret=r_t_est_cycle,q=q_hat,
Exp_Short="Yes",B=1,K=K,mv_m=epu_mv_cycle,X=rvol_5_est_cycle,
out_of_sample=lstep),
error=function(e) return(NULL)
)

if (!is.null(fit_mfx_es_rvol_5$coef_mat) && all(!is.na(fit_mfx_es_rvol_5$coef_mat[,1]))){
break
}

}

VaR_oos[tt,m]<-fit_mfx_es_rvol_5$VaR_oos
ES_oos[tt,m]<-fit_mfx_es_rvol_5$ES_oos

VaR_MFX_rvol_5_in_s<-fit_mfx_es_rvol_5$VaR
ES_MFX_rvol_5_in_s<-fit_mfx_es_rvol_5$ES


if (tt == 1) {
      # first step, only in-sample
      VaR_training_data_mod[, m, tt] <- VaR_MFX_rvol_5_in_s
      ES_training_data_mod[, m, tt] <- ES_MFX_rvol_5_in_s
      
    } else if (tt <= Tin) {
      # Mixed window: first in-sample, last oos
      n_in <- Tin - tt + 1
      n_oos <- (n_in+1):Tin
      
      VaR_training_data_mod[1:n_in, m, tt] <- VaR_MFX_rvol_5_in_s[1:n_in]
      VaR_training_data_mod[n_oos, m, tt] <- VaR_oos[1:(tt-1), m]
      
      ES_training_data_mod[1:n_in, m, tt] <- ES_MFX_rvol_5_in_s[1:n_in]
      ES_training_data_mod[n_oos, m, tt] <- ES_oos[1:(tt-1), m]
      
    } else {
      # Only oos 
      VaR_training_data_mod[, m, tt] <- VaR_oos[(tt - Tin):(tt-1), m]
      ES_training_data_mod[, m, tt] <- ES_oos[(tt - Tin):(tt-1), m]
    }
  
############################################ 
############################################ MF-QR-X (EPU+RB-SS)
############################################ M31

m<-31

repeat {
fit_mfx_es_rb_ss<-NULL

fit_mfx_es_rb_ss<-tryCatch(
uqfit(model="lARCHMIDASX", tau, daily_ret=r_t_est_cycle,q=q_hat,
Exp_Short="Yes",B=1,K=K,mv_m=epu_mv_cycle,X=rb_ss_est_cycle,
out_of_sample=lstep),
error=function(e) return(NULL)
)

if (!is.null(fit_mfx_es_rb_ss$coef_mat) && all(!is.na(fit_mfx_es_rb_ss$coef_mat[,1]))){
break
}

}

VaR_oos[tt,m]<-fit_mfx_es_rb_ss$VaR_oos
ES_oos[tt,m]<-fit_mfx_es_rb_ss$ES_oos

VaR_MFX_rb_ss_in_s<-fit_mfx_es_rb_ss$VaR
ES_MFX_rb_ss_in_s<-fit_mfx_es_rb_ss$ES


if (tt == 1) {
      # first step, only in-sample
      VaR_training_data_mod[, m, tt] <- VaR_MFX_rb_ss_in_s
      ES_training_data_mod[, m, tt] <- ES_MFX_rb_ss_in_s
      
    } else if (tt <= Tin) {
      # Mixed window: first in-sample, last oos
      n_in <- Tin - tt + 1
      n_oos <- (n_in+1):Tin
      
      VaR_training_data_mod[1:n_in, m, tt] <- VaR_MFX_rb_ss_in_s[1:n_in]
      VaR_training_data_mod[n_oos, m, tt] <- VaR_oos[1:(tt-1), m]
      
      ES_training_data_mod[1:n_in, m, tt] <- ES_MFX_rb_ss_in_s[1:n_in]
      ES_training_data_mod[n_oos, m, tt] <- ES_oos[1:(tt-1), m]
      
    } else {
      # Only oos 
      VaR_training_data_mod[, m, tt] <- VaR_oos[(tt - Tin):(tt-1), m]
      ES_training_data_mod[, m, tt] <- ES_oos[(tt - Tin):(tt-1), m]
    }
  



############################################ 
############################################ MF-QR-X (EPU+RK)
############################################ M32

m<-32

repeat {
fit_mfx_es_rk<-NULL

fit_mfx_es_rk<-tryCatch(
uqfit(model="lARCHMIDASX", tau, daily_ret=r_t_est_cycle,q=q_hat,
Exp_Short="Yes",B=1,K=K,mv_m=epu_mv_cycle,X=rk_est_cycle,
out_of_sample=lstep),
error=function(e) return(NULL)
)

if (!is.null(fit_mfx_es_rk$coef_mat) && all(!is.na(fit_mfx_es_rk$coef_mat[,1]))){
break
}

}

VaR_oos[tt,m]<-fit_mfx_es_rk$VaR_oos
ES_oos[tt,m]<-fit_mfx_es_rk$ES_oos

VaR_MFX_rk_in_s<-fit_mfx_es_rk$VaR
ES_MFX_rk_in_s<-fit_mfx_es_rk$ES

if (tt == 1) {
      # first step, only in-sample
      VaR_training_data_mod[, m, tt] <- VaR_MFX_rk_in_s
      ES_training_data_mod[, m, tt] <- ES_MFX_rk_in_s
      
    } else if (tt <= Tin) {
      # Mixed window: first in-sample, last oos
      n_in <- Tin - tt + 1
      n_oos <- (n_in+1):Tin
      
      VaR_training_data_mod[1:n_in, m, tt] <- VaR_MFX_rk_in_s[1:n_in]
      VaR_training_data_mod[n_oos, m, tt] <- VaR_oos[1:(tt-1), m]
      
      ES_training_data_mod[1:n_in, m, tt] <- ES_MFX_rk_in_s[1:n_in]
      ES_training_data_mod[n_oos, m, tt] <- ES_oos[1:(tt-1), m]
      
    } else {
      # Only oos 
      VaR_training_data_mod[, m, tt] <- VaR_oos[(tt - Tin):(tt-1), m]
      ES_training_data_mod[, m, tt] <- ES_oos[(tt - Tin):(tt-1), m]
    }
  

############################################
############################################
############################################ STORE DATA 
############################################
############################################


## 
if (tt %% 10 == 0) message("tt: ", tt, " out of ", nstep)


if (tt %% 50 == 0| tt == nstep){
save(
tt,
Tin,
tau,
list_of_models,
N_model,
nstep,
VaR_oos,
ES_oos,
r_t_oos_full,
r_t_in_s_matrix,
VaR_training_data_mod,
ES_training_data_mod,
file=filename)


}

} # 



############################################
############################################
############################################ Final storing data
############################################
############################################

r_t_oos<-coredata(r_t_oos_full)

save(
Tin,
tau,
list_of_models,
N_model,
nstep,
VaR_oos,
ES_oos,
r_t_oos_full,
r_t_oos,
r_t_in_s_matrix,
VaR_training_data_mod,
ES_training_data_mod,
file=filename)


