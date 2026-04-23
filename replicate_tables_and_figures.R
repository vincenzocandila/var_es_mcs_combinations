
#### Analysis performed using R version 4.4.3

########################################################################
######  REPLICATION CODE FOR		
######  'Combining Value-at-Risk and Expected Shortfall forecasts
######  via the Model Confidence Set'
########################################################################

###################################
#### Load the libraries
###################################

library(xts)
library(fBasics)
library(rugarch)
library(np)
library(esback)
library(GAS)
library(zoo)
library(DT)
library(htmltools)

###################################
#### Load the functions
###################################

source("functions/main_functions.R")

###################################
#### Set the period of interest
###################################

period_of_interest<-"2013-01/2022-06"

#####################################################################
#### Load data (subset from Oxford-Man Institute's Realized Library)
#####################################################################

load("data/SP500_data.RData") 

######################################################################    
####	SP500_data includes:
## 	r_t	  	: daily close-to-close log returns; 
## 	rb_ss 	: realized volatility computed using realized bipower variation with subsampling;
## 	rk		: realized kernel volatility;
## 	rvol5	: realized volatility based on 5-minute intervals.
######################################################################    

sp500_summary<-rbind(
sum_f(r_t[period_of_interest]),
sum_f(rvol5[period_of_interest]),
sum_f(rb_ss[period_of_interest]),
sum_f(rk[period_of_interest]))

load("data/Shanghai_comp_data.RData") ## Shanghai Comp. data

######################################################################    
####	Shanghai_comp_data includes:
## 	r_t	  	: daily close-to-close log returns; 
## 	rb_ss 	: realized volatility computed using realized bipower variation with subsampling;
## 	rk		: realized kernel volatility;
## 	rvol5	: realized volatility based on 5-minute intervals.
###################################################################### 

shanghai_summary<-rbind(
sum_f(r_t[period_of_interest]),
sum_f(rvol5[period_of_interest]),
sum_f(rb_ss[period_of_interest]),
sum_f(rk[period_of_interest]))

###### Economic Policy Uncertainty data collected from:
###### https://www.policyuncertainty.com/global_monthly.html
###### We use the column labeled "GEPU_current"
###### The EPU indices were downloaded on May 26, 2023. 
###### Note that the series may be subject to revisions over time

epu_full<-read.csv(file="data/Global_Policy_Uncertainty_Data.csv",sep=";",dec=",")
epu_full_i <- make_xts(epu_full, value_col = names(epu_full)[3])

MV<-diff(epu_full_i)/lag(epu_full_i)
MV[1]<-0
MV_row<-sum_f(MV['2010/2022-06'])

###################################
#### Replication of Table 2
###################################

tab_summary<-rbind(
sp500_summary,
shanghai_summary,
MV_row)

colnames(tab_summary)<-c("Obs.","Min.","Max.","Mean","SD","Skew.","Kurt.")

summary_tab_dt_f(
  tab_summary,
  title = "Table 2. Summary statistics"
)

##########################################
#### Load data for plots 
##########################################

load("data/res_shanghai_comp_tau_0.025.RData")

#############################################################################
#### res_shanghai_comp_tau_0.025 includes:
## r_t_oos_full_plot				: daily log-returns for the training plots;
## r_t_oos_full					: daily log-returns for the out-of-sample period;
## nstep							: total number of model re-estimations (i.e., length of the out-of-sample period);
## lab_full						: model labels, including 32 models, 4 benchmarks (EW-Comb, Median-Comb, RS-Comb, MS-Comb),
##           							and 6 proposed combined predictors (MCS-Comb, WL-MCS-Comb, MCS-RS-Comb, WL-MCS-RS-Comb, 
##           							MCS-MS-Comb, WL-MCS-MS-Comb);
## tau							: coverage level;
## Backtesting_in_sample_pvalues	: p-values of the six backtests (see Table 4 in the paper), by model (32) and over nstep;
## db_loss_oos					: matrix of FZ losses for the out-of-sample period for the 42 models and combined predictors;
## VaR_oos_ev and ES_oos_ev		: matrices of VaR and ES for the out-of-sample period
##                           			for the 42 models and combined predictors;
## MCS_est_included				: list of models included at each step (nstep) in the Set of Superior Models (SSM)
##                  					from the Model Confidence Set (MCS), using the unweighted FZ loss;
## MCS_est_included_w				: list of models included at each step (nstep) in the Set of Superior Models (SSM)
##                     				from the Model Confidence Set (MCS), using the weighted FZ loss.
#############################################################################

###################################
#### Replication of Figure 2(a)
###################################

plot_inclusion(
  X = Backtesting_in_sample_pvalues,
  nstep = nstep,
  lab_full = lab_full,
  r_t_plot = r_t_oos_full_plot,
  N_model = 32,
  square_col="green",
  input_type = "Backtesting"
)

###################################
#### Replication of Figure 2(b)
###################################

plot_inclusion(
  X = MCS_est_included,
  nstep = nstep,
  lab_full = lab_full,
  r_t_plot = r_t_oos_full_plot,
  N_model = 32,
  square_col="blue",
  input_type = "MCS"
)

###################################
#### Replication of Figure 2(c)
###################################

plot_inclusion(
  X = MCS_est_included_w,
  nstep = nstep,
  lab_full = lab_full,
  r_t_plot = r_t_oos_full_plot,
  N_model = 32,
  square_col="black",
  input_type = "MCS"
)

###################################
#### Code for replication of Table 5 
###################################

load("data/res_sp500_tau_0.025.RData")

#############################################################################
#### res_sp500_tau_0.025 includes:
## r_t_oos_full_plot				: daily log-returns for the training plots;
## r_t_oos_full					: daily log-returns for the out-of-sample period;
## nstep							: total number of model re-estimations (i.e., length of the out-of-sample period);
## lab_full						: model labels, including 32 models, 4 benchmarks (EW-Comb, Median-Comb, RS-Comb, MS-Comb),
##           							and 6 proposed combined predictors (MCS-Comb, WL-MCS-Comb, MCS-RS-Comb, WL-MCS-RS-Comb, 
##           							MCS-MS-Comb, WL-MCS-MS-Comb);
## tau							: coverage level;
## Backtesting_in_sample_pvalues	: p-values of the six backtests (see Table 4 in the paper), by model (32) and over nstep;
## db_loss_oos					: matrix of FZ losses for the out-of-sample period for the 42 models and combined predictors;
## VaR_oos_ev and ES_oos_ev		: matrices of VaR and ES for the out-of-sample period
##                           			for the 42 models and combined predictors;
## MCS_est_included				: list of models included at each step (nstep) in the Set of Superior Models (SSM)
##                  					from the Model Confidence Set (MCS), using the unweighted FZ loss;
## MCS_est_included_w				: list of models included at each step (nstep) in the Set of Superior Models (SSM)
##                     				from the Model Confidence Set (MCS), using the weighted FZ loss.
#############################################################################

#### Setting of the MCS

set.seed(123)

B <- 1000                  		# number of bootstrap replicates for the MCS
alpha <- 0.25              		# significance level for the MCS
N_model_full <- length(lab_full)   # total number of models and combined predictors

#### Store the averages

col_means<-colMeans(db_loss_oos)

#### Evaluation MCS

MCS_est_with_c<-MCS_f(db_loss_oos,B=B,alpha=alpha)

col_mcs<-cbind(
round(col_means,3),
ifelse(1:N_model_full %in% MCS_est_with_c$includedR,1,0)
)

colnames(col_mcs)<-c("FZLoss","SSM")

#### VaR backtesting

cols_var<-matrix(rep(NA,3*N_model_full),ncol=3)

for (i in 1:N_model_full){
cols_var[i,1]<-as.numeric(BacktestVaR(r_t_oos_full, VaR_oos_ev[,i], tau)$LRuc[2])
cols_var[i,2]<-as.numeric(BacktestVaR(r_t_oos_full, VaR_oos_ev[,i], tau)$LRcc[2])
cols_var[i,3]<-as.numeric(BacktestVaR(r_t_oos_full, VaR_oos_ev[,i], tau, Lags = 2)$DQ[2])
}

cols_var<-round(cols_var,3)

colnames(cols_var)<-c("UC","CC","DQ")

#### ES backtesting

cols_es<-matrix(rep(NA,3*N_model_full),ncol=3)

r_t_oos<-zoo::coredata(r_t_oos_full)

for (i in 1:N_model_full){
# Bayer & Dimitriadis (2020) test version 1
cols_es[i,1]<-as.numeric(esr_backtest(r_t_oos, VaR_oos_ev[,i], ES_oos_ev[,i],tau, version = 1)$pvalue_twosided_asymptotic)
# Bayer & Dimitriadis (2020) test version 2
cols_es[i,2]<-as.numeric(esr_backtest(r_t_oos, VaR_oos_ev[,i], ES_oos_ev[,i],tau, version = 2)$pvalue_twosided_asymptotic)
# Bayer & Dimitriadis (2020) test version 3
cols_es[i,3]<-as.numeric(esr_backtest(r_t_oos, VaR_oos_ev[,i], ES_oos_ev[,i],tau, version = 3)$pvalue_twosided_asymptotic)
message(i)
}

colnames(cols_es)<-c("BD-1","BD-2","BD-3")
cols_es<-round(cols_es,3)

final_tab<-cbind(cols_var,cols_es,col_mcs)

### 1 if all backtests are passed (p-values >= 0.05), 0 otherwise
Backtest_all <- as.integer(rowSums(final_tab[, 1:6] >= 0.05) == 6)

final_tab<-cbind(final_tab,Backtest=Backtest_all)

######################################
#### Replication of Table 5
######################################

# To reproduce exactly the p-values of the Bayer & Dimitriadis (2020) tests reported in the paper, 
# the code should be run under R version 4.4.3.
# When using different R versions or computing environments, small numerical differences may arise;
# these affect only the p-values of these tests and do not alter the main conclusions.

tab_dt_f(final_tab, title = "Table 5: S&P 500 out-of-sample evaluation")

######################################
#### Code for replication of Table 6 
######################################

## To ensure full reproducibility of the esr_backtest results, 
## each index should be evaluated in a separate fresh R session after setting the random seed.

load("data/res_shanghai_comp_tau_0.025.RData")

#############################################################################
#### res_shanghai_comp_tau_0.025 includes:
## r_t_oos_full_plot				: daily log-returns for the training plots;
## r_t_oos_full					: daily log-returns for the out-of-sample period;
## nstep							: total number of model re-estimations (i.e., length of the out-of-sample period);
## lab_full						: model labels, including 32 models, 4 benchmarks (EW-Comb, Median-Comb, RS-Comb, MS-Comb),
##           							and 6 proposed combined predictors (MCS-Comb, WL-MCS-Comb, MCS-RS-Comb, WL-MCS-RS-Comb, 
##           							MCS-MS-Comb, WL-MCS-MS-Comb);
## tau							: coverage level;
## Backtesting_in_sample_pvalues	: p-values of the six backtests (see Table 4 in the paper), by model (32) and over nstep;
## db_loss_oos					: matrix of FZ losses for the out-of-sample period for the 42 models and combined predictors;
## VaR_oos_ev and ES_oos_ev		: matrices of VaR and ES for the out-of-sample period
##                           			for the 42 models and combined predictors;
## MCS_est_included				: list of models included at each step (nstep) in the Set of Superior Models (SSM)
##                  					from the Model Confidence Set (MCS), using the unweighted FZ loss;
## MCS_est_included_w				: list of models included at each step (nstep) in the Set of Superior Models (SSM)
##                     				from the Model Confidence Set (MCS), using the weighted FZ loss.
#############################################################################


#### Setting of the MCS

set.seed(123)
B <- 1000                  		# number of bootstrap replicates for the MCS
alpha <- 0.25              		# significance level for the MCS
N_model_full <- length(lab_full)   # total number of models and combined predictors

#### Store the averages

col_means<-colMeans(db_loss_oos)

#### Evaluation MCS

MCS_est_with_c<-MCS_f(db_loss_oos,B=B,alpha=alpha)

col_mcs<-cbind(
round(col_means,3),
ifelse(1:N_model_full %in% MCS_est_with_c$includedR,1,0)
)

colnames(col_mcs)<-c("FZLoss","SSM")

#### VaR backtesting

cols_var<-matrix(rep(NA,3*N_model_full),ncol=3)

for (i in 1:N_model_full){
cols_var[i,1]<-as.numeric(BacktestVaR(r_t_oos_full, VaR_oos_ev[,i], tau)$LRuc[2])
cols_var[i,2]<-as.numeric(BacktestVaR(r_t_oos_full, VaR_oos_ev[,i], tau)$LRcc[2])
cols_var[i,3]<-as.numeric(BacktestVaR(r_t_oos_full, VaR_oos_ev[,i], tau, Lags = 2)$DQ[2])
}

cols_var<-round(cols_var,3)

colnames(cols_var)<-c("UC","CC","DQ")

#### ES backtesting

cols_es<-matrix(rep(NA,3*N_model_full),ncol=3)

r_t_oos<-zoo::coredata(r_t_oos_full)

for (i in 1:N_model_full){
# Bayer & Dimitriadis (2020) test version 1
cols_es[i,1]<-as.numeric(esr_backtest(r_t_oos, VaR_oos_ev[,i], ES_oos_ev[,i],tau, version = 1)$pvalue_twosided_asymptotic)
# Bayer & Dimitriadis (2020) test version 2
cols_es[i,2]<-as.numeric(esr_backtest(r_t_oos, VaR_oos_ev[,i], ES_oos_ev[,i],tau, version = 2)$pvalue_twosided_asymptotic)
# Bayer & Dimitriadis (2020) test version 3
cols_es[i,3]<-as.numeric(esr_backtest(r_t_oos, VaR_oos_ev[,i], ES_oos_ev[,i],tau, version = 3)$pvalue_twosided_asymptotic)
message(i)
}

colnames(cols_es)<-c("BD-1","BD-2","BD-3")
cols_es<-round(cols_es,3)

final_tab<-cbind(cols_var,cols_es,col_mcs)

### 1 if all backtests are passed (p-values >= 0.05), 0 otherwise
Backtest_all <- as.integer(rowSums(final_tab[, 1:6] >= 0.05) == 6)

final_tab<-cbind(final_tab,Backtest=Backtest_all)

###################################
#### Replication of Table 6
###################################

# As for Table 5, to reproduce exactly the p-values of the Bayer & Dimitriadis (2020) tests reported in the paper, 
# the code should be run under R version 4.4.3.
# When using different R versions or computing environments, small numerical differences may arise;
# these affect only the p-values of these tests and do not alter the main conclusions.

tab_dt_f(final_tab, title = "Table 6: Shanghai Composite out-of-sample evaluation")

###################################
#### Replication of Table 7
###################################

load("data/all_indices_tau_0.025.RData")

#######################################################################################
#### all_indices_tau_0.025 includes:
## list_of_tabs	: list of tables (one per index) with models in rows and
##               		backtesting p-values (UC, CC, DQ, BD-1, BD-2, BD-3) and FZLoss in columns;
##               		LaTeX markers (e.g., cellcolor{...}) indicate whether models pass all backtests at the 5% level
##               		and belong to the SSM of the MCS at the 25% significance level; 
## lab_full		: model labels, including 32 models, 4 benchmarks (EW-Comb, Median-Comb, RS-Comb, MS-Comb),
##           			and 6 proposed combined predictors (MCS-Comb, WL-MCS-Comb, MCS-RS-Comb, WL-MCS-RS-Comb, 
##          			MCS-MS-Comb, WL-MCS-MS-Comb).
## VaR and ES are estimated at the 2.5% level.
#######################################################################################

res_final_0.025 <- build_summary_from_list(
  tabs_list = list_of_tabs,
  lab_full = lab_full,
  as_percent = FALSE
)

rm(list_of_tabs)

load("data/all_indices_tau_0.01.RData")

#######################################################################################
#### all_indices_tau_0.01 includes:
## list_of_tabs	: list of tables (one per index) with models in rows and
##               		backtesting p-values (UC, CC, DQ, BD-1, BD-2, BD-3) and FZLoss in columns;
##               		LaTeX markers (e.g., cellcolor{...}) indicate whether models pass all backtests at the 5% level
##               		and belong to the SSM of the MCS at the 25% significance level; 
## lab_full		: model labels, including 32 models, 4 benchmarks (EW-Comb, Median-Comb, RS-Comb, MS-Comb),
##           			and 6 proposed combined predictors (MCS-Comb, WL-MCS-Comb, MCS-RS-Comb, WL-MCS-RS-Comb, 
##          			MCS-MS-Comb, WL-MCS-MS-Comb).
## VaR and ES are estimated at the 1% level.
#######################################################################################

res_final_0.01 <- build_summary_from_list(
  tabs_list = list_of_tabs,
  lab_full = lab_full,
  as_percent = FALSE
)

tab_7<-cbind(res_final_0.025,res_final_0.01)
tab_7

### Table 7 in hmtl
tab_7_dt_f(tab_7)


