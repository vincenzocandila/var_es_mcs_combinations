
###########################################################################################################################
### This code replicates Table 2 of the paper:                                             						   		###
### 'Combining Value-at-Risk and Expected Shortfall Forecasts via the Model Confidence Set'						  		###
###########################################################################################################################

#### Analysis performed using R version 4.4.3

###################################
#### Load the libraries
###################################

library(xts)				# Version used: 0.14.1
library(fBasics)			# Version used: 4041.97
library(zoo)				# Version used: 1.8.13
library(DT)					# Version used: 0.33
library(htmltools)			# Version used: 0.5.8.1

###################################
#### Load the functions
###################################

source("functions/main_functions.R")

###################################
#### Set the period of interest
###################################

period_of_interest<-"2013-01/2022-06"

#####################################################################
#### Load SP500 (subset from Oxford-Man Institute's Realized Library)
#####################################################################

load("data/raw/sp500_data.RData")  # S&P500 data

######################################################################    
####	sp500_data includes:
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

#####################################################################
#### Load Shanghai Comp. (subset from Oxford-Man Institute's Realized Library)
#####################################################################

load("data/raw/shanghai_comp_data.RData") ## Shanghai Comp. data

######################################################################    
####	shanghai_comp_data includes:
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

epu_full<-read.csv(file="data/raw/Global_Policy_Uncertainty_Data.csv",sep=";",dec=",")
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
  title = "Table 2: Summary statistics"
)

