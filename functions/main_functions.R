
################################## Summary function

###############################################################################################
# Computes the summary statistics reported in Table 2
# Input  : numeric vector of the variable of interest
# Output : matrix containing summary statistics (N, min, max, mean, sd, skewness, kurtosis)
###############################################################################################

sum_f<-function(x){

cbind(
length(x),
min(x),
max(x),
mean(x),
sd(x),
fBasics::skewness(x),
fBasics::kurtosis(x))
}

################################ EPU importing function

###############################################################################################
# Converts a data frame (downloaded from https://www.policyuncertainty.com/global_monthly.html) 
#  with Year and Month columns into an xts time series
# Input  : db         -> data frame containing the data
#          value_col  -> name of the column to be converted into a time series
#          year_col   -> name of the year column (default: "Year")
#          month_col  -> name of the month column (default: "Month")
# Output : xts object indexed by monthly dates
###############################################################################################


make_xts <- function(db, value_col, year_col = "Year", month_col = "Month") {
  
dates <- as.Date(sprintf("%04d-%02d-01", db[[year_col]], db[[month_col]]))
xts::xts(db[[value_col]], order.by = dates)

}


################################## Exponential smoothing

###############################################################################################
# Computes the exponentially weighted loss as in Eq. (7) of the paper
# Input  : theta -> smoothing parameter (lambda in Eq. (7))
#          loss  -> matrix (T x M) of FZ losses, where T is the time dimension and M the number of models
# Output : matrix (T x M) of exponentially weighted FZ losses
###############################################################################################

ewma_v <- function(theta,loss) {
TT<-nrow(loss)
M<-ncol(loss)
wloss<-matrix(0,ncol=M,nrow=TT)

for (m in 1:M){
for (tt in 2:(TT)){
wloss[tt,m]<-theta*loss[tt,m]+(1-theta)*wloss[tt-1,m]
}  
}
return(wloss)
}

################################## Fix zeros 

###############################################################################################
# Replaces zero values in a numeric vector using the nearest non-zero neighbours
# Input  : x -> numeric vector possibly containing zeros
# Output : numeric vector where zeros are replaced as follows:
#          - if non-zero values exist on both sides, the zero is replaced by their average
#          - if only a left non-zero value exists, the zero is replaced by that value
#          - if only a right non-zero value exists, the zero is replaced by that value
###############################################################################################


fix_zeros <- function(x) {
  nz <- which(x != 0)
  if (length(nz) == length(x)) return(x)  # nessuno zero, ritorna invariato
  for (i in which(x == 0)) {
    left  <- if (any(nz < i)) max(nz[nz < i]) else NA
    right <- if (any(nz > i)) min(nz[nz > i]) else NA
    if (!is.na(left) && !is.na(right)) {
      x[i] <- mean(c(x[left], x[right]))
    } else if (!is.na(left)) {
      x[i] <- x[left]
    } else if (!is.na(right)) {
      x[i] <- x[right]
    }
  }
  x
}

##################################### Cornish-Fisher VaR and ES

###############################################################################################
# Computes standardized Cornish-Fisher VaR and ES estimates
# Input  : ret -> vector of returns
#          ht  -> vector of conditional standard deviations (volatility forecasts)
#          tau -> probability level for VaR/ES, e.g. 0.01 or 0.025
# Output : matrix containing:
#          VaR_cf -> Cornish-Fisher adjusted quantile
#          ES_cf  -> Cornish-Fisher adjusted Expected Shortfall
# Notes  : the procedure computes skewness and excess kurtosis of standardized
#          residuals and applies the Cornish-Fisher expansion to correct the
#          Gaussian quantile. ES is obtained using the adjustment proposed in
#          Maillard's "A User's Guide to the Cornish Fisher Expansion".
###############################################################################################

Quantile_CF_f<-function(ret,ht,tau){

residuals_std <- ret / ht

# Step 1: skewness and excess kurtosis
skew <- mean((residuals_std - mean(residuals_std))^3) / sd(residuals_std)^3
kurt <- mean((residuals_std - mean(residuals_std))^4) / sd(residuals_std)^4 - 3

# Step 2: normal quantile 
z_norm <- qnorm(tau)

# Step 3: CF-approx quantile
z_CF <- z_norm +
  (1/6)*(z_norm^2 - 1)*skew +
  (1/24)*(z_norm^3 - 3*z_norm)*kurt -
  (1/36)*(2*z_norm^3 - 5*z_norm)*skew^2

# Step 4: correction from "A User’s Guide to the Cornish Fisher Expansion" by Didier MAILLARD

Wa <- 1 + z_norm*(skew/6) + (1-2*z_norm^2)*(skew^2/36)+(-1+z_norm^2)*(kurt/24) 

# Step 5: quantile to be returned

VaR_cf <-  z_CF
ES_cf <-  (-dnorm(z_norm) / tau * Wa) 

res<-cbind(VaR_cf=VaR_cf,ES_cf=ES_cf)

return(res)

}



################################## Historical simulation 

###############################################################################################
# Computes VaR and ES forecasts using the Historical Simulation approach
# Input  : rt           -> xts object of returns
#          w            -> rolling window size used to compute the empirical quantile
#          tau          -> probability level for VaR, e.g. 0.01 or 0.025
#          begin_period -> index identifying the first out-of-sample observation
#          end_period   -> index identifying the last out-of-sample observation
# Output : list containing:
#          VaR_HS 	  -> xts object of Historical Simulation VaR forecasts
#          ES_HS        -> xts object of Historical Simulation ES forecasts
###############################################################################################

hs_f<-function(rt,w,tau,begin_period,end_period){

myDate_b <- as.Date(time(first(rt[begin_period])))
myDate_e <- as.Date(time(last(rt[end_period])))

beg_sample<-which(as.Date(time(rt))==myDate_b) - w # new begin

beg_sample_f<-which(as.Date(time(rt))==myDate_b)

end_sample<-which(as.Date(time(rt))==myDate_e)   # end


r_t_hs<-coredata(rt[beg_sample:end_sample])

VaR_HS<-list()
ES_HS<-list()

for(tt in (w+1):length(r_t_hs)){

ret<-r_t_hs[(tt-w):(tt-1)]

VaR_HS[[tt]]<-quantile(ret,tau)
Ind<-ifelse(ret<=rep(VaR_HS[[tt]],w),1,0)

ES_HS[[tt]]<-sum(ret*Ind)/sum(Ind)

}

VaR_HS<-unlist(VaR_HS)
VaR_HS<-as.xts(VaR_HS,time(rt[beg_sample_f:end_sample]))

ES_HS<-unlist(ES_HS)
ES_HS<-as.xts(ES_HS,time(rt[beg_sample_f:end_sample]))

### final results


res<-list(
VaR_HS=VaR_HS,
ES_HS=ES_HS)

return(res)

}


################################## Riskmetrics

###############################################################################################
# Computes conditional variances using the RiskMetrics EWMA model
# Input  : ret    -> vector or xts object of returns
#          lambda -> RiskMetrics decay factor (typically 0.94 for daily data)
# Output : ht     -> vector of conditional variance estimates
###############################################################################################

riskmetrics_f<-function(ret,lambda){

ret<-coredata(ret)
N<-length(ret)
ht<-rep(var(ret),N)

for(i in 2:N){

ht[i]<-lambda*ht[i-1]+(1-lambda)*ret[i-1]^2

}

return(ht)
}


################################ Plot Backtesting and MCS inclusion

###############################################################################################
# Plots model inclusion over time as in Figures 2(a), 2(b), or 2(c) of the paper
# Main input 		: X -> object containing inclusion information
#                  - list of selected models over time if input_type = "MCS"
#                  - array of backtesting p-values if input_type = "Backtesting"
# Main output		: inclusion plot;
###############################################################################################

plot_inclusion <- function(X,
                           nstep,
                           lab_full,
                           r_t_plot,
                           N_model = 32,
                           input_type = c("MCS", "Backtesting"),
                           main = "",
                           inclusion_threshold = 0.95,
                           red_col = "red",
                           square_col = "black",
                           grid_col = "grey") {
  
  input_type <- match.arg(input_type)
  
  point_plot_1 <- matrix(NA, ncol = nstep, nrow = N_model)
  point_plot_2 <- matrix(NA, ncol = nstep, nrow = N_model)
  
  rownames(point_plot_1) <- rownames(point_plot_2) <- lab_full[1:N_model]
  
  for (tt in 1:nstep) {
    
    if (input_type == "MCS") {
      point_plot_1[, tt] <- ifelse(1:N_model %in% X[[tt]], 1, 0)
    }
    
    if (input_type == "Backtesting") {
      step1 <- rowSums(ifelse(X[, , tt] > 0.05, 1, 0))
      point_plot_1[, tt] <- ifelse(step1 == 6, 1, 0)
    }
    
    point_plot_2[, tt] <- ifelse(point_plot_1[, tt] == 1, N_model:1, -20)
  }
  
  lab <- lab_full[1:N_model]
  if (N_model >= 26) {
    lab[24] <- "CAViaR-SAV"
    lab[25] <- "CAViaR-AS"
    lab[26] <- "CAViaR-IG"
  }
  
  Days <- time(r_t_plot)
  Months <- format(as.Date(Days), "%m/%Y")
  first_in_month <- !duplicated(Months)
  
  ticks  <- which(first_in_month)
  labels <- Months[first_in_month]
  
  if (length(ticks) > 1) {
    ticks  <- ticks[-1]
    labels <- labels[-1]
  }
  
  models_over_threshold <- which(rowSums(point_plot_1) > floor(nstep * inclusion_threshold))
  
  par(mar = c(5.1, 11.1, 4.1, 3.1), mgp = c(3, 0.6, 0), xpd = NA)
  
  plot(x = NULL,
       y = NULL,
       xlim = c(1, nstep),
       ylim = c(1, N_model),
       xlab = "",
       ylab = "",
       xaxt = "n",
       yaxt = "n",
       main = main)
  
  axis(2, at = 1:N_model, labels = rev(lab), las = 2, cex.axis = 1.1)
  axis(1, at = ticks, labels = labels, las = 3, cex.axis = 0.9)
  
  segments(x0 = 1, y0 = 1:N_model, x1 = nstep, y1 = 1:N_model, col = grid_col)
  
  for (i in 1:nstep) {
    points(rep(i, N_model), rev(point_plot_2[, i]), pch = 15, cex = 1, col = square_col)
  }
  
  for (i in 1:nstep) {
    for (m in models_over_threshold) {
      y_val <- N_model - m + 1
      if (point_plot_2[m, i] == -20) {
        points(i, y_val, pch = 19, cex = 1, col = red_col)
      }
    }
  }
  
  invisible(list(
    point_plot_1 = point_plot_1,
    point_plot_2 = point_plot_2,
    models_over_threshold = models_over_threshold
  ))
}

####################################### MCS function

MCS_f<-function(LOSS,B,alpha){

###############################################################################################
# Performs the Model Confidence Set (MCS) procedure through rugarch 
# Input  : LOSS  -> matrix (T x M) of loss values, where T is the time dimension and M the number of models
#          B     -> number of bootstrap replications
#          alpha -> significance level of the MCS procedure
# Output : object returned by rugarch::mcsTest containing the selected model set and test statistics
###############################################################################################

set.seed(123)

	k_len			<-	max(np::b.star(LOSS)[,2])  # automatic block length
					# as in Patton, Politis and White (2009)

	k_len			<- ifelse(k_len<1.99,2,k_len)
	
	temp 		<-	rugarch::mcsTest(LOSS, 
						alpha = alpha, nboot = B, nblock = k_len,
  						boot = c("block"))

return(temp)
}

####################################### First html function for summary table (Table 2)

###############################################################################################
# Builds the data frame used to create the HTML version of Table 2
# Input  : tab_summary      -> matrix or data frame of summary statistics
#          sp500_period     -> sample period reported for the S&P500 series
#          shanghai_period  -> sample period reported for the Shanghai Composite series
#          epu_period       -> sample period reported for the EPU series
# Output : data frame formatted for HTML display, including group header rows
###############################################################################################

summary_tab_html_f <- function(tab_summary,
                               sp500_period = "from 2013-01-02 to 2022-06-28",
                               shanghai_period = "from 2013-01-04 to 2022-06-28",
                               epu_period = "from 01-01-2010 to 01-06-2022") {
  
  #### Convert to data frame
  table_data <- as.data.frame(tab_summary, stringsAsFactors = FALSE)
  
  #### Add the series names
  table_data <- cbind(
    Series = c("log-returns", "RVOL5", "RB-SS", "RK",
               "log-returns", "RVOL5", "RB-SS", "RK",
               "\u0394EPU"),
    table_data,
    stringsAsFactors = FALSE
  )
  
  #### Assign column names
  colnames(table_data) <- c("Series", "Obs.", "Min.", "Max.", "Mean", "SD", "Skew.", "Kurt.")
  
  #### Format numeric columns
  table_data[, "Obs."] <- formatC(as.numeric(table_data[, "Obs."]), format = "f", digits = 0)
  
  for (j in c("Min.", "Max.", "Mean", "SD", "Skew.", "Kurt.")) {
    table_data[, j] <- formatC(as.numeric(table_data[, j]), format = "f", digits = 3)
  }
  
  #### Build group header rows
  group_row_1 <- data.frame(
    Series = paste0("S&P500; Sample period: ", sp500_period),
    "Obs." = "", "Min." = "", "Max." = "", "Mean" = "",
    "SD" = "", "Skew." = "", "Kurt." = "",
    stringsAsFactors = FALSE
  )
  
  group_row_2 <- data.frame(
    Series = paste0("Shanghai Comp.; Sample period: ", shanghai_period),
    "Obs." = "", "Min." = "", "Max." = "", "Mean" = "",
    "SD" = "", "Skew." = "", "Kurt." = "",
    stringsAsFactors = FALSE
  )
  
  group_row_3 <- data.frame(
    Series = paste0("EPU; Sample period: ", epu_period),
    "Obs." = "", "Min." = "", "Max." = "", "Mean" = "",
    "SD" = "", "Skew." = "", "Kurt." = "",
    stringsAsFactors = FALSE
  )
  
  #### Combine everything
  table_html <- rbind(
    group_row_1,
    table_data[1:4, ],
    group_row_2,
    table_data[5:8, ],
    group_row_3,
    table_data[9, , drop = FALSE]
  )
  
  rownames(table_html) <- NULL
  
  return(table_html)
}


####################################### Second html function for summary table (Table 2)

###############################################################################################
# Creates the DT::datatable version of Table 2
# Input  : tab_summary      -> matrix or data frame of summary statistics
#          title            -> table title
#          sp500_period     -> sample period reported for the S&P500 series
#          shanghai_period  -> sample period reported for the Shanghai Composite series
#          epu_period       -> sample period reported for the EPU series
# Output : DT::datatable object reproducing Table 2 in HTML format
###############################################################################################


summary_tab_dt_f <- function(tab_summary,
                             title = "Table 2. Summary statistics",
                             sp500_period = "from 2013-01-02 to 2022-06-28",
                             shanghai_period = "from 2013-01-04 to 2022-06-28",
                             epu_period = "from 01-01-2010 to 01-06-2022") {
  
  table_html <- summary_tab_html_f(
    tab_summary = tab_summary,
    sp500_period = sp500_period,
    shanghai_period = shanghai_period,
    epu_period = epu_period
  )
  
  datatable(
    table_html,
    rownames = FALSE,
    escape = FALSE,
    caption = if (!is.null(title)) {
      htmltools::tags$caption(
        style = "caption-side: top; text-align: center; font-weight: bold; font-size: 16px;",
        title
      )
    },
    options = list(
      dom = "t",
      paging = FALSE,
      ordering = FALSE,
      autoWidth = TRUE,
      rowCallback = JS(
        "function(row, data, index) {",
        "  if (index == 0 || index == 5 || index == 10) {",
        "    $('td', row).css({'font-style': 'italic',",
        "                       'border-top': '1px solid black',",
        "                       'border-bottom': '1px solid black'});",
        "    $('td:gt(0)', row).html('');",
        "  }",
        "}"
      )
    ),
    class = "compact"
  )
}

####################################### First html function for evaluation tables (Tables 5 and 6)

###############################################################################################
# Builds the HTML-formatted version of Tables 5 and 6
# Input  : x -> matrix or data frame containing backtesting results, FZLoss, and auxiliary indicators (for the colors)
# Output : data frame formatted for HTML display, with conditional shading applied to selected cells
###############################################################################################

tab_html_f <- function(x) {
  
  #### Convert to data frame
  table_data <- as.data.frame(x, stringsAsFactors = FALSE)
  
  #### Add model names as first column
  table_data <- cbind(Model = rownames(table_data), table_data, stringsAsFactors = FALSE)
  
  #### Assign column names
  colnames(table_data) <- c("Model", "UC", "CC", "DQ", "BD-1", "BD-2", "BD-3",
                            "FZLoss", "SSM", "Backtest")
  
  #### Helper function to wrap a cell in HTML with optional background color
  make_cell <- function(value, bg = NULL, digits = NULL) {
    
    numeric_value <- suppressWarnings(as.numeric(value))
    
    if (!is.na(numeric_value) && !is.null(digits)) {
      formatted_value <- formatC(numeric_value, format = "f", digits = digits)
    } else {
      formatted_value <- as.character(value)
    }
    
    if (is.null(bg)) {
      return(formatted_value)
    }
    
    paste0(
      '<div style="background-color:', bg,
      '; padding:2px 4px; box-sizing:border-box; width:100%;">',
      formatted_value,
      '</div>'
    )
  }
  
  #### Create a copy to be styled in HTML
  table_html <- table_data
  
  #### Shade model names in light gray when Backtest = 1
  for (i in 1:nrow(table_data)) {
    backtest_flag <- suppressWarnings(as.numeric(table_data[i, "Backtest"]))
    
    if (!is.na(backtest_flag) && backtest_flag == 1) {
      table_html[i, "Model"] <- make_cell(table_data[i, "Model"], bg = "#E6E6E6")
    } else {
      table_html[i, "Model"] <- as.character(table_data[i, "Model"])
    }
  }
  
  #### Format the six backtesting columns without additional shading
  for (j in c("UC", "CC", "DQ", "BD-1", "BD-2", "BD-3")) {
    table_html[, j] <- sapply(table_data[, j], function(x) make_cell(x, digits = 3))
  }
  
  #### Shade FZLoss in dark gray when SSM = 1
  for (i in 1:nrow(table_data)) {
    ssm_flag <- suppressWarnings(as.numeric(table_data[i, "SSM"]))
    
    if (!is.na(ssm_flag) && ssm_flag == 1) {
      table_html[i, "FZLoss"] <- make_cell(table_data[i, "FZLoss"], bg = "#A6A6A6", digits = 3)
    } else {
      table_html[i, "FZLoss"] <- make_cell(table_data[i, "FZLoss"], digits = 3)
    }
  }
  
  #### Remove the auxiliary columns from the final HTML table
  table_html <- table_html[, !(colnames(table_html) %in% c("SSM", "Backtest")), drop = FALSE]
  
  rownames(table_html) <- NULL
  
  return(table_html)
}


####################################### Second html function for evaluation tables (Tables 5 and 6)

###############################################################################################
# Creates the DT::datatable version of Tables 5 and 6
# Input  : x     -> matrix or data frame containing backtesting results, FZLoss, and auxiliary indicators
#          title -> optional table title
# Output : DT::datatable object reproducing the backtesting table in HTML format
###############################################################################################

tab_dt_f <- function(x, title = NULL) {
  
  table_html <- tab_html_f(x)
  
  datatable(
    table_html,
    rownames = FALSE,
    escape = FALSE,
    caption = if (!is.null(title)) {
      htmltools::tags$caption(
        style = "caption-side: top; text-align: center; font-weight: bold; font-size: 16px;",
        title
      )
    },
    options = list(
      dom = "t",
      paging = FALSE,
      ordering = FALSE,
      autoWidth = TRUE,
      
      rowCallback = JS("
        function(row, data, index) {
          // index parte da 0
          if (index == 17 || index == 22 || index == 31|| index == 33 || index == 35) {
            $('td', row).css('border-bottom', '2px dashed #666');
          }
        }
      ")
      
    ),
    class = "compact"
  )
}

################################################ build summary function for all the indices

###############################################################################################
# Builds a summary table from a list of evaluation tables
# Input  : tabs_list   -> list of tables containing backtesting p-values and FZLoss information
#          lab_full    -> vector of model names used to reorder and align rows across tables
#          as_percent  -> logical value; if TRUE, returns percentages, otherwise raw counts
# Output : data frame reporting, for each model, the proportion or count of cases in which:
#          all six backtests are passed, the model belongs to the MCS, and both conditions hold
###############################################################################################

build_summary_from_list <- function(tabs_list,
                                    lab_full,
                                    as_percent = TRUE) {
  stopifnot(length(tabs_list) >= 1)

  pat_mcs <- "cellcolor\\{gray!50\\}"

  bt_count     <- NULL
  mcs_count    <- NULL
  both_count   <- NULL

  clean_num <- function(x) {
    x2 <- gsub("cellcolor\\{gray!25\\}|cellcolor\\{gray!50\\}", "", x)
    x2 <- trimws(x2)
    suppressWarnings(as.numeric(x2))
  }

  clean_lab <- function(x) {
    x2 <- gsub("cellcolor\\{gray!25\\}|cellcolor\\{gray!50\\}", "", x)
    trimws(x2)
  }

  for (tab in tabs_list) {

    if (is.null(rownames(tab))) {
      stop("Each table in 'tabs_list' must have row names.")
    }

    #### Reorder rows according to lab_full
    tab_labs_clean <- clean_lab(rownames(tab))
    ord <- match(lab_full, tab_labs_clean)

    if (any(is.na(ord))) {
      stop("Some models in 'lab_full' were not found in one of the tables.")
    }

    tab <- tab[ord, , drop = FALSE]

    tab_chr <- as.matrix(tab)
    mode(tab_chr) <- "character"

    #### Backtests: first 6 columns must all be > 0.05
    bt_block <- apply(tab_chr[, 1:6, drop = FALSE], 2, clean_num)
    if (is.null(dim(bt_block))) bt_block <- matrix(bt_block, ncol = 1)

    has_bt <- apply(bt_block, 1, function(r) all(r > 0.05, na.rm = FALSE))

    #### MCS: detect gray!50 in FZLoss
    cn <- colnames(tab)
    fz_cols_any <- which(grepl("FZ\\s*Loss", cn, ignore.case = TRUE))
    if (length(fz_cols_any) < 1) stop("FZLoss column not found.")

    has_mcs <- apply(
      tab_chr[, fz_cols_any, drop = FALSE],
      1,
      function(r) any(grepl(pat_mcs, r, perl = TRUE))
    )

    #### Initialize accumulators
    if (is.null(bt_count)) {
      nmod <- length(lab_full)
      bt_count   <- integer(nmod)
      mcs_count  <- integer(nmod)
      both_count <- integer(nmod)
    }

    bt_count   <- bt_count   + as.integer(has_bt)
    mcs_count  <- mcs_count  + as.integer(has_mcs)
    both_count <- both_count + as.integer(has_bt & has_mcs)
  }

  denom <- length(tabs_list)

  if (as_percent) {
    out <- data.frame(
      `BT_pass(%)` = 100 * bt_count / denom,
      `MCS_in(%)`  = 100 * mcs_count / denom,
      `BT&MCS(%)`  = 100 * both_count / denom,
      row.names = lab_full,
      check.names = FALSE
    )
  } else {
    out <- data.frame(
      BT_pass    = bt_count,
      MCS_in     = mcs_count,
      BT_and_MCS = both_count,
      row.names = lab_full,
      check.names = FALSE
    )
  }

  out
}

################################## htlm version of Table 7 

###############################################################################################
# Creates the DT::datatable version of Table 7
# Input  : x -> matrix or data frame containing model-level summary measures for the two
#               quantile levels considered in Table 7
# Output : DT::datatable object reproducing Table 7 in HTML format, with conditional formatting
#          highlighting the best and second-best values in each column
###############################################################################################

tab_7_dt_f <- function(x, title=NULL) {
 
  # Assign column names
  colnames(x) <- c(
    "BT_pass_0.025", "MCS_in_0.025", "BT_and_MCS_0.025",
    "BT_pass_0.01",  "MCS_in_0.01",  "BT_and_MCS_0.01"
  )
  
  # Compute best and second-best values by column
  max_vals <- apply(x, 2, max, na.rm = TRUE)
  
  # 
  second_vals <- apply(x, 2, function(x) {
    ux <- sort(unique(x), decreasing = TRUE)
    if (length(ux) >= 2) ux[2] else NA
  })
  
  tab_dt <- datatable(
    x,
    rownames = TRUE,
    class = "compact",
    caption = if (!is.null(title)) {
      htmltools::tags$caption(
        style = "caption-side: top; text-align: center; font-weight: bold; font-size: 16px;",
        title
      )
    },
    options = list(
      dom = "t",
      ordering = FALSE,
      paging = FALSE,
      searching = FALSE,
      info = FALSE,
      autoWidth = TRUE,
      rowCallback = JS(
        sprintf("
        function(row, data, index) {
          
          var max_vals = [%s];
          var second_vals = [%s];
          
          if (index == 17 || index == 22 || index == 31 || index == 33 || index == 35) {
            $('td', row).css('border-bottom', '2px dashed #666');
          }
          
          for (var i = 0; i < max_vals.length; i++) {
            
            var val = parseFloat(data[i + 1]);
            var cell = $('td:eq(' + (i + 1) + ')', row);
            
            // best 
            if (!isNaN(val) && val === max_vals[i]) {
              cell.css('background-color', '#39e639');
            }
            
            // second-best 
            else if (!isNaN(val) && val === second_vals[i]) {
              cell.css('background-color', '#d9f7d9');
            }
            
            // bold 
            if (!isNaN(val) && val === 9) {
              cell.css('font-weight', 'bold');
            }
          }
        }
        ",
        paste(max_vals, collapse = ","),
        paste(second_vals, collapse = ","))
      )
    )
  )
  
  return(tab_dt)
}

#################################### QL

###############################################################################################
# Computes the quantile loss (tick loss) function
# Input  : ret   -> vector of returns
#          VAR   -> Value-at-Risk forecasts
#          alpha -> quantile level
# Output : vector of quantile loss values
###############################################################################################

QL<-function(ret,VAR,alpha){

Ind<-ifelse(ret<=VAR,1,0)

return((ret-VAR)*(alpha-Ind))

}

#################################### Print method for "rqmidas" class

###############################################################################################
# Summary method for objects of class "rqmidas"
# Input  : object -> fitted rqmidas object
# Output : model summary including coefficients, loss and VaR violations
###############################################################################################

print.rqmidas <- function(x, ...) {

cat(utils::capture.output(x$coef_mat),  sep = '\n')

}

#################################### Summary method for "rqmidas" class

###############################################################################################
# Summary method for objects of class "rqmidas"
# Input  : object -> fitted rqmidas object
# Output : model summary including coefficients, loss and VaR violations
###############################################################################################

summary.rqmidas <- function(object, ...) {

model<-object$model
mat_coef<-object$coef_mat
Obs<-object$obs
Period<-object$period
loss<-round(100*mean(object$loss_in_s),6)
var_viol<-round(object$hit_in_s*100,6)


Period<-paste(substr(Period[1], 1, 10),"/",
substr(Period[2], 1, 10),sep="")


cat(
cat("\n"),
cat("Coefficients:\n"),
cat(utils::capture.output(mat_coef),  sep = '\n'),
cat("--- \n"),
cat("Obs.:", paste(Obs, ".",sep=""), "Sample Period:", Period, "\n"),
cat("QL(%): ", loss, "\n"),
cat("VaR Viol.(%): ", var_viol, "\n"),
cat("\n"))


}

#################################### Sequential test for finding the optimal number of daily returns' lags 

###############################################################################################
# Performs a sequential likelihood-ratio test to determine the optimal number
# of autoregressive lags in the MF-QR-X models
# More details on Candila et al. (Annals of Operational Research, 2023)
# Link: https://link.springer.com/article/10.1007/s10479-023-05370-x
# Input  : Y    -> return series
#          tau  -> quantile level
#          mv_m -> MIDAS variable matrix (optional)
#          K    -> number of MIDAS lags (optional)
#          X    -> exogenous variable (optional)
# Output : list of p-values associated with sequential LR tests
###############################################################################################

seq_test<-function(Y,tau,mv_m=NULL,K=NULL,X=NULL){

N<-length(Y)

res<-list(NA,NA,NA,NA)

res<-setNames(res,c("Linear ARCH","Linear ARCH-MIDAS",
"Linear ARCH-X","Linear ARCH-MIDAS-X"))
######################################################

X_t_test<-cbind(
beta_0=rep(1,N),
beta_1=lag(abs(Y),k=1),
beta_2=lag(abs(Y),k=2),
beta_3=lag(abs(Y),k=3),
beta_4=lag(abs(Y),k=4),
beta_5=lag(abs(Y),k=5),
beta_6=lag(abs(Y),k=6),
beta_7=lag(abs(Y),k=7),
beta_8=lag(abs(Y),k=8),
beta_9=lag(abs(Y),k=9),
beta_10=lag(abs(Y),k=10)
)

X_t_test<-replace(X_t_test,is.na(X_t_test),0)

#################################### 
#################################### no MV and X variables
#################################### 

if(missing(mv_m)&missing(X)){

seq_test<-matrix(rep(NA,(ncol(X_t_test)-2)),ncol=(ncol(X_t_test)-2))

colnames(seq_test)<-paste("Beta","_",2:(ncol(X_t_test)-1),"=0",sep="")

#################################### unrestricted and restrict. model

for (COL in (ncol(X_t_test)-1):2){

X_t_res		<-	X_t_test[,2:(ncol(X_t_test)+1-COL)]	
X_t_unres	<-	X_t_test[,2:(ncol(X_t_test)+2-COL)]

res_mod<-quantreg::rq(formula = Y ~ X_t_res, tau = tau)
unres_mod<-quantreg::rq(formula = Y ~ X_t_unres, tau = tau)

v_tau<-unres_mod$rho			# rho uncostrained
v_tilde<-res_mod$rho			# rho costrained

h_n<-quantreg::bandwidth.rq(tau,N)

x_bar_red<-c(1,colMeans(X_t_res)) #restricted

Q_tau_plus_h_n<-rbind(x_bar_red)%*%rq(formula = Y ~ X_t_res, 
tau = tau+h_n)$coefficients

Q_tau_minus_h_n<-rbind(x_bar_red)%*%rq(formula = Y ~ X_t_res, 
tau = tau-h_n)$coefficients

s_t<-(Q_tau_plus_h_n-Q_tau_minus_h_n)/(2*h_n)

LR<-2*(v_tilde-v_tau)/(tau*(1-tau)*s_t)

p_value<-pchisq(LR,df=1,lower.tail=FALSE)

seq_test[,c(0,(ncol(X_t_test)-2):1)[COL]]<-p_value
}

res[[1]]<-seq_test

}
#################################### 
#################################### only MV variable
#################################### 

if(!missing(mv_m)&missing(X)){

seq_test<-matrix(rep(NA,(ncol(X_t_test)-2)),ncol=(ncol(X_t_test)-2))

colnames(seq_test)<-paste("Beta","_",2:(ncol(X_t_test)-1),"=Theta","=0",sep="")

for (COL in (ncol(X_t_test)-1):2){

X_t_res		<-	X_t_test[,2:(ncol(X_t_test)+1-COL)]	
Q<-length(ncol(X_t_res))+1

res_mod<-quantreg::rq(formula = Y ~ X_t_res, tau = tau)
unres_mod<-uqfit(model="lARCHMIDAS",tau=tau,daily_ret=Y,mv_m=mv_m,K=K,q=Q,B=1)

v_tau<-sum(unres_mod$loss_in_s)			# rho uncostrained
v_tilde<-res_mod$rho			# rho costrained

h_n<-quantreg::bandwidth.rq(tau,N)

x_bar_red<-c(1,colMeans(X_t_res)) #restricted

Q_tau_plus_h_n<-rbind(x_bar_red)%*%rq(formula = Y ~ X_t_res, 
tau = tau+h_n)$coefficients

Q_tau_minus_h_n<-rbind(x_bar_red)%*%rq(formula = Y ~ X_t_res, 
tau = tau-h_n)$coefficients

s_t<-(Q_tau_plus_h_n-Q_tau_minus_h_n)/(2*h_n)

LR<-2*(v_tilde-v_tau)/(tau*(1-tau)*s_t)

p_value<-pchisq(LR,df=2,lower.tail=FALSE)

seq_test[,c(0,(ncol(X_t_test)-2):1)[COL]]<-p_value
}

res[[2]]<-seq_test
}

#################################### 
#################################### only X variable
#################################### 

if(missing(mv_m)&!missing(X)){

seq_test<-matrix(rep(NA,(ncol(X_t_test)-2)),ncol=(ncol(X_t_test)-2))

colnames(seq_test)<-paste("Beta","_",2:(ncol(X_t_test)-1),"=z","=0",sep="")

for (COL in (ncol(X_t_test)-1):2){

X_t_res		<-	X_t_test[,2:(ncol(X_t_test)+1-COL)]	
Q<-length(ncol(X_t_res))+1

res_mod<-quantreg::rq(formula = Y ~ X_t_res, tau = tau)
unres_mod<-uqfit(model="lARCHX",tau=tau,daily_ret=Y,q=Q,B=1,X=X)

v_tau<-sum(unres_mod$loss_in_s)			# rho uncostrained
v_tilde<-res_mod$rho			# rho costrained

h_n<-quantreg::bandwidth.rq(tau,N)

x_bar_red<-c(1,colMeans(X_t_res)) #restricted

Q_tau_plus_h_n<-rbind(x_bar_red)%*%rq(formula = Y ~ X_t_res, 
tau = tau+h_n)$coefficients

Q_tau_minus_h_n<-rbind(x_bar_red)%*%rq(formula = Y ~ X_t_res, 
tau = tau-h_n)$coefficients

s_t<-(Q_tau_plus_h_n-Q_tau_minus_h_n)/(2*h_n)

LR<-2*(v_tilde-v_tau)/(tau*(1-tau)*s_t)

p_value<-pchisq(LR,df=2,lower.tail=FALSE)

seq_test[,c(0,(ncol(X_t_test)-2):1)[COL]]<-p_value
}

res[[3]]<-seq_test
}

#################################### 
#################################### MV + X variables
#################################### 

if(!missing(mv_m)&!missing(X)){

seq_test<-matrix(rep(NA,(ncol(X_t_test)-2)),ncol=(ncol(X_t_test)-2))

colnames(seq_test)<-paste("Beta","_",2:(ncol(X_t_test)-1),"=Theta=z","=0",sep="")

for (COL in (ncol(X_t_test)-1):2){

X_t_res		<-	X_t_test[,2:(ncol(X_t_test)+1-COL)]	
Q<-length(ncol(X_t_res))+1

res_mod<-quantreg::rq(formula = Y ~ X_t_res, tau = tau)
unres_mod<-uqfit(model="lARCHMIDASX",tau=tau,daily_ret=Y,mv_m=mv_m,K=K,q=Q,B=1,X=X)

v_tau<-sum(unres_mod$loss_in_s)			# rho uncostrained
v_tilde<-res_mod$rho			# rho costrained

h_n<-quantreg::bandwidth.rq(tau,N)

x_bar_red<-c(1,colMeans(X_t_res)) #restricted

Q_tau_plus_h_n<-rbind(x_bar_red)%*%rq(formula = Y ~ X_t_res, 
tau = tau+h_n)$coefficients

Q_tau_minus_h_n<-rbind(x_bar_red)%*%rq(formula = Y ~ X_t_res, 
tau = tau-h_n)$coefficients

s_t<-(Q_tau_plus_h_n-Q_tau_minus_h_n)/(2*h_n)

LR<-2*(v_tilde-v_tau)/(tau*(1-tau)*s_t)

p_value<-pchisq(LR,df=3,lower.tail=FALSE)

seq_test[,c(0,(ncol(X_t_test)-2):1)[COL]]<-p_value
}

res[[4]]<-seq_test
}


return(res)

}


######################################## Indirect GARCH estimation function

###############################################################################################
# Computes the QL for the Indirect GARCH model
# Input  : param -> vector of model parameters
#          ret   -> return series
#          tau   -> quantile level
# Output : vector of QL contributions
###############################################################################################


ig_fun <- function(param, ret, tau){

b1 <- param[1] 
b2 <- param[2] 
b3 <- param[3]

N<-length(ret)
VaR_t <- rep(NA,N) 
VaR_t[1] <- quantile(ret,tau)

    for(i in 2:N) {
      VaR_t[i] <- c(- (b1 + b2 * VaR_t[i-1]^2 + b3 * ret[i-1]^2)^(0.5))
    }

    est<- -QL(ret,VaR_t,tau)
    est
  }

######################################## Indirect GARCH estimation function via ALD with ES

###############################################################################################
# Computes the ALD log-likelihood for the Indirect GARCH CAViaR model with Expected Shortfall
# Input  : param -> vector of model parameters
#          ret   -> return series
#          tau   -> quantile level
# Output : vector of ALD contributions
###############################################################################################

ig_fun_es <- function(param, ret, tau){

b1 <- param[1] 
b2 <- param[2] 
b3 <- param[3]
theta_gamma<-param[4]

N<-length(ret)
VaR_t<- rep(NA,N)
VaR_t[1] <- quantile(ret,tau)

 for(i in 2:N) {
      VaR_t[i] <-  c(- (b1 + b2 * VaR_t[i-1]^2 + b3 * ret[i-1]^2)^(0.5))
    }

es<-(1+exp(theta_gamma))*VaR_t
es<- -abs(es)
num_llk<-QL(ret,VaR_t,tau)
den_llk<-tau*es

llk<-(tau-1)*(es)^-1*exp(num_llk/den_llk)

return(log(llk))

}

######################################## Indirect GARCH function to obtain the Value-at-Risk

###############################################################################################
# Computes Value-at-Risk forecasts from the Indirect GARCH CAViaR model
# Input  : param -> vector of model parameters
#          ret   -> return series
#          tau   -> quantile level
# Output : vector of Value-at-Risk 
###############################################################################################

ig_var <- function(param, ret, tau){

b1 <- param[1] 
b2 <- param[2] 
b3 <- param[3]

N<-length(ret)
VaR_t <- rep(NA,N) 
VaR_t[1] <- quantile(ret,tau)

    for(i in 2:N) {
      VaR_t[i] <- c(- (b1 + b2 * VaR_t[i-1]^2 + b3 * ret[i-1]^2)^(0.5))
    }

    return(VaR_t)

}



################################# Symmetric Absolute Value estimation function

###############################################################################################
# Computes the QL for the SAV model 
# Input  : param -> vector of model parameters
#          ret   -> return series
#          tau   -> quantile level
# Output : vector of QL contributions
###############################################################################################

sav_fun <- function(param, ret, tau){

b0 <- param[1] 
b1 <- param[2] 
b2 <- param[3]

N<-length(ret)
VaR_t<- rep(NA,N)
VaR_t[1] <- quantile(ret,tau)

 for(i in 2:N) {
      VaR_t[i] <- b0 + b1 * VaR_t[i-1]  + b2*abs(ret[i-1])
    }

    est<- -QL(ret,VaR_t,tau)
    est
  }

########################### Symmetric Absolute Value estimation function via ALD with ES

###############################################################################################
# Computes the ALD log-likelihood for the SAV CAViaR model with Expected Shortfall
# Input  : param -> vector of model parameters
#          ret   -> return series
#          tau   -> quantile level
# Output : vector of ALD contributions
###############################################################################################

sav_fun_es <- function(param, ret, tau){

b0 <- param[1] 
b1 <- param[2] 
b2 <- param[3]
theta_gamma<-param[4]

N<-length(ret)
VaR_t<- rep(NA,N)
VaR_t[1] <- quantile(ret,tau)

 for(i in 2:N) {
      VaR_t[i] <- b0 + b1 * VaR_t[i-1]  + b2*abs(ret[i-1])
    }

es<-(1+exp(theta_gamma))*VaR_t
es<- -abs(es)
num_llk<-QL(ret,VaR_t,tau)
den_llk<-tau*es

llk<-(tau-1)*(es)^-1*exp(num_llk/den_llk)

return(log(llk))

}



# Symmetric Absolute Value function to obtain the Value-at-Risk

###############################################################################################
# Computes Value-at-Risk forecasts from the SAV CAViaR model
# Input  : param -> vector of model parameters
#          ret   -> return series
#          tau   -> quantile level
# Output : vector of Value-at-Risk 
###############################################################################################

sav_var <- function(param, ret, tau){

b0 <- param[1] 
b1 <- param[2] 
b2 <- param[3]

N<-length(ret)
VaR_t<- rep(NA,N)
VaR_t[1] <- quantile(ret,tau)

 for(i in 2:N) {
      VaR_t[i] <- b0 + b1 * VaR_t[i-1]  + b2*abs(ret[i-1])
    }

    return(VaR_t)

  }

############################## Symmetric Absolute Value with X variable estimation function


###############################################################################################
# Computes the QL for the SAV-X model 
# Input  : param -> vector of model parameters
#          ret   -> return series
#          tau   -> quantile level
# Output : vector of QL contributions
###############################################################################################

sav_fun_x <- function(param, ret, tau, X){

b0 <- param[1] 
b1 <- param[2] 
b2 <- param[3]


N<-length(ret)
VaR_t<- rep(NA,N)
VaR_t[1] <- quantile(ret,tau)

 for(i in 2:N) {
      VaR_t[i] <- b0 + b1 * VaR_t[i-1]  + b2*X[i-1]
    }

    est<- -QL(ret,VaR_t,tau)
    est
  }

############################## SAV-X estimation function via ALD with ES

###############################################################################################
# Computes the ALD log-likelihood for the SAV-X model with Expected Shortfall
# Input  : param -> vector of model parameters
#          ret   -> return series
#          tau   -> quantile level
# Output : vector of ALD contributions
###############################################################################################

sav_fun_es_x <- function(param, ret, tau, X){

b0 <- param[1] 
b1 <- param[2] 
b2 <- param[3]
theta_gamma<-param[4]

N<-length(ret)
VaR_t<- rep(NA,N)
VaR_t[1] <- quantile(ret,tau)

 for(i in 2:N) {
      VaR_t[i] <- b0 + b1 * VaR_t[i-1] + b2*X[i-1]
    }

es<-(1+exp(theta_gamma))*VaR_t
es<- -abs(es)
num_llk<-QL(ret,VaR_t,tau)
den_llk<-tau*es

llk<-(tau-1)*(es)^-1*exp(num_llk/den_llk)

return(log(llk))

}

########################################### SAV-X function to obtain the Value-at-Risk

###############################################################################################
# Computes Value-at-Risk forecasts from the SAV-X model
# Input  : param -> vector of model parameters
#          ret   -> return series
#          tau   -> quantile level
# Output : vector of Value-at-Risk 
###############################################################################################


sav_var_x <- function(param, ret, tau, X){

b0 <- param[1] 
b1 <- param[2] 
b2 <- param[3]

N<-length(ret)
VaR_t<- rep(NA,N)
VaR_t[1] <- quantile(ret,tau)

 for(i in 2:N) {
      VaR_t[i] <- b0 + b1 * VaR_t[i-1] + b2*X[i-1]
    }

    return(VaR_t)

  }

# Asymmetric Slope estimation function

as_fun <- function(param, ret, tau){

b0 <- param[1] 
b1 <- param[2] 
b2 <- param[3]
b3 <- param[4]

N<-length(ret)
VaR_t<- rep(NA,N); 
VaR_t[1] <- quantile(ret,tau)

ind_pos<-ifelse(ret>0,1,0)
ind_neg<-ifelse(ret<0,1,0)

    for(i in 2:N) {
      VaR_t[i] <- b0 + b1 * VaR_t[i-1] + (
		b2*ind_pos[i-1] + b3*ind_neg[i-1]
		)*abs(ret[i-1])
    }


    est<- -QL(ret,VaR_t,tau)
    #est<- mean(est)
    est
  }


########################################### Asymmetric Slope estimation function to obtain the Value-at-Risk

###############################################################################################
# Computes Value-at-Risk forecasts from the AS CAViaR model
# Input  : param -> vector of model parameters
#          ret   -> return series
#          tau   -> quantile level
# Output : vector of Value-at-Risk 
###############################################################################################


as_var <- function(param, ret, tau){

b0 <- param[1] 
b1 <- param[2] 
b2 <- param[3]
b3 <- param[4]

N<-length(ret)
VaR_t<- rep(NA,N)
VaR_t[1] <- quantile(ret,tau)

ind_pos<-ifelse(ret>0,1,0)
ind_neg<-ifelse(ret<0,1,0)

    for(i in 2:N) {
      VaR_t[i] <- b0 + b1 * VaR_t[i-1] + (
		b2*ind_pos[i-1] + b3*ind_neg[i-1]
		)*abs(ret[i-1])
    }


   return(VaR_t)

  }

########################################### Asymmetric Slope estimation function via ALD with ES

###############################################################################################
# Computes the ALD log-likelihood for the AS CAViaR model with Expected Shortfall
# Input  : param -> vector of model parameters
#          ret   -> return series
#          tau   -> quantile level
# Output : vector of ALD contributions
###############################################################################################

as_fun_es <- function(param, ret, tau){

b0 <- param[1] 
b1 <- param[2] 
b2 <- param[3]
b3 <- param[4]
theta_gamma<-param[5]

N<-length(ret)
VaR_t<- rep(NA,N) 
VaR_t[1] <- quantile(ret,tau)

ind_pos<-ifelse(ret>0,1,0)
ind_neg<-ifelse(ret<0,1,0)

    for(i in 2:N) {
      VaR_t[i] <- b0 + b1 * VaR_t[i-1] + (
		b2*ind_pos[i-1] + b3*ind_neg[i-1]
		)*abs(ret[i-1])
    }

es<-(1+exp(theta_gamma))*VaR_t
es<- -abs(es)
num_llk<-QL(ret,VaR_t,tau)
den_llk<-tau*es

llk<-(tau-1)*(es)^-1*exp(num_llk/den_llk)

return(log(llk))
  
}

########################################### Methods for estimating the CAViaR models

###############################################################################################
# Estimates a variety of CAViaR-type models (only VaR or VaR and ES)
# Input  : model         -> model specification ("SAV", "SAVX", "AS", "IG")
#          tau           -> quantile level
#          daily_ret     -> return series
#          X             -> exogenous variable (optional)
#          out_of_sample -> out-of-sample size (optional)
#          std_err       -> standard error estimation method
#          B             -> bootstrap replications
#          R             -> random starting values
#          Exp_Short     -> whether Expected Shortfall is estimated
# Output : object of class "rqmidas"
###############################################################################################

ucaviarfit<-function(
model,
tau,
daily_ret,
X=NULL,
out_of_sample=NULL,
std_err="bootstrap",
B=100,
R=100,
Exp_Short='No'){

############################# check on valid choices
if(tau<0|tau>1) { stop(cat("#Warning:\n Parameter 'tau' must be a number strictly included between 0 and 1 \n"))}

if(!is.null(X)&length(daily_ret)!=length(X)) { stop(cat("#Warning:\n 'daily_ret' and 'X' must have the same length \n"))}

if (!(model %in% c("SAV", "SAVX", "AS", "IG"))) {
  stop("#Warning:\n Valid choices for the parameter 'model' are currently 'SAV', 'SAVX', 'AS', and 'IG'")
}

N<-length(daily_ret)

VaR_oos<-NA
ES_oos<-NA

if (missing(out_of_sample)){ # out_of_sample parameter is missing

r_t_in_s<-daily_ret
r_t_in_s_est<-zoo::coredata(daily_ret)
X_in_s_est<-zoo::coredata(X)


} else {					# out_of_sample parameter is present

r_t_in_s<-daily_ret[1:(N-out_of_sample)]
X_in_s_est<-X[1:(N-out_of_sample)]
X_in_s_est<-zoo::coredata(X_in_s_est)
r_t_in_s_est<-zoo::coredata(r_t_in_s)

r_t_oos<-daily_ret[(N-out_of_sample+1):N]
#r_t_oos<-zoo::coredata(r_t_oos)
X_oos<-X[(N-out_of_sample+1):N]
#X_oos<-zoo::coredata(X_oos)

N<-length(r_t_in_s)

}


############################################### Estimation

if (model=="SAV"&Exp_Short=="No"){

begin_val<-matrix(NA,ncol=3,nrow=R)
begin_val[,1]<-runif(R,-0.3,0.1)
begin_val[,2]<-runif(R,0.5,0.98)
begin_val[,3]<-runif(R,-0.3,0.1)

colnames(begin_val)<-c("beta_0","beta_1","beta_2")
which_row<-rep(NA,R)

for(i in 1:nrow(begin_val)){
which_row[i]<-sum(sav_fun(begin_val[i,],r_t_in_s_est,tau))
}

ui<-ci<-NULL
ui<-rbind(c(0,-1,0))		## b1 < 0.9999
ci<-c(0.9999)

est<-suppressWarnings(maxLik(
logLik=sav_fun,
start=begin_val[which.max(which_row),],
ret=r_t_in_s_est,
tau=tau,
constraints=list(ineqA=ui, ineqB=ci),
iterlim=1000,
method="BFGS"))

VaR_est <- sav_var(stats::coef(est),r_t_in_s_est,tau)

if (!missing(out_of_sample)) {
VaR_oos <- sav_var(stats::coef(est),r_t_oos,tau)
VaR_oos <- xts::as.xts(VaR_oos,stats::time(r_t_oos))
}

} else if (model=="SAV"&Exp_Short=="Yes"){


begin_val<-matrix(NA,ncol=3,nrow=R)
begin_val[,1]<-runif(R,-0.3,0.1)
begin_val[,2]<-runif(R,0.5,0.98)
begin_val[,3]<-runif(R,-0.3,0.1)

colnames(begin_val)<-c("beta_0","beta_1","beta_2")
which_row<-rep(NA,R)

for(i in 1:nrow(begin_val)){
which_row[i]<-sum(sav_fun(begin_val[i,],r_t_in_s_est,tau))
}

ui<-ci<-NULL
ui<-rbind(c(0,-1,0))		## b1 < 0.9999
ci<-c(0.9999)

est_0<-suppressWarnings(maxLik(
logLik=sav_fun,
start=begin_val[which.max(which_row),],
ret=r_t_in_s_est,
tau=tau,
constraints=list(ineqA=ui, ineqB=ci),
iterlim=1000,
method="BFGS"))

max_start_value<-matrix(NA,ncol=4,nrow=R)
sim_llk<-rep(NA,R)

for(i in 1:R){
max_start_value[i,]<-c(stats::coef(est_0),runif(1,-1,1))
sim_llk[i]<-sum(sav_fun_es(max_start_value[i,],r_t_in_s_est,tau))
}

start_val<-max_start_value[which.max(sim_llk),]

names(start_val)<-c("beta_0","beta_1","beta_2","gamma")

ui<-ci<-NULL
ui<-rbind(c(0,-1,0,0))		## b1 < 0.9999
ci<-c(0.9999)

est<- maxLik(sav_fun_es, 
start=start_val, 
ret=r_t_in_s_est,
tau=tau,
constraints=list(ineqA=ui, ineqB=ci),
iterlim=1000,
method="BFGS")

VaR_est <- sav_var(stats::coef(est)[1:3],r_t_in_s_est,tau)
ES_est <- (1+exp(stats::coef(est)[4]))*VaR_est

if (!missing(out_of_sample)) {
VaR_oos <- sav_var(stats::coef(est)[1:3],r_t_oos,tau)
VaR_oos <- xts::as.xts(VaR_oos,stats::time(r_t_oos))
ES_oos <- (1+exp(stats::coef(est)[4]))*VaR_oos
}


} else if (model=="SAVX"&Exp_Short=="No"){

begin_val<-matrix(NA,ncol=3,nrow=R)
begin_val[,1]<-runif(R,-0.3,0.1)
begin_val[,2]<-runif(R,0.5,0.98)
begin_val[,3]<-runif(R,-0.3,0.1)

colnames(begin_val)<-c("beta_0","beta_1","beta_X")
which_row<-rep(NA,R)

for(i in 1:nrow(begin_val)){
which_row[i]<-sum(sav_fun_x(begin_val[i,],r_t_in_s_est,tau,X_in_s_est))
}

ui<-ci<-NULL
ui<-rbind(c(0,-1,0))		## b1 < 0.9999
ci<-c(0.9999)

est<-suppressWarnings(maxLik(
logLik=sav_fun_x,
start=begin_val[which.max(which_row),],
ret=r_t_in_s_est,
tau=tau,
X=X_in_s_est,
constraints=list(ineqA=ui, ineqB=ci),
iterlim=1000,
method="BFGS"))

VaR_est <- sav_var_x(stats::coef(est),r_t_in_s_est,tau,X_in_s_est)

if (!missing(out_of_sample)) {
VaR_oos <- sav_var_x(stats::coef(est),r_t_oos,tau,X_oos)
VaR_oos <- xts::as.xts(VaR_oos,stats::time(r_t_oos))
} 

} else if (model=="SAVX"&Exp_Short=="Yes"){

begin_val<-matrix(NA,ncol=3,nrow=R)
begin_val[,1]<-runif(R,-0.3,0.1)
begin_val[,2]<-runif(R,0.5,0.98)
begin_val[,3]<-runif(R,-0.3,0.1)

colnames(begin_val)<-c("beta_0","beta_1","beta_X")
which_row<-rep(NA,R)

for(i in 1:nrow(begin_val)){
which_row[i]<-sum(sav_fun_x(begin_val[i,],r_t_in_s_est,tau,X_in_s_est))
}

ui<-ci<-NULL
ui<-rbind(c(0,-1,0))		## b1 < 0.9999
ci<-c(0.9999)

est_0<-suppressWarnings(maxLik(
logLik=sav_fun_x,
start=begin_val[which.max(which_row),],
ret=r_t_in_s_est,
tau=tau,
X=X_in_s_est,
method="BFGS"))

max_start_value<-matrix(NA,ncol=4,nrow=R)
sim_llk<-rep(NA,R)

for(i in 1:R){
max_start_value[i,]<-c(stats::coef(est_0),runif(1,-1,1))
sim_llk[i]<-sum(sav_fun_es_x(max_start_value[i,],r_t_in_s_est,tau,X_in_s_est))
}

start_val<-max_start_value[which.max(sim_llk),]

names(start_val)<-c("beta_0","beta_1","beta_X","gamma")


ui<-ci<-NULL
ui<-rbind(c(0,-1,0,0))		## b1 < 0.9999
ci<-c(0.9999)

est<- maxLik(sav_fun_es_x, 
start=start_val, 
ret=r_t_in_s_est,
tau=tau,
X=X_in_s_est,
constraints=list(ineqA=ui, ineqB=ci),
iterlim=1000,
method="BFGS")

VaR_est <- sav_var_x(stats::coef(est)[1:3],r_t_in_s_est,tau,X_in_s_est)
ES_est <- (1+exp(stats::coef(est)[4]))*VaR_est

if (!missing(out_of_sample)) {
VaR_oos <- sav_var_x(stats::coef(est)[1:3],r_t_oos,tau,X_oos)
VaR_oos <- xts::as.xts(VaR_oos,stats::time(r_t_oos))
ES_oos <- (1+exp(stats::coef(est)[4]))*VaR_oos
}


} else if (model=="AS"&Exp_Short=="No"){

begin_val<-matrix(NA,ncol=4,nrow=R)
begin_val[,1]<-runif(R,-0.3,0.1)
begin_val[,2]<-runif(R,0.5,0.98)
begin_val[,3]<-runif(R,-0.3,0.1)
begin_val[,4]<-runif(R,-0.3,0.1)

colnames(begin_val)<-c("beta_0","beta_1","beta_2","beta_3")
which_row<-rep(NA,R)

for(i in 1:nrow(begin_val)){
which_row[i]<-sum(as_fun(begin_val[i,],r_t_in_s_est,tau))
}

ui<-ci<-NULL
ui<-rbind(c(0,-1,0,0))		## b1 < 0.9999
ci<-c(0.9999)

est<- suppressWarnings(maxLik(
logLik=as_fun,
start=begin_val[which.max(which_row),],
ret=r_t_in_s_est,
tau=tau,
constraints=list(ineqA=ui, ineqB=ci),
iterlim=1000,
method="BFGS"))

VaR_est <- as_var(stats::coef(est),r_t_in_s_est,tau)

if (!missing(out_of_sample)) {
VaR_oos <- as_var(stats::coef(est),r_t_oos,tau)
VaR_oos <- xts::as.xts(VaR_oos,stats::time(r_t_oos))
}

} else if (model=="AS"&Exp_Short=="Yes"){

begin_val<-matrix(NA,ncol=4,nrow=R)
begin_val[,1]<-runif(R,-0.3,0.1)
begin_val[,2]<-runif(R,0.5,0.98)
begin_val[,3]<-runif(R,-0.3,0.1)
begin_val[,4]<-runif(R,-0.3,0.1)

colnames(begin_val)<-c("beta_0","beta_1","beta_2","beta_3")
which_row<-rep(NA,R)

for(i in 1:nrow(begin_val)){
which_row[i]<-sum(as_fun(begin_val[i,],r_t_in_s_est,tau))
}

ui<-ci<-NULL
ui<-rbind(c(0,-1,0,0))		## b1 < 0.9999
ci<-c(0.9999)

est_0<- suppressWarnings(maxLik(
logLik=as_fun,
start=begin_val[which.max(which_row),],
ret=r_t_in_s_est,
tau=tau,
constraints=list(ineqA=ui, ineqB=ci),
iterlim=1000,
method="BFGS"))

max_start_value<-matrix(NA,ncol=5,nrow=R)
sim_llk<-rep(NA,R)

for(i in 1:R){
max_start_value[i,]<-c(stats::coef(est_0),runif(1,-1,1))
sim_llk[i]<-sum(as_fun_es(max_start_value[i,],r_t_in_s_est,tau))
}

start_val<-max_start_value[which.max(sim_llk),]

names(start_val)<-c("beta_0","beta_1","beta_2","beta_3","gamma")

ui<-ci<-NULL
ui<-rbind(c(0,-1,0,0,0))		## b1 < 0.9999
ci<-c(0.9999)

est<- maxLik(as_fun_es, 
start=start_val, 
ret=r_t_in_s_est,
tau=tau,
constraints=list(ineqA=ui, ineqB=ci),
iterlim=1000,
method="BFGS")

VaR_est <- as_var(stats::coef(est)[1:4],r_t_in_s_est,tau)
ES_est <- (1+exp(stats::coef(est)[5]))*VaR_est

if (!missing(out_of_sample)) {
VaR_oos <- as_var(stats::coef(est)[1:4],r_t_oos,tau)
VaR_oos <- xts::as.xts(VaR_oos,stats::time(r_t_oos))
ES_oos <- (1+exp(stats::coef(est)[5]))*VaR_oos
}


} else if (model=="IG"&Exp_Short=="No"){
  
  begin_val<-matrix(NA,ncol=3,nrow=R)
  begin_val[,1]<-runif(R,0.01,0.3)
  begin_val[,2]<-runif(R,0.4,0.98)
  begin_val[,3]<-runif(R,0.01,0.35)
  
  colnames(begin_val)<-c("beta_0","beta_1","beta_2")
  which_row<-rep(NA,R)
  
  for(i in 1:nrow(begin_val)){
    which_row[i]<-sum(ig_fun(begin_val[i,],r_t_in_s_est,tau))
  }
  
  ui<-ci<-NULL
  
  ui<-rbind(
    c(-1,0,0),					## b0 < 0.9999
    c(0,-1,0),				 	## b1 < 0.9999
    c(0,0,-1))					## b2 < 0.9999
  
  ci<-c(0.9999,0.9999,0.9999)
  
  est<- suppressWarnings(maxLik(
    logLik=ig_fun,
    start=begin_val[which.max(which_row),],
    ret=r_t_in_s_est,
    tau=tau,
    constraints=list(ineqA=ui, ineqB=ci),
    iterlim=1000,
    method="BFGS"))
  
  VaR_est <- ig_var(stats::coef(est),r_t_in_s_est,tau)
  
  if (!missing(out_of_sample)) {
    VaR_oos <- ig_var(stats::coef(est),r_t_oos,tau)
    VaR_oos <- xts::as.xts(VaR_oos,stats::time(r_t_oos))
  }
  
} else if (model=="IG"&Exp_Short=="Yes"){

begin_val<-matrix(NA,ncol=3,nrow=R)
begin_val[,1]<-runif(R,0.01,0.3)
begin_val[,2]<-runif(R,0.4,0.98)
begin_val[,3]<-runif(R,0.01,0.35)

colnames(begin_val)<-c("beta_1","beta_2","beta_3")
which_row<-rep(NA,R)

for(i in 1:nrow(begin_val)){
which_row[i]<-sum(ig_fun(begin_val[i,],r_t_in_s_est,tau))
}

ui<-ci<-NULL

ui<-rbind(
c(-1,0,0),					## b0 < 0.9999
c(0,-1,0),				 	## b1 < 0.9999
c(0,0,-1))					## b2 < 0.9999

ci<-c(0.9999,0.9999,0.9999)

est_0<- suppressWarnings(maxLik(
logLik=ig_fun,
start=begin_val[which.max(which_row),],
ret=r_t_in_s_est,
tau=tau,
constraints=list(ineqA=ui, ineqB=ci),
iterlim=1000,
method="BFGS"))

max_start_value<-matrix(NA,ncol=4,nrow=R)
sim_llk<-rep(NA,R)

for(i in 1:R){
max_start_value[i,]<-c(stats::coef(est_0),runif(1,-1,1))
sim_llk[i]<-sum(ig_fun_es(max_start_value[i,],r_t_in_s_est,tau))
}

start_val<-max_start_value[which.max(sim_llk),]

names(start_val)<-c("beta_1","beta_2","beta_3","gamma")

ui<-ci<-NULL

ui<-rbind(
c(-1,0,0,0),					## b0 < 0.9999
c(0,-1,0,0),				 	## b1 < 0.9999
c(0,0,-1,0))					## b2 < 0.9999

ci<-c(0.9999,0.9999,0.9999)


est<- maxLik(ig_fun_es, 
start=start_val, 
ret=r_t_in_s_est,
tau=tau,
constraints=list(ineqA=ui, ineqB=ci),
iterlim=1000,
method="BFGS")

VaR_est <- ig_var(stats::coef(est)[1:3],r_t_in_s_est,tau)
ES_est <- (1+exp(stats::coef(est)[4]))*VaR_est

if (!missing(out_of_sample)) {
VaR_oos <- ig_var(stats::coef(est)[1:3],r_t_oos,tau)
VaR_oos <- xts::as.xts(VaR_oos,stats::time(r_t_oos))
ES_oos <- (1+exp(stats::coef(est)[4]))*VaR_oos
}

} 

N_coef<-length(stats::coef(est))

################################################# 
################################################# residuals
################################################# 

res<-r_t_in_s_est-VaR_est
stand_res<-(r_t_in_s_est-VaR_est)/abs(VaR_est)

N_est<-length(res)
mat_coef<-data.frame(rep(NA,N_coef),rep(NA,N_coef),rep(NA,N_coef),rep(NA,N_coef))
colnames(mat_coef)<-c("Estimate","Std. Error","t value","Pr(>|t|)")
rownames(mat_coef)<-names(stats::coef(est))

#### var-cov-matrix

if (std_err!="bootstrap"){
k_t_hat<-stats::mad(res)
m_t_hat<-quantreg::bandwidth.rq(tau,N_est)
c_t_hat<-k_t_hat*(qnorm(tau+m_t_hat)-qnorm(tau-m_t_hat))

GR<-est$gradientObs

A_T_hat<-N_est^-1*(tau)*(1-tau)*t(GR)%*%GR
IND<-ifelse(abs(r_t_in_s_est-VaR_est)<c_t_hat,1,0)
D_T_hat<-(2*N_est*c_t_hat)^(-1)*sum(IND*rowSums(GR^2))

sd_est<-(diag((N_est^0.5*A_T_hat%^%(-1/2)*D_T_hat))^-1)

} else { ## bootstrap standard errors

mat_boot_results<-matrix(rep(NA),nrow=B,ncol=N_coef)

pb <- utils::txtProgressBar(0, B, style = 3)

for(b in 1:B){

##### first step: compute the boostrapped series of residuals

res_boot<-sample(stand_res,N_est,replace=T)

##### second step: compute the bootstrapped series of returns and VaRs

r_t_boot<-r_t_in_s_est
VaR_boot<-VaR_est

if(model=="SAV"&Exp_Short=="No"){ ############################################# model SAV

b0<-stats::coef(est)[1]
b1<-stats::coef(est)[2]
b2<-stats::coef(est)[3]

for(i in 2:N_est){
VaR_boot[i]<-b0 + b1 * VaR_boot[i-1]  + b2*abs(r_t_boot[i-1])
r_t_boot[i]<-VaR_boot[i]+abs(VaR_boot[i])*res_boot[i]
}

est_boot<- suppressWarnings(maxLik(
logLik=sav_fun,
start=stats::coef(est),
ret=r_t_boot,
tau=tau,
method="NM"))

} else if (model=="SAV"&Exp_Short=="Yes"){ ############################################# model SAV + ES

b0<-stats::coef(est)[1]
b1<-stats::coef(est)[2]
b2<-stats::coef(est)[3]

for(i in 2:N_est){
VaR_boot[i]<-b0 + b1 * VaR_boot[i-1]  + b2*abs(r_t_boot[i-1])
r_t_boot[i]<-VaR_boot[i]+abs(VaR_boot[i])*res_boot[i]
}

est_boot<- suppressWarnings(maxLik(
logLik=sav_fun_es,
start=stats::coef(est),
ret=r_t_boot,
tau=tau,
method="NM"))

} else if(model=="SAVX"&Exp_Short=="No"){ ############################################# model SAV-X

b0<-stats::coef(est)[1]
b1<-stats::coef(est)[2]
bx<-stats::coef(est)[3]

for(i in 2:N_est){
VaR_boot[i]<-b0 + b1 * VaR_boot[i-1] + bx*X_in_s_est[i-1] 
r_t_boot[i]<-VaR_boot[i]+abs(VaR_boot[i])*res_boot[i]
}

est_boot<- suppressWarnings(maxLik(
logLik=sav_fun_x,
start=stats::coef(est),
ret=r_t_boot,
tau=tau,
X=X_in_s_est,
method="NM"))

} else if (model=="SAVX"&Exp_Short=="Yes"){ ############################################# model SAV-X + ES

  b0 <- stats::coef(est)[1]
  b1 <- stats::coef(est)[2]
  bx <- stats::coef(est)[3]

for(i in 2:N_est){
VaR_boot[i]<-b0 + b1 * VaR_boot[i-1]  + bx*X_in_s_est[i-1] 
r_t_boot[i]<-VaR_boot[i]+abs(VaR_boot[i])*res_boot[i]
}

est_boot<- suppressWarnings(maxLik(
logLik=sav_fun_es_x,
start=stats::coef(est),
ret=r_t_boot,
tau=tau,
X=X_in_s_est,
method="NM"))

} else if(model=="AS" & Exp_Short=="No"){ ############################################# model AS

b0<-stats::coef(est)[1]
b1<-stats::coef(est)[2]
b2<-stats::coef(est)[3]
b3<-stats::coef(est)[4]

for(i in 2:N_est){

ind_pos<-ifelse(r_t_boot[i-1]>0,1,0)
ind_neg<-ifelse(r_t_boot[i-1]<0,1,0)

VaR_boot[i]<- b0 + b1 * VaR_boot[i-1] + (
                b2*ind_pos + b3*ind_neg
                )*abs(r_t_boot[i-1])
r_t_boot[i]<-VaR_boot[i]+abs(VaR_boot[i])*res_boot[i]
}

est_boot<- suppressWarnings(maxLik(
logLik=as_fun,
start=stats::coef(est),
ret=r_t_boot,
tau=tau,
method="NM"))

} else if(model=="AS" & Exp_Short=="Yes"){ ############################################# model AS

b0<-stats::coef(est)[1]
b1<-stats::coef(est)[2]
b2<-stats::coef(est)[3]
b3<-stats::coef(est)[4]

for(i in 2:N_est){

ind_pos<-ifelse(r_t_boot[i-1]>0,1,0)
ind_neg<-ifelse(r_t_boot[i-1]<0,1,0)

VaR_boot[i]<- b0 + b1 * VaR_boot[i-1] + (
                b2*ind_pos + b3*ind_neg
                )*abs(r_t_boot[i-1])
r_t_boot[i]<-VaR_boot[i]+abs(VaR_boot[i])*res_boot[i]
}

est_boot<- suppressWarnings(maxLik(
logLik=as_fun_es,
start=stats::coef(est),
ret=r_t_boot,
tau=tau,
method="BFGS"))

} else if(model=="IG" & Exp_Short=="No"){ ############################################# model IG

b1<-stats::coef(est)[1]
b2<-stats::coef(est)[2]
b3<-stats::coef(est)[3]

for(i in 2:N_est){
VaR_boot[i]<-c(- (b1 + b2 * VaR_boot[i-1]^2 + b3 * r_t_boot[i-1]^2)^(0.5))
r_t_boot[i]<-VaR_boot[i]+abs(VaR_boot[i])*res_boot[i]
}

est_boot<- suppressWarnings(maxLik(
logLik=ig_fun,
start=stats::coef(est),
ret=r_t_boot,
tau=tau,
method="NM"))

} else if(model=="IG" & Exp_Short=="Yes"){ ############################################# model IG
  
  b1<-stats::coef(est)[1]
  b2<-stats::coef(est)[2]
  b3<-stats::coef(est)[3]
  
  for(i in 2:N_est){
    VaR_boot[i]<-c(- (b1 + b2 * VaR_boot[i-1]^2 + b3 * r_t_boot[i-1]^2)^(0.5))
    r_t_boot[i]<-VaR_boot[i]+abs(VaR_boot[i])*res_boot[i]
  }
  
  est_boot<- suppressWarnings(maxLik(
    logLik=ig_fun_es,
    start=stats::coef(est),
    ret=r_t_boot,
    tau=tau,
    method="NM"))
  
} 

##### third step: compute the new estimates (based on bootstrapped returns)

mat_boot_results[b,]<-stats::coef(est_boot)

utils::setTxtProgressBar(pb, b)
}  

Sys.sleep(1)
close(pb)

## compute the standard errors 

sd_est<-apply(mat_boot_results,2,stats::sd)

}




mat_coef[,1]<-round(stats::coef(est),6)
mat_coef[,2]<-round(sd_est,6)
mat_coef[,3]<-round(stats::coef(est)/sd_est,6)
mat_coef[,4]<-round(apply(rbind(stats::coef(est)/sd_est),1,function(x) 2*(1-stats::pnorm(abs(x)))),6)


###

if (Exp_Short=="No"){
res_f<-list(
model=model,
coef_mat=round(mat_coef,5),
obs=N,
period=range(stats::time(r_t_in_s)),
VaR=xts::as.xts(VaR_est,stats::time(r_t_in_s)),
hit_in_s=sum(ifelse(r_t_in_s<=VaR_est,1,0))/N,
loss_in_s=QL(r_t_in_s,VaR_est,tau),
VaR_oos=VaR_oos
)} else {
res_f<-list(
model=model,
coef_mat=round(mat_coef,5),
obs=N,
period=range(stats::time(r_t_in_s)),
VaR=xts::as.xts(VaR_est,stats::time(r_t_in_s)),
ES=xts::as.xts(ES_est,stats::time(r_t_in_s)),
hit_in_s=sum(ifelse(r_t_in_s<=VaR_est,1,0))/N,
loss_in_s=QL(r_t_in_s,VaR_est,tau),
VaR_oos=VaR_oos,
ES_oos=ES_oos
)

}

class(res_f)<-c("rqmidas")

return(res_f)
print.rqmidas(res_f)


}



########################################### Mixed-Frequency Quantile Regression with ES via ALD 

###############################################################################################
# Computes the ALD log-likelihood for the MF-QR model with Expected Shortfall
# Input  : param -> vector of model parameters
#          ret   -> return series
#          q     -> number of autoregressive lags
#          mv    -> MIDAS variable matrix
#          K     -> number of MIDAS lags
#          tau   -> quantile level
# Output : vector of ALD contributions
###############################################################################################


AR_lARCHMIDAS_fun_es <- function(param, ret, q, mv, K, tau){

N<-length(ret)

########################## parameter setting

theta<-param[q+2]
w1<-1
w2<-param[q+3]
z<-param[q+4]

########################## X_t

X_t<-cbind(
rep(1,N),
stats::lag(abs(ret),k=1:q)) 
X_t<-replace(X_t,is.na(X_t),0)

ret<-coredata(ret)
X_t<-coredata(X_t)

colnames(X_t)<-paste("beta","_",0:q,sep="")

########################## MIDAS part

betas<-c(rev(rumidas::beta_function(1:(K+1),(K+1),w1,w2))[2:(K+1)],0)
tau_d<- suppressWarnings(roll::roll_sum(mv_m, c(K+1),weights = betas))
tau_d<-abs(tau_d[(K+1),])	

######################### Cycle definition

VaR_t<- rep(NA,N)
VaR_t[1:q] <- quantile(ret,tau)

 for(i in (q+1):N) {
      VaR_t[i] <- X_t[i,]%*%cbind(param[1:(q+1)])+
			    theta*tau_d[i]+
                   z*VaR_t[i-1]
    }

es<-(1+exp(param[q+5]))*VaR_t
es<- -abs(es)
num_llk<-QL(ret,VaR_t,tau)
den_llk<-tau*es

llk<-(tau-1)*(es)^-1*exp(num_llk/den_llk)
llk<-llk[(q+1):N]
llk<-ifelse(is.na(llk)|llk<0|llk==-Inf|llk==0,0.000001,llk)

return(-sum(log(llk)))

}

########################################### VaR and ES for the Mixed-Frequency Quantile Regression  

###############################################################################################
# Computes Value-at-Risk and ES for the MF-QR model
# Input  : param -> vector of model parameters
#          ret   -> return series
#          q     -> number of autoregressive lags
#          mv    -> MIDAS variable matrix
#          K     -> number of MIDAS lags
#          tau   -> quantile level
# Output : matrix containing Value-at-Risk and Expected Shortfall forecasts
###############################################################################################


AR_lARCHMIDAS_var_es <- function(param, ret, q, mv, K, tau){

N<-length(ret)

########################## parameter setting

theta<-param[q+2]
w1<-1
w2<-param[q+3]
z<-param[q+4]

########################## X_t

X_t<-cbind(
rep(1,N),
stats::lag(abs(ret),k=1:q)) 
X_t<-replace(X_t,is.na(X_t),0)

ret<-coredata(ret)
X_t<-coredata(X_t)

colnames(X_t)<-paste("beta","_",0:q,sep="")

########################## MIDAS part

betas<-c(rev(rumidas::beta_function(1:(K+1),(K+1),w1,w2))[2:(K+1)],0)
tau_d<- suppressWarnings(roll::roll_sum(mv_m, c(K+1),weights = betas))
tau_d<-abs(tau_d[(K+1),])	

######################### Cycle definition

VaR_t<- rep(NA,N)
VaR_t[1:q] <- quantile(ret,tau)

 for(i in (q+1):N) {
      VaR_t[i] <- X_t[i,]%*%cbind(param[1:(q+1)])+
			    theta*tau_d[i]+
                   z*VaR_t[i-1]
    }

es<-(1+exp(param[q+5]))*VaR_t
es<- -abs(es)

res_fin<-cbind(VaR_t,es)
colnames(res_fin)<-c("VaR","ES")

return(res_fin)

}


##################################################### QR model via ALD with ES

###############################################################################################
# Computes the ALD log-likelihood for the QR model with Expected Shortfall
# Input  : param -> vector of model parameters
#          ret   -> return series
#          X_t   -> matrix of regressors
#          tau   -> quantile level
# Output : vector of ALD contributions
###############################################################################################


lARCH_fun_es <- function(param, ret, X_t, tau){

q<-ncol(X_t)-1

N<-length(ret)

VaR_t<- rep(NA,N)
VaR_t[1:q] <- quantile(ret,tau)

 for(i in (q+1):N) {
      VaR_t[i] <- X_t[i,]%*%cbind(param[1:(q+1)])
    }

es<-(1+exp(param[q+2]))*VaR_t
es<- -abs(es)
num_llk<-QL(ret,VaR_t,tau)
den_llk<-tau*es

llk<-(tau-1)*(es)^-1*exp(num_llk/den_llk)
llk<-llk[(q+1):N]
llk<-ifelse(is.na(llk)|llk<0|llk==-Inf|llk==0,0.000001,llk)
return(log(llk))

}

############################################### QR model VaR and ES

###############################################################################################
# Computes Value-at-Risk and ES for the QR model
# Input  : param -> vector of model parameters
#          ret   -> return series
#          X_t   -> matrix of regressors
#          tau   -> quantile level
# Output : matrix containing Value-at-Risk and Expected Shortfall forecasts
###############################################################################################

lARCH_var_es <- function(param, ret, X_t, tau){

q<-ncol(X_t)-1

N<-length(ret)

VaR_t<- rep(NA,N)
VaR_t[1:q] <- quantile(ret,tau)

 for(i in (q+1):N) {
      VaR_t[i] <- X_t[i,]%*%cbind(param[1:(q+1)])
    }

return(VaR_t)

}

####################################### Methods for estimating (and evaluating) a variety of quantile regression based models

###############################################################################################
# Estimates and evaluates QR and MF-QR models,
# optionally including MIDAS and exogenous variables
# Input  : model         -> model specification ("lARCH", "lARCHX", "lARCHMIDAS", or "lARCHMIDASX", respectively corresponding to QR, QR-X, MF-QR, and MF-QR-X)
#          tau           -> quantile level
#          daily_ret     -> return series
#          q             -> number of autoregressive lags
#          Exp_Short     -> whether Expected Shortfall is estimated
#          B             -> bootstrap replications
#          mv_m          -> MIDAS variable matrix (optional)
#          K             -> number of MIDAS lags (optional)
#          X             -> exogenous variable (optional)
#          out_of_sample -> out-of-sample size (optional)
# Output : object containing estimation and forecasting results
###############################################################################################

uqfit<-function(
model,
tau,
daily_ret,
q,
Exp_Short="No",
B=100,
mv_m=NULL,
K=NULL,
X=NULL,
out_of_sample=NULL,
...){

############################# check on valid choices
if(tau<0|tau>1) { stop(cat("#Warning:\n Parameter 'tau' must be a number strictly included between 0 and 1 \n"))}
if((model != "lARCH")&(model != "lARCHX")&(model != "lARCHMIDAS")&(model != "lARCHMIDASX")) { stop(cat("#Warning:\n Valid choices for the parameter 'model' 
are currently 'lARCH', 'lARCHX', 'lARCHMIDAS', and 'lARCHMIDASX' \n"))}
if(missing(q)) { stop(cat("#Warning:\n Parameter 'q' must be provided \n"))}
if(q<1) { stop(cat("#Warning:\n Parameter 'q' must be positive \n"))}

N<-length(daily_ret)

cond_r_t<- class(daily_ret)[1]
cond_mv_m<- class(mv_m)[1]

if(cond_r_t != "xts") { stop(
cat("#Warning:\n Parameter 'daily_ret' must be an xts object. Please provide it in the correct form \n")
)}

if(q > 10) { stop(
cat("#Warning:\n Parameter 'q' must be at most equal to 10 \n")
)}

if((model=="lARCHMIDAS"|model=="lARCHMIDASX")&cond_mv_m != "matrix") { stop(
cat("#Warning:\n Parameter 'mv_m' must be a matrix. Please provide it in the correct form \n")
)}

if((model=="lARCHX"|model=="lARCHMIDASX")& missing(X)) { stop(
cat("#Warning:\n Parameter 'X' must be provided. \n")
)}


if((model=="lARCHMIDAS"|model=="lARCHMIDASX")& missing(K)) { stop(
cat("#Warning:\n Parameter 'K' must be provided. \n")
)}


############## check if the X term is provided and if it has the same time span of daily_ret

if(!missing(X)){ 
if(any(range(time(daily_ret))!=range(time(X)))){
stop(
cat("#Warning:\n The vector 'X' has to be observed during the same time span of 'daily_ret' \n")
)}}

if (missing(out_of_sample)){ # out_of_sample parameter is missing

r_t_in_s<-daily_ret
mv_m_in_s<-mv_m
X_in_s<-X


} else {					# out_of_sample parameter is present

r_t_in_s<-daily_ret[1:(N-out_of_sample)]
mv_m_in_s<-mv_m[,1:(N-out_of_sample)]
X_in_s<-X[1:(N-out_of_sample)]

r_t_oos<-daily_ret[(N-out_of_sample+1):N]
mv_m_oos<-mv_m[,(N-out_of_sample+1):N]
X_oos<-X[(N-out_of_sample+1):N]

N<-length(r_t_in_s)

}



############################################### Estimation

if (model=="lARCH"&Exp_Short=="No"){
ES<-ES_oos<-NA
############################ construction of X_t matrix

X_t<-cbind(
rep(1,N),
stats::lag(abs(r_t_in_s),k=1:q)) 

X_t<-replace(X_t,is.na(X_t),0)

colnames(X_t)<-paste("beta","_",0:q,sep="")

## Est.

X_t<-zoo::coredata(cbind(X_t[,2:(q+1)]))
Y<-zoo::coredata(r_t_in_s)
MAT<-data.frame(Y,X_t)
colnames(MAT)[1]<-"Y_t"
res<-rq(formula = Y_t ~ ., tau = tau, data=MAT)

VaR<-stats::fitted.values(res)

if (!missing(out_of_sample)) {

N_full<-length(daily_ret)
X_t_full<-cbind(
stats::lag(abs(daily_ret),k=1:q)) 

X_t_full<-replace(X_t_full,is.na(X_t_full),0)
X_t_oos<-zoo::coredata(X_t_full[(N_full-out_of_sample*5-q+1):N_full,])
X_t_oos<-as.data.frame(X_t_oos)
colnames(X_t_oos)<-colnames(X_t)
VaR_oos<-stats::predict(res,newdata=X_t_oos)[(length(VaR_oos)-out_of_sample+1):length(VaR_oos)]

}

} else if (model=="lARCH"&Exp_Short=="Yes"){
############################ construction of X_t matrix

X_t<-cbind(
rep(1,N),
stats::lag(abs(r_t_in_s),k=1:q)) 

X_t<-replace(X_t,is.na(X_t),0)

colnames(X_t)<-paste("beta","_",0:q,sep="")

## Est.

X_t<-zoo::coredata(cbind(X_t[,2:(q+1)]))
Y<-zoo::coredata(r_t_in_s)
MAT<-data.frame(Y,X_t)
colnames(MAT)[1]<-"Y_t"
res<-rq(formula = Y_t ~ ., tau = tau, data=MAT)

start_val<-c(stats::coef(res),0)

names(start_val)[(q+2)]<-"theta_gamma"

res<- maxLik(lARCH_fun_es, 
start=start_val, 
ret=zoo::coredata(r_t_in_s),
X_t=cbind(1,X_t),
tau=tau,
method="BFGS")

VaR<-lARCH_var_es(stats::coef(res)[1:(q+1)],
zoo::coredata(r_t_in_s),cbind(1,X_t),tau=tau)
ES<-(1+exp(stats::coef(res)[q+2]))*VaR
ES<-as.xts(ES,stats::time(r_t_in_s))

################################ oos

if (!missing(out_of_sample)) {

N_full<-length(daily_ret)
X_t_full<-cbind(
stats::lag(abs(daily_ret),k=1:q)) 

X_t_full<-replace(X_t_full,is.na(X_t_full),0)
X_t_oos<-zoo::coredata(X_t_full[(N_full-out_of_sample*5-q+1):N_full,])
X_t_oos<-as.matrix(X_t_oos)
colnames(X_t_oos)<-colnames(X_t)

r_t_oos_2<-daily_ret[(N_full-out_of_sample*5-q+1):N_full]

VaR_oos<-lARCH_var_es(stats::coef(res)[1:(q+1)],
zoo::coredata(r_t_oos_2),cbind(1,X_t_oos),tau=tau)

VaR_oos<-VaR_oos[(length(VaR_oos)-out_of_sample+1):length(VaR_oos)]

ES_oos<-(1+exp(stats::coef(res)[q+2]))*VaR_oos
ES_oos<-as.xts(ES_oos,stats::time(r_t_oos))

}

} else if (model=="lARCHX"&Exp_Short=="No"){
ES<-ES_oos<-NA
############################ construction of X_t matrix

X_t<-cbind(
rep(1,N),
stats::lag(abs(r_t_in_s),k=1:q),
stats::lag(abs(X_in_s),k=1)) 

X_t<-replace(X_t,is.na(X_t),0)

colnames(X_t)[1:(q+1)]<-paste("beta","_",0:q,sep="")
colnames(X_t)[(q+2)]<-paste("beta","_","X",sep="")

## Est.

X_t<-zoo::coredata(cbind(X_t[,2:(q+2)]))
Y<-zoo::coredata(r_t_in_s)
MAT<-data.frame(Y,X_t)
colnames(MAT)[1]<-"Y_t"

res<-rq(formula = Y_t ~ ., tau = tau, data=MAT)

VaR<-stats::fitted.values(res)

################################ oos

if (!missing(out_of_sample)) {

N_full<-length(daily_ret)
X_t_full<-cbind(
stats::lag(abs(daily_ret),k=1:q),
stats::lag(abs(X),k=1)) 

X_t_full<-replace(X_t_full,is.na(X_t_full),0)
X_t_oos<-zoo::coredata(X_t_full[(N_full-out_of_sample*5-q+1):N_full,])
X_t_oos<-as.data.frame(X_t_oos)
colnames(X_t_oos)<-colnames(X_t)
VaR_oos<-stats::predict(res,newdata=X_t_oos)[(length(VaR_oos)-out_of_sample+1):length(VaR_oos)]

}


} else if (model=="lARCHX"&Exp_Short=="Yes"){
############################ construction of X_t matrix

X_t<-cbind(
rep(1,N),
stats::lag(abs(r_t_in_s),k=1:q),
stats::lag(abs(X_in_s),k=1)) 

X_t<-replace(X_t,is.na(X_t),0)

colnames(X_t)[1:(q+1)]<-paste("beta","_",0:q,sep="")
colnames(X_t)[(q+2)]<-paste("beta","_","X",sep="")

## Est.

X_t<-zoo::coredata(cbind(X_t[,2:(q+2)]))
Y<-zoo::coredata(r_t_in_s)
MAT<-data.frame(Y,X_t)
colnames(MAT)[1]<-"Y_t"
res<-rq(formula = Y_t ~ ., tau = tau, data=MAT)

start_val<-c(stats::coef(res),0)

names(start_val)[(q+3)]<-"theta_gamma"

if(start_val["beta_X"]<0.001){
start_val["beta_X"]<-0.01
}


res<- maxLik(lARCH_fun_es, 
start=start_val, 
ret=zoo::coredata(r_t_in_s),
X_t=cbind(1,X_t),
tau=tau,
method="BFGS")

VaR<-lARCH_var_es(stats::coef(res)[1:(q+2)],
zoo::coredata(r_t_in_s),cbind(1,X_t),tau=tau)
ES<-(1+exp(stats::coef(res)[q+3]))*VaR
ES<-as.xts(ES,stats::time(r_t_in_s))

################################ oos

if (!missing(out_of_sample)) {

N_full<-length(daily_ret)
X_t_full<-cbind(
stats::lag(abs(daily_ret),k=1:q),
stats::lag(abs(X),k=1)) 

X_t_full<-replace(X_t_full,is.na(X_t_full),0)
X_t_oos<-zoo::coredata(X_t_full[(N_full-out_of_sample*5-q+1):N_full,])
X_t_oos<-as.matrix(X_t_oos)
colnames(X_t_oos)<-colnames(X_t)

r_t_oos_2<-daily_ret[(N_full-out_of_sample*5-q+1):N_full]

VaR_oos<-lARCH_var_es(stats::coef(res)[1:(q+2)],
zoo::coredata(r_t_oos_2),cbind(1,X_t_oos),tau=tau)

VaR_oos<-VaR_oos[(length(VaR_oos)-out_of_sample+1):length(VaR_oos)]

ES_oos<-(1+exp(stats::coef(res)[q+3]))*VaR_oos
ES_oos<-as.xts(ES_oos,stats::time(r_t_oos))

}
######################################################
} else if (model=="lARCHMIDAS"&Exp_Short=="No"){
ES<-ES_oos<-NA

X_t<-cbind(
rep(1,N),
stats::lag(abs(r_t_in_s),k=1:q)) 

X_t<-replace(X_t,is.na(X_t),0)

colnames(X_t)<-paste("beta","_",0:q,sep="")

## MIDAS variable inclusion with profiling

res_first_step<-list()
rho_est<-list()

w1<-1
w2_seq<-seq(1.0,10,0.1)

for(r in 1:length(w2_seq)){

w2<-w2_seq[r]

betas<-c(rev(rumidas::beta_function(1:(K+1),(K+1),w1,w2))[2:(K+1)],0)
tau_d<- suppressWarnings(roll::roll_sum(mv_m_in_s, c(K+1),weights = betas))
tau_d<-tau_d[(K+1),]	 

X_t_2<-NULL
X_t_2<-zoo::coredata(cbind(X_t,abs(tau_d)))[,2:(q+2)]
colnames(X_t_2)[ncol(X_t_2)]<-c("theta")
Y<-coredata(r_t_in_s)
MAT<-data.frame(Y,X_t_2)
colnames(MAT)[1]<-"Y_t"

res_first_step[[r]]<-suppressWarnings(rq(formula = Y_t~ ., tau = tau,  data=MAT))
rho_est[[r]]<-res_first_step[[r]]$rho
}

rho_est<-unlist(rho_est)

################ find the optimal value of w2 and hence the optimal estimates

w2_star<-w2_seq[which.min(rho_est)]
r_star<-which.min(rho_est)

res<-(res_first_step[r_star])[[1]]

VaR<-stats::fitted.values(res)

################################ oos

if (!missing(out_of_sample)) {

N_full<-length(daily_ret)
X_t_full<-cbind(
stats::lag(abs(daily_ret),k=1:q)) 

X_t_full<-replace(X_t_full,is.na(X_t_full),0)
X_t_oos<-zoo::coredata(X_t_full[(N_full-out_of_sample*5-q+1):N_full,])

mv_m_oos_2<-mv_m[,(N_full-out_of_sample*5-q+1):N_full]

betas<-c(rev(rumidas::beta_function(1:(K+1),(K+1),w1,w2_star))[2:(K+1)],0)
tau_d<- suppressWarnings(roll::roll_sum(mv_m_oos_2, c(K+1),weights = betas))
tau_d<-tau_d[(K+1),]	 

X_t_oos<-as.data.frame(cbind(X_t_oos,abs(tau_d)))

colnames(X_t_oos)<-colnames(X_t_2)

VaR_oos<-stats::predict(res,newdata=X_t_oos)[(length(VaR_oos)-out_of_sample+1):length(VaR_oos)]

}
########################################################
} else if (model=="lARCHMIDAS"&Exp_Short=="Yes"){

X_t<-cbind(
rep(1,N),
stats::lag(abs(r_t_in_s),k=1:q)) 

X_t<-replace(X_t,is.na(X_t),0)

colnames(X_t)<-paste("beta","_",0:q,sep="")

## MIDAS variable inclusion with profiling

res_first_step<-list()
rho_est<-list()

w1<-1
w2_seq<-seq(1.0,10,0.1)

for(r in 1:length(w2_seq)){

w2<-w2_seq[r]

betas<-c(rev(rumidas::beta_function(1:(K+1),(K+1),w1,w2))[2:(K+1)],0)
tau_d<- suppressWarnings(roll::roll_sum(mv_m_in_s, c(K+1),weights = betas))
tau_d<-tau_d[(K+1),]	 

X_t_2<-NULL
X_t_2<-zoo::coredata(cbind(X_t,abs(tau_d)))[,2:(q+2)]
colnames(X_t_2)[ncol(X_t_2)]<-c("theta")
Y<-coredata(r_t_in_s)
MAT<-data.frame(Y,X_t_2)
colnames(MAT)[1]<-"Y_t"

res_first_step[[r]]<-suppressWarnings(rq(formula = Y_t~ ., tau = tau,  data=MAT))
rho_est[[r]]<-res_first_step[[r]]$rho
}

rho_est<-unlist(rho_est)

################ find the optimal value of w2 and hence the optimal estimates

w2_star<-w2_seq[which.min(rho_est)]
r_star<-which.min(rho_est)

res<-(res_first_step[r_star])[[1]]

############### Estimation via ALD

betas<-c(rev(rumidas::beta_function(1:(K+1),(K+1),w1,w2_star))[2:(K+1)],0)
tau_d<- suppressWarnings(roll::roll_sum(mv_m_in_s, c(K+1),weights = betas))
tau_d<-tau_d[(K+1),]	 

X_t_2<-NULL
X_t_2<-zoo::coredata(cbind(X_t,abs(tau_d)))[,2:(q+2)]
colnames(X_t_2)[ncol(X_t_2)]<-c("theta")
X_t<-X_t_2

start_val<-c(stats::coef(res),0)

names(start_val)[(q+3)]<-"theta_gamma"

res<- maxLik(lARCH_fun_es, 
start=start_val, 
ret=zoo::coredata(r_t_in_s),
X_t=cbind(1,X_t),
tau=tau,
method="BFGS")

VaR<-lARCH_var_es(stats::coef(res)[1:(q+2)],
zoo::coredata(r_t_in_s),cbind(1,X_t),tau=tau)
ES<-(1+exp(stats::coef(res)[q+3]))*VaR
ES<-as.xts(ES,stats::time(r_t_in_s))

################################ oos

if (!missing(out_of_sample)) {

N_full<-length(daily_ret)
X_t_full<-cbind(
stats::lag(abs(daily_ret),k=1:q)) 

X_t_full<-replace(X_t_full,is.na(X_t_full),0)
X_t_oos<-zoo::coredata(X_t_full[(N_full-out_of_sample*5-q+1):N_full,])

mv_m_oos_2<-mv_m[,(N_full-out_of_sample*5-q+1):N_full]

betas<-c(rev(rumidas::beta_function(1:(K+1),(K+1),w1,w2_star))[2:(K+1)],0)
tau_d<- suppressWarnings(roll::roll_sum(mv_m_oos_2, c(K+1),weights = betas))
tau_d<-tau_d[(K+1),]	 

X_t_oos<-as.data.frame(cbind(X_t_oos,abs(tau_d)))

colnames(X_t_oos)<-colnames(X_t_2)
X_t_oos<-as.matrix(X_t_oos)

r_t_oos_2<-daily_ret[(N_full-out_of_sample*5-q+1):N_full]

VaR_oos<-lARCH_var_es(stats::coef(res)[1:(q+2)],
zoo::coredata(r_t_oos_2),cbind(1,X_t_oos),tau=tau)

VaR_oos<-VaR_oos[(length(VaR_oos)-out_of_sample+1):length(VaR_oos)]

ES_oos<-(1+exp(stats::coef(res)[q+3]))*VaR_oos
ES_oos<-as.xts(ES_oos,stats::time(r_t_oos))

}
######################################################
} else if (model=="lARCHMIDASX"&Exp_Short=="No"){
ES<-ES_oos<-NA
############################ construction of X_t matrix

X_t<-cbind(
rep(1,N),
stats::lag(abs(r_t_in_s),k=1:q),
stats::lag(abs(X_in_s),k=1)) 

X_t<-replace(X_t,is.na(X_t),0)

colnames(X_t)[1:(q+1)]<-paste("beta","_",0:q,sep="")
colnames(X_t)[(q+2)]<-paste("beta","_","X",sep="")

## MIDAS variable inclusion with profiling

res_first_step<-list()
rho_est<-list()

w1<-1
w2_seq<-seq(1.0,10,0.1)

for(r in 1:length(w2_seq)){

w2<-w2_seq[r]

betas<-c(rev(rumidas::beta_function(1:(K+1),(K+1),w1,w2))[2:(K+1)],0)
tau_d<- suppressWarnings(roll::roll_sum(mv_m_in_s, c(K+1),weights = betas)) 
tau_d<-tau_d[(K+1),]	 

X_t_2<-NULL
X_t_2<-coredata(cbind(X_t,abs(tau_d)))[,2:(q+3)]
colnames(X_t_2)[ncol(X_t_2)]<-c("theta")
Y<-coredata(r_t_in_s)
MAT<-data.frame(Y,X_t_2)
colnames(MAT)[1]<-"Y_t"

res_first_step[[r]]<-suppressWarnings(rq(formula = Y_t~ ., tau = tau,  data=MAT))

rho_est[[r]]<-res_first_step[[r]]$rho
}

rho_est<-unlist(rho_est)

################ find the optimal value of w2 and hence the optimal estimates

w2_star<-w2_seq[which.min(rho_est)]
r_star<-which.min(rho_est)

res<-(res_first_step[r_star])[[1]]

VaR<-stats::fitted.values(res)

 
################################ oos

if (!missing(out_of_sample)) {

N_full<-length(daily_ret)
X_t_full<-cbind(
stats::lag(abs(daily_ret),k=1:q),
stats::lag(abs(X),k=1)) 

X_t_full<-replace(X_t_full,is.na(X_t_full),0)
X_t_oos<-zoo::coredata(X_t_full[(N_full-out_of_sample*5-q+1):N_full,])

mv_m_oos_2<-mv_m[,(N_full-out_of_sample*5-q+1):N_full]

betas<-c(rev(rumidas::beta_function(1:(K+1),(K+1),w1,w2_star))[2:(K+1)],0)
tau_d<- suppressWarnings(roll::roll_sum(mv_m_oos_2, c(K+1),weights = betas))
tau_d<-tau_d[(K+1),]	 

X_t_oos<-as.data.frame(cbind(X_t_oos,tau_d))

colnames(X_t_oos)<-colnames(X_t_2)

VaR_oos<-stats::predict(res,newdata=X_t_oos)[(length(VaR_oos)-out_of_sample+1):length(VaR_oos)]

}
##################################################
} else if (model=="lARCHMIDASX"&Exp_Short=="Yes"){

X_t<-cbind(
rep(1,N),
stats::lag(abs(r_t_in_s),k=1:q),
stats::lag(abs(X_in_s),k=1)) 

X_t<-replace(X_t,is.na(X_t),0)

colnames(X_t)[1:(q+1)]<-paste("beta","_",0:q,sep="")
colnames(X_t)[(q+2)]<-paste("beta","_","X",sep="")

## MIDAS variable inclusion with profiling

res_first_step<-list()
rho_est<-list()

w1<-1
w2_seq<-seq(1.0,10,0.1)

for(r in 1:length(w2_seq)){

w2<-w2_seq[r]

betas<-c(rev(rumidas::beta_function(1:(K+1),(K+1),w1,w2))[2:(K+1)],0)
tau_d<- suppressWarnings(roll::roll_sum(mv_m_in_s, c(K+1),weights = betas))
tau_d<-tau_d[(K+1),]	 

X_t_2<-NULL
X_t_2<-zoo::coredata(cbind(X_t,abs(tau_d)))[,2:(q+3)]
colnames(X_t_2)[ncol(X_t_2)]<-c("theta")
Y<-coredata(r_t_in_s)
MAT<-data.frame(Y,X_t_2)
colnames(MAT)[1]<-"Y_t"

res_first_step[[r]]<-suppressWarnings(rq(formula = Y_t~ ., tau = tau,  data=MAT))
rho_est[[r]]<-res_first_step[[r]]$rho
}

rho_est<-unlist(rho_est)

################ find the optimal value of w2 and hence the optimal estimates

w2_star<-w2_seq[which.min(rho_est)]
r_star<-which.min(rho_est)

res<-(res_first_step[r_star])[[1]]

############### Estimation via ALD

betas<-c(rev(rumidas::beta_function(1:(K+1),(K+1),w1,w2_star))[2:(K+1)],0)
tau_d<- suppressWarnings(roll::roll_sum(mv_m_in_s, c(K+1),weights = betas))
tau_d<-tau_d[(K+1),]	 

X_t_2<-NULL
X_t_2<-zoo::coredata(cbind(X_t,abs(tau_d)))[,2:(q+3)]
colnames(X_t_2)[ncol(X_t_2)]<-c("theta")
X_t<-X_t_2

start_val<-c(stats::coef(res),0)

names(start_val)[(q+4)]<-"theta_gamma"

#if(start_val["beta_X"]<0.001){
#start_val["beta_X"]<-0.01
#}

res<- maxLik(lARCH_fun_es, 
start=start_val, 
ret=zoo::coredata(r_t_in_s),
X_t=cbind(1,X_t),
tau=tau,
method="BFGS")

VaR<-lARCH_var_es(stats::coef(res)[1:(q+3)],
zoo::coredata(r_t_in_s),cbind(1,X_t),tau=tau)
ES<-(1+exp(stats::coef(res)[q+4]))*VaR
ES<-as.xts(ES,stats::time(r_t_in_s))

################################ oos

if (!missing(out_of_sample)) {

N_full<-length(daily_ret)
X_t_full<-cbind(
stats::lag(abs(daily_ret),k=1:q),
stats::lag(abs(X),k=1)) 

X_t_full<-replace(X_t_full,is.na(X_t_full),0)
X_t_oos<-zoo::coredata(X_t_full[(N_full-out_of_sample*5-q+1):N_full,])

mv_m_oos_2<-mv_m[,(N_full-out_of_sample*5-q+1):N_full]

betas<-c(rev(rumidas::beta_function(1:(K+1),(K+1),w1,w2_star))[2:(K+1)],0)
tau_d<- suppressWarnings(roll::roll_sum(mv_m_oos_2, c(K+1),weights = betas))
tau_d<-tau_d[(K+1),]	 

X_t_oos<-as.data.frame(cbind(X_t_oos,abs(tau_d)))

r_t_oos_2<-daily_ret[(N_full-out_of_sample*5-q+1):N_full]

colnames(X_t_oos)<-colnames(X_t_2)
X_t_oos<-as.matrix(X_t_oos)

VaR_oos<-lARCH_var_es(stats::coef(res)[1:(q+3)],
zoo::coredata(r_t_oos_2),cbind(1,X_t_oos),tau=tau)

VaR_oos<-VaR_oos[(length(VaR_oos)-out_of_sample+1):length(VaR_oos)]

ES_oos<-(1+exp(stats::coef(res)[q+4]))*VaR_oos
ES_oos<-as.xts(ES_oos,stats::time(r_t_oos))

}
######################################################
}

###################################################### 
###################################################### standard errors
######################################################
 
N_coef<-length(stats::coef(res))

# standardized residuals

stand_resids<-Y/abs(VaR)
N_est<-length(stand_resids)

mat_coef<-data.frame(rep(NA,N_coef),rep(NA,N_coef),rep(NA,N_coef),rep(NA,N_coef))
colnames(mat_coef)<-c("Estimate","Std. Error","t value","Pr(>|t|)")
rownames(mat_coef)<-names(stats::coef(res))

mat_boot_results<-matrix(rep(NA),nrow=B,ncol=N_coef)

pb <- utils::txtProgressBar(0, B, style = 3)

for(b in 1:B){

##### first step: compute the boostrapped series of residuals

resids_boot<-sample(stand_resids,N_est,replace=T)

##### second step: compute the bootstraped series of returns and VaRs

r_t_boot<-Y
VaR_boot<-VaR

if(model=="lARCH"){ ############################################# model Linear ARCH

if(q==1){
b0<-stats::coef(res)[1]
b1<-stats::coef(res)[2]

for(i in (q+1):N_est){
VaR_boot[i]<-b0 + b1*abs(r_t_boot[i-1]) 
r_t_boot[i]<-abs(VaR_boot[i])*resids_boot[i]
}
} else if(q==2){

b0<-stats::coef(res)[1]
b1<-stats::coef(res)[2]
b2<-stats::coef(res)[3]

for(i in (q+1):N_est){
VaR_boot[i]<-b0 + b1*abs(r_t_boot[i-1])  + b2*abs(r_t_boot[i-2])
r_t_boot[i]<-abs(VaR_boot[i])*resids_boot[i]
}
} else if(q==3){

b0<-stats::coef(res)[1]
b1<-stats::coef(res)[2]
b2<-stats::coef(res)[3]
b3<-stats::coef(res)[4]

for(i in (q+1):N_est){
VaR_boot[i]<-b0 + b1*abs(r_t_boot[i-1])  + b2*abs(r_t_boot[i-2])+ b3*abs(r_t_boot[i-3])
r_t_boot[i]<-abs(VaR_boot[i])*resids_boot[i]
}
} else if(q==4){

b0<-stats::coef(res)[1]
b1<-stats::coef(res)[2]
b2<-stats::coef(res)[3]
b3<-stats::coef(res)[4]
b4<-stats::coef(res)[5]

for(i in (q+1):N_est){
VaR_boot[i]<-b0 + b1*abs(r_t_boot[i-1])  + b2*abs(r_t_boot[i-2])+ b3*abs(r_t_boot[i-3])+
b4*abs(r_t_boot[i-4])
r_t_boot[i]<-abs(VaR_boot[i])*resids_boot[i]
}
} else if(q==5){

b0<-stats::coef(res)[1]
b1<-stats::coef(res)[2]
b2<-stats::coef(res)[3]
b3<-stats::coef(res)[4]
b4<-stats::coef(res)[5]
b5<-stats::coef(res)[6]

for(i in (q+1):N_est){
VaR_boot[i]<-b0 + b1*abs(r_t_boot[i-1])  + b2*abs(r_t_boot[i-2])+ b3*abs(r_t_boot[i-3])+
b4*abs(r_t_boot[i-4])+b5*abs(r_t_boot[i-5])
r_t_boot[i]<-abs(VaR_boot[i])*resids_boot[i]
}
} else if(q==6){

b0<-stats::coef(res)[1]
b1<-stats::coef(res)[2]
b2<-stats::coef(res)[3]
b3<-stats::coef(res)[4]
b4<-stats::coef(res)[5]
b5<-stats::coef(res)[6]
b6<-stats::coef(res)[7]

for(i in (q+1):N_est){
VaR_boot[i]<-b0 + b1*abs(r_t_boot[i-1])  + b2*abs(r_t_boot[i-2])+ b3*abs(r_t_boot[i-3])+
b4*abs(r_t_boot[i-4])+b5*abs(r_t_boot[i-5])+b6*abs(r_t_boot[i-6])
r_t_boot[i]<-abs(VaR_boot[i])*resids_boot[i]
}
} else if(q==7){

b0<-stats::coef(res)[1]
b1<-stats::coef(res)[2]
b2<-stats::coef(res)[3]
b3<-stats::coef(res)[4]
b4<-stats::coef(res)[5]
b5<-stats::coef(res)[6]
b6<-stats::coef(res)[7]
b7<-stats::coef(res)[8]

for(i in (q+1):N_est){
VaR_boot[i]<-b0 + b1*abs(r_t_boot[i-1]) + b2*abs(r_t_boot[i-2])+ b3*abs(r_t_boot[i-3])+
b4*abs(r_t_boot[i-4])+b5*abs(r_t_boot[i-5])+b6*abs(r_t_boot[i-6])+b7*abs(r_t_boot[i-7])
r_t_boot[i]<-abs(VaR_boot[i])*resids_boot[i]
}
} else if(q==8){

b0<-stats::coef(res)[1]
b1<-stats::coef(res)[2]
b2<-stats::coef(res)[3]
b3<-stats::coef(res)[4]
b4<-stats::coef(res)[5]
b5<-stats::coef(res)[6]
b6<-stats::coef(res)[7]
b7<-stats::coef(res)[8]
b8<-stats::coef(res)[9]

for(i in (q+1):N_est){
VaR_boot[i]<-b0 + b1*abs(r_t_boot[i-1]) + b2*abs(r_t_boot[i-2])+ b3*abs(r_t_boot[i-3])+
b4*abs(r_t_boot[i-4])+b5*abs(r_t_boot[i-5])+b6*abs(r_t_boot[i-6])+b7*abs(r_t_boot[i-7])+
b8*abs(r_t_boot[i-8])
r_t_boot[i]<-abs(VaR_boot[i])*resids_boot[i]
} 
} else if(q==9){

b0<-stats::coef(res)[1]
b1<-stats::coef(res)[2]
b2<-stats::coef(res)[3]
b3<-stats::coef(res)[4]
b4<-stats::coef(res)[5]
b5<-stats::coef(res)[6]
b6<-stats::coef(res)[7]
b7<-stats::coef(res)[8]
b8<-stats::coef(res)[9]
b9<-stats::coef(res)[10]

for(i in (q+1):N_est){
VaR_boot[i]<-b0 + b1*abs(r_t_boot[i-1]) + b2*abs(r_t_boot[i-2])+ b3*abs(r_t_boot[i-3])+
b4*abs(r_t_boot[i-4])+b5*abs(r_t_boot[i-5])+b6*abs(r_t_boot[i-6])+b7*abs(r_t_boot[i-7])+
b8*abs(r_t_boot[i-8])+b9*abs(r_t_boot[i-9])
r_t_boot[i]<-abs(VaR_boot[i])*resids_boot[i]
}

} else if(q==10){

b0<-stats::coef(res)[1]
b1<-stats::coef(res)[2]
b2<-stats::coef(res)[3]
b3<-stats::coef(res)[4]
b4<-stats::coef(res)[5]
b5<-stats::coef(res)[6]
b6<-stats::coef(res)[7]
b7<-stats::coef(res)[8]
b8<-stats::coef(res)[9]
b9<-stats::coef(res)[10]
b10<-stats::coef(res)[11]

for(i in (q+1):N_est){
VaR_boot[i]<-b0 + b1*abs(r_t_boot[i-1])  + b2*abs(r_t_boot[i-2])+ b3*abs(r_t_boot[i-3])+
b4*abs(r_t_boot[i-4])+b5*abs(r_t_boot[i-5])+b6*abs(r_t_boot[i-6])+b7*abs(r_t_boot[i-7])+
b8*abs(r_t_boot[i-8])+b9*abs(r_t_boot[i-9])+b10*abs(r_t_boot[i-10])
r_t_boot[i]<-abs(VaR_boot[i])*resids_boot[i]
}

} 

######################################### est

r_t_boot<-as.xts(r_t_boot,time(r_t_in_s))

X_t_boot<-cbind(
rep(1,N_est),
stats::lag(abs(r_t_boot),k=1:q)) 

X_t_boot<-replace(X_t_boot,is.na(X_t_boot),0)

colnames(X_t_boot)<-paste("beta","_",0:q,sep="")

## Estimation

X_t_boot<-zoo::coredata(cbind(X_t_boot[,2:(q+1)]))
Y_boot<-zoo::coredata(r_t_boot)
MAT_boot<-data.frame(Y_boot,X_t_boot)
colnames(MAT_boot)[1]<-"Y_t"

if(Exp_Short=="No"){
res_boot<-rq(formula = Y_t ~ ., tau = tau, data=MAT_boot)
} else if (Exp_Short=="Yes") {

res_boot_start<-rq(formula = Y_t ~ ., tau = tau, data=MAT_boot)
start_val_boot<-c(stats::coef(res_boot_start),0)
names(start_val_boot)[(q+2)]<-"theta_gamma"


res_boot<- maxLik(lARCH_fun_es, 
start=start_val_boot, 
ret=Y_boot,
X_t=cbind(1,X_t_boot),
tau=tau,
method="BFGS")

}
############################################ end Linear ARCH
} else if(model=="lARCHX"){ ############################################# model Linear ARCH-X

X_boot<-abs(X_in_s)

if(q==1){
b0<-stats::coef(res)[1]
b1<-stats::coef(res)[2]
z<-stats::coef(res)[3]

for(i in (q+1):N_est){
VaR_boot[i]<-b0 + b1*abs(r_t_boot[i-1]) + z*X_boot[i-1]
r_t_boot[i]<-abs(VaR_boot[i])*resids_boot[i]
}
} else if(q==2){

b0<-stats::coef(res)[1]
b1<-stats::coef(res)[2]
b2<-stats::coef(res)[3]
z<-stats::coef(res)[4]

for(i in (q+1):N_est){
VaR_boot[i]<-b0 + b1*abs(r_t_boot[i-1])  + b2*abs(r_t_boot[i-2]) + z*X_boot[i-1]
r_t_boot[i]<-abs(VaR_boot[i])*resids_boot[i]
}
} else if(q==3){

b0<-stats::coef(res)[1]
b1<-stats::coef(res)[2]
b2<-stats::coef(res)[3]
b3<-stats::coef(res)[4]
z<-stats::coef(res)[5]

for(i in (q+1):N_est){
VaR_boot[i]<-b0 + b1*abs(r_t_boot[i-1])  + b2*abs(r_t_boot[i-2])+ b3*abs(r_t_boot[i-3]) + z*X_boot[i-1]
r_t_boot[i]<-abs(VaR_boot[i])*resids_boot[i]
}
} else if(q==4){

b0<-stats::coef(res)[1]
b1<-stats::coef(res)[2]
b2<-stats::coef(res)[3]
b3<-stats::coef(res)[4]
b4<-stats::coef(res)[5]
z<-stats::coef(res)[6]

for(i in (q+1):N_est){
VaR_boot[i]<-b0 + b1*abs(r_t_boot[i-1])  + b2*abs(r_t_boot[i-2])+ b3*abs(r_t_boot[i-3])+
b4*abs(r_t_boot[i-4]) + z*X_boot[i-1]
r_t_boot[i]<-abs(VaR_boot[i])*resids_boot[i]
}
} else if(q==5){

b0<-stats::coef(res)[1]
b1<-stats::coef(res)[2]
b2<-stats::coef(res)[3]
b3<-stats::coef(res)[4]
b4<-stats::coef(res)[5]
b5<-stats::coef(res)[6]
z<-stats::coef(res)[7]

for(i in (q+1):N_est){
VaR_boot[i]<-b0 + b1*abs(r_t_boot[i-1])  + b2*abs(r_t_boot[i-2])+ b3*abs(r_t_boot[i-3])+
b4*abs(r_t_boot[i-4])+b5*abs(r_t_boot[i-5])+ z*X_boot[i-1]
r_t_boot[i]<-abs(VaR_boot[i])*resids_boot[i]
}
} else if(q==6){

b0<-stats::coef(res)[1]
b1<-stats::coef(res)[2]
b2<-stats::coef(res)[3]
b3<-stats::coef(res)[4]
b4<-stats::coef(res)[5]
b5<-stats::coef(res)[6]
b6<-stats::coef(res)[7]
z<-stats::coef(res)[8]

for(i in (q+1):N_est){
VaR_boot[i]<-b0 + b1*abs(r_t_boot[i-1])  + b2*abs(r_t_boot[i-2])+ b3*abs(r_t_boot[i-3])+
b4*abs(r_t_boot[i-4])+b5*abs(r_t_boot[i-5])+b6*abs(r_t_boot[i-6])+ z*X_boot[i-1]
r_t_boot[i]<-abs(VaR_boot[i])*resids_boot[i]
}
} else if(q==7){

b0<-stats::coef(res)[1]
b1<-stats::coef(res)[2]
b2<-stats::coef(res)[3]
b3<-stats::coef(res)[4]
b4<-stats::coef(res)[5]
b5<-stats::coef(res)[6]
b6<-stats::coef(res)[7]
b7<-stats::coef(res)[8]
z<-stats::coef(res)[9]

for(i in (q+1):N_est){
VaR_boot[i]<-b0 + b1*abs(r_t_boot[i-1]) + b2*abs(r_t_boot[i-2])+ b3*abs(r_t_boot[i-3])+
b4*abs(r_t_boot[i-4])+b5*abs(r_t_boot[i-5])+b6*abs(r_t_boot[i-6])+b7*abs(r_t_boot[i-7])+ z*X_boot[i-1]
r_t_boot[i]<-abs(VaR_boot[i])*resids_boot[i]
}
} else if(q==8){

b0<-stats::coef(res)[1]
b1<-stats::coef(res)[2]
b2<-stats::coef(res)[3]
b3<-stats::coef(res)[4]
b4<-stats::coef(res)[5]
b5<-stats::coef(res)[6]
b6<-stats::coef(res)[7]
b7<-stats::coef(res)[8]
b8<-stats::coef(res)[9]
z<-stats::coef(res)[10]

for(i in (q+1):N_est){
VaR_boot[i]<-b0 + b1*abs(r_t_boot[i-1]) + b2*abs(r_t_boot[i-2])+ b3*abs(r_t_boot[i-3])+
b4*abs(r_t_boot[i-4])+b5*abs(r_t_boot[i-5])+b6*abs(r_t_boot[i-6])+b7*abs(r_t_boot[i-7])+
b8*abs(r_t_boot[i-8])+ z*X_boot[i-1]
r_t_boot[i]<-abs(VaR_boot[i])*resids_boot[i]
} 
} else if(q==9){

b0<-stats::coef(res)[1]
b1<-stats::coef(res)[2]
b2<-stats::coef(res)[3]
b3<-stats::coef(res)[4]
b4<-stats::coef(res)[5]
b5<-stats::coef(res)[6]
b6<-stats::coef(res)[7]
b7<-stats::coef(res)[8]
b8<-stats::coef(res)[9]
b9<-stats::coef(res)[10]
z<-stats::coef(res)[11]

for(i in (q+1):N_est){
VaR_boot[i]<-b0 + b1*abs(r_t_boot[i-1]) + b2*abs(r_t_boot[i-2])+ b3*abs(r_t_boot[i-3])+
b4*abs(r_t_boot[i-4])+b5*abs(r_t_boot[i-5])+b6*abs(r_t_boot[i-6])+b7*abs(r_t_boot[i-7])+
b8*abs(r_t_boot[i-8])+b9*abs(r_t_boot[i-9])+ z*X_boot[i-1]
r_t_boot[i]<-abs(VaR_boot[i])*resids_boot[i]
}

} else if(q==10){

b0<-stats::coef(res)[1]
b1<-stats::coef(res)[2]
b2<-stats::coef(res)[3]
b3<-stats::coef(res)[4]
b4<-stats::coef(res)[5]
b5<-stats::coef(res)[6]
b6<-stats::coef(res)[7]
b7<-stats::coef(res)[8]
b8<-stats::coef(res)[9]
b9<-stats::coef(res)[10]
b10<-stats::coef(res)[11]
z<-stats::coef(res)[12]

for(i in (q+1):N_est){
VaR_boot[i]<-b0 + b1*abs(r_t_boot[i-1])  + b2*abs(r_t_boot[i-2])+ b3*abs(r_t_boot[i-3])+
b4*abs(r_t_boot[i-4])+b5*abs(r_t_boot[i-5])+b6*abs(r_t_boot[i-6])+b7*abs(r_t_boot[i-7])+
b8*abs(r_t_boot[i-8])+b9*abs(r_t_boot[i-9])+b10*abs(r_t_boot[i-10])+ z*X_boot[i-1]
r_t_boot[i]<-abs(VaR_boot[i])*resids_boot[i]
}

} 

######################################### est

r_t_boot<-as.xts(r_t_boot,time(r_t_in_s))

X_t_boot<-cbind(
rep(1,N_est),
stats::lag(abs(r_t_boot),k=1:q),
stats::lag(abs(X_in_s),k=1)) 

X_t_boot<-replace(X_t_boot,is.na(X_t_boot),0)

colnames(X_t_boot)[1:(q+1)]<-paste("beta","_",0:q,sep="")
colnames(X_t_boot)[(q+2)]<-paste("beta","_","X",sep="")

## Estimation

X_t_boot<-zoo::coredata(cbind(X_t_boot[,2:(q+2)]))
Y_boot<-zoo::coredata(r_t_boot)
MAT_boot<-data.frame(Y_boot,X_t_boot)
colnames(MAT_boot)[1]<-"Y_t"

if(Exp_Short=="No"){
res_boot<-rq(formula = Y_t ~ ., tau = tau, data=MAT_boot)
} else if (Exp_Short=="Yes") {

res_boot_start<-rq(formula = Y_t ~ ., tau = tau, data=MAT_boot)
start_val_boot<-c(stats::coef(res_boot_start),0)
names(start_val_boot)[(q+3)]<-"theta_gamma"

if(start_val_boot["beta_X"]<0.001){
start_val_boot["beta_X"]<-0.01
}

res_boot<- maxLik(lARCH_fun_es, 
start=start_val_boot, 
ret=Y_boot,
X_t=cbind(1,X_t_boot),
tau=tau,
method="BFGS")

}
###### end Linear ARCH-X 
} else if(model=="lARCHMIDAS"){ ##################################### model Linear ARCH MIDAS

betas<-c(rev(rumidas::beta_function(1:(K+1),(K+1),w1,w2_star))[2:(K+1)],0)
tau_d<- suppressWarnings(roll::roll_sum(mv_m_in_s, c(K+1),weights = betas))
tau_d<-abs(tau_d[(K+1),])	

if(q==1){
b0<-stats::coef(res)[1]
b1<-stats::coef(res)[2]
theta<-stats::coef(res)[3]

for(i in (q+1):N_est){
VaR_boot[i]<-b0 + b1*abs(r_t_boot[i-1]) + theta*tau_d[i]
r_t_boot[i]<-abs(VaR_boot[i])*resids_boot[i]
}
} else if(q==2){

b0<-stats::coef(res)[1]
b1<-stats::coef(res)[2]
b2<-stats::coef(res)[3]
theta<-stats::coef(res)[4]

for(i in (q+1):N_est){
VaR_boot[i]<-b0 + b1*abs(r_t_boot[i-1])  + b2*abs(r_t_boot[i-2])
+ theta*tau_d[i]
r_t_boot[i]<-abs(VaR_boot[i])*resids_boot[i]
}
} else if(q==3){

b0<-stats::coef(res)[1]
b1<-stats::coef(res)[2]
b2<-stats::coef(res)[3]
b3<-stats::coef(res)[4]
theta<-stats::coef(res)[5]

for(i in (q+1):N_est){
VaR_boot[i]<-b0 + b1*abs(r_t_boot[i-1])  + b2*abs(r_t_boot[i-2])+ 
b3*abs(r_t_boot[i-3])+ theta*tau_d[i]
r_t_boot[i]<-abs(VaR_boot[i])*resids_boot[i]
}
} else if(q==4){

b0<-stats::coef(res)[1]
b1<-stats::coef(res)[2]
b2<-stats::coef(res)[3]
b3<-stats::coef(res)[4]
b4<-stats::coef(res)[5]
theta<-stats::coef(res)[6]

for(i in (q+1):N_est){
VaR_boot[i]<-b0 + b1*abs(r_t_boot[i-1])  + b2*abs(r_t_boot[i-2])+ b3*abs(r_t_boot[i-3])+
b4*abs(r_t_boot[i-4])+ theta*tau_d[i]
r_t_boot[i]<-abs(VaR_boot[i])*resids_boot[i]
}
} else if(q==5){

b0<-stats::coef(res)[1]
b1<-stats::coef(res)[2]
b2<-stats::coef(res)[3]
b3<-stats::coef(res)[4]
b4<-stats::coef(res)[5]
b5<-stats::coef(res)[6]
theta<-stats::coef(res)[7]

for(i in (q+1):N_est){
VaR_boot[i]<-b0 + b1*abs(r_t_boot[i-1])  + b2*abs(r_t_boot[i-2])+ b3*abs(r_t_boot[i-3])+
b4*abs(r_t_boot[i-4])+b5*abs(r_t_boot[i-5])+ theta*tau_d[i]
r_t_boot[i]<-abs(VaR_boot[i])*resids_boot[i]
}
} else if(q==6){

b0<-stats::coef(res)[1]
b1<-stats::coef(res)[2]
b2<-stats::coef(res)[3]
b3<-stats::coef(res)[4]
b4<-stats::coef(res)[5]
b5<-stats::coef(res)[6]
b6<-stats::coef(res)[7]
theta<-stats::coef(res)[8]

for(i in (q+1):N_est){
VaR_boot[i]<-b0 + b1*abs(r_t_boot[i-1])  + b2*abs(r_t_boot[i-2])+ b3*abs(r_t_boot[i-3])+
b4*abs(r_t_boot[i-4])+b5*abs(r_t_boot[i-5])+b6*abs(r_t_boot[i-6])+ theta*tau_d[i]
r_t_boot[i]<-abs(VaR_boot[i])*resids_boot[i]
}
} else if(q==7){

b0<-stats::coef(res)[1]
b1<-stats::coef(res)[2]
b2<-stats::coef(res)[3]
b3<-stats::coef(res)[4]
b4<-stats::coef(res)[5]
b5<-stats::coef(res)[6]
b6<-stats::coef(res)[7]
b7<-stats::coef(res)[8]
theta<-stats::coef(res)[9]

for(i in (q+1):N_est){
VaR_boot[i]<-b0 + b1*abs(r_t_boot[i-1]) + b2*abs(r_t_boot[i-2])+ b3*abs(r_t_boot[i-3])+
b4*abs(r_t_boot[i-4])+b5*abs(r_t_boot[i-5])+b6*abs(r_t_boot[i-6])+
b7*abs(r_t_boot[i-7])+ theta*tau_d[i]
r_t_boot[i]<-abs(VaR_boot[i])*resids_boot[i]
}
} else if(q==8){

b0<-stats::coef(res)[1]
b1<-stats::coef(res)[2]
b2<-stats::coef(res)[3]
b3<-stats::coef(res)[4]
b4<-stats::coef(res)[5]
b5<-stats::coef(res)[6]
b6<-stats::coef(res)[7]
b7<-stats::coef(res)[8]
b8<-stats::coef(res)[9]
theta<-stats::coef(res)[10]

for(i in (q+1):N_est){
VaR_boot[i]<-b0 + b1*abs(r_t_boot[i-1]) + b2*abs(r_t_boot[i-2])+ b3*abs(r_t_boot[i-3])+
b4*abs(r_t_boot[i-4])+b5*abs(r_t_boot[i-5])+b6*abs(r_t_boot[i-6])+
b7*abs(r_t_boot[i-7])+
b8*abs(r_t_boot[i-8])+ theta*tau_d[i]
r_t_boot[i]<-abs(VaR_boot[i])*resids_boot[i]
} 
} else if(q==9){

b0<-stats::coef(res)[1]
b1<-stats::coef(res)[2]
b2<-stats::coef(res)[3]
b3<-stats::coef(res)[4]
b4<-stats::coef(res)[5]
b5<-stats::coef(res)[6]
b6<-stats::coef(res)[7]
b7<-stats::coef(res)[8]
b8<-stats::coef(res)[9]
b9<-stats::coef(res)[10]
theta<-stats::coef(res)[11]

for(i in (q+1):N_est){
VaR_boot[i]<-b0 + b1*abs(r_t_boot[i-1]) + b2*abs(r_t_boot[i-2])+ b3*abs(r_t_boot[i-3])+
b4*abs(r_t_boot[i-4])+b5*abs(r_t_boot[i-5])+b6*abs(r_t_boot[i-6])+b7*abs(r_t_boot[i-7])+
b8*abs(r_t_boot[i-8])+b9*abs(r_t_boot[i-9])+ theta*tau_d[i]
r_t_boot[i]<-abs(VaR_boot[i])*resids_boot[i]
}

} else if(q==10){

b0<-stats::coef(res)[1]
b1<-stats::coef(res)[2]
b2<-stats::coef(res)[3]
b3<-stats::coef(res)[4]
b4<-stats::coef(res)[5]
b5<-stats::coef(res)[6]
b6<-stats::coef(res)[7]
b7<-stats::coef(res)[8]
b8<-stats::coef(res)[9]
b9<-stats::coef(res)[10]
b10<-stats::coef(res)[11]
theta<-stats::coef(res)[12]

for(i in (q+1):N_est){
VaR_boot[i]<-b0 + b1*abs(r_t_boot[i-1])  + b2*abs(r_t_boot[i-2])+ b3*abs(r_t_boot[i-3])+
b4*abs(r_t_boot[i-4])+b5*abs(r_t_boot[i-5])+b6*abs(r_t_boot[i-6])+b7*abs(r_t_boot[i-7])+
b8*abs(r_t_boot[i-8])+b9*abs(r_t_boot[i-9])+b10*abs(r_t_boot[i-10])+ theta*tau_d[i]
r_t_boot[i]<-abs(VaR_boot[i])*resids_boot[i]
}

} 

######################################### est

r_t_boot<-as.xts(r_t_boot,time(r_t_in_s))

X_t_boot<-cbind(
rep(1,N_est),
stats::lag(abs(r_t_boot),k=1:q),
tau_d) 

X_t_boot<-replace(X_t_boot,is.na(X_t_boot),0)

colnames(X_t_boot)[1:(q+1)]<-paste("beta","_",0:q,sep="")
colnames(X_t_boot)[ncol(X_t_boot)]<-c("theta")

## Estimation

X_t_boot<-zoo::coredata(cbind(X_t_boot[,2:(q+2)]))
Y_boot<-zoo::coredata(r_t_boot)
MAT_boot<-data.frame(Y_boot,X_t_boot)
colnames(MAT_boot)[1]<-"Y_t"

if(Exp_Short=="No"){
res_boot<-rq(formula = Y_t ~ ., tau = tau, data=MAT_boot)
} else if (Exp_Short=="Yes") {

res_boot_start<-rq(formula = Y_t ~ ., tau = tau, data=MAT_boot)
start_val_boot<-c(stats::coef(res_boot_start),0)
names(start_val_boot)[(q+3)]<-"theta_gamma"


res_boot<- maxLik(lARCH_fun_es, 
start=start_val_boot, 
ret=Y_boot,
X_t=cbind(1,X_t_boot),
tau=tau,
method="BFGS")

}
############################################ end Linear ARCH MIDAS
} else if(model=="lARCHMIDASX"){ ##################################### model Linear ARCH MIDAS X

betas<-c(rev(rumidas::beta_function(1:(K+1),(K+1),w1,w2_star))[2:(K+1)],0)
tau_d<- suppressWarnings(roll::roll_sum(mv_m_in_s, c(K+1),weights = betas))
tau_d<-abs(tau_d[(K+1),])	

X_boot<-abs(X_in_s)

if(q==1){
b0<-stats::coef(res)[1]
b1<-stats::coef(res)[2]
z<-stats::coef(res)[3]
theta<-stats::coef(res)[4]

for(i in (q+1):N_est){
VaR_boot[i]<-b0 + b1*abs(r_t_boot[i-1]) + theta*tau_d[i] + z*X_boot[i-1]
r_t_boot[i]<-abs(VaR_boot[i])*resids_boot[i]
}
} else if(q==2){

b0<-stats::coef(res)[1]
b1<-stats::coef(res)[2]
b2<-stats::coef(res)[3]
z<-stats::coef(res)[4]
theta<-stats::coef(res)[5]


for(i in (q+1):N_est){
VaR_boot[i]<-b0 + b1*abs(r_t_boot[i-1])  + b2*abs(r_t_boot[i-2])
+ theta*tau_d[i] + z*X_boot[i-1]
r_t_boot[i]<-abs(VaR_boot[i])*resids_boot[i]
}
} else if(q==3){

b0<-stats::coef(res)[1]
b1<-stats::coef(res)[2]
b2<-stats::coef(res)[3]
b3<-stats::coef(res)[4]
z<-stats::coef(res)[5]
theta<-stats::coef(res)[6]


for(i in (q+1):N_est){
VaR_boot[i]<-b0 + b1*abs(r_t_boot[i-1])  + b2*abs(r_t_boot[i-2])+ 
b3*abs(r_t_boot[i-3])+ theta*tau_d[i] + z*X_boot[i-1]
r_t_boot[i]<-abs(VaR_boot[i])*resids_boot[i]
}
} else if(q==4){

b0<-stats::coef(res)[1]
b1<-stats::coef(res)[2]
b2<-stats::coef(res)[3]
b3<-stats::coef(res)[4]
b4<-stats::coef(res)[5]
z<-stats::coef(res)[6]
theta<-stats::coef(res)[7]


for(i in (q+1):N_est){
VaR_boot[i]<-b0 + b1*abs(r_t_boot[i-1])  + b2*abs(r_t_boot[i-2])+ b3*abs(r_t_boot[i-3])+
b4*abs(r_t_boot[i-4])+ theta*tau_d[i] + z*X_boot[i-1]
r_t_boot[i]<-abs(VaR_boot[i])*resids_boot[i]
}
} else if(q==5){

b0<-stats::coef(res)[1]
b1<-stats::coef(res)[2]
b2<-stats::coef(res)[3]
b3<-stats::coef(res)[4]
b4<-stats::coef(res)[5]
b5<-stats::coef(res)[6]
z<-stats::coef(res)[7]
theta<-stats::coef(res)[8]


for(i in (q+1):N_est){
VaR_boot[i]<-b0 + b1*abs(r_t_boot[i-1])  + b2*abs(r_t_boot[i-2])+ b3*abs(r_t_boot[i-3])+
b4*abs(r_t_boot[i-4])+b5*abs(r_t_boot[i-5])+ theta*tau_d[i] + z*X_boot[i-1]
r_t_boot[i]<-abs(VaR_boot[i])*resids_boot[i]
}
} else if(q==6){

b0<-stats::coef(res)[1]
b1<-stats::coef(res)[2]
b2<-stats::coef(res)[3]
b3<-stats::coef(res)[4]
b4<-stats::coef(res)[5]
b5<-stats::coef(res)[6]
b6<-stats::coef(res)[7]
z<-stats::coef(res)[8]
theta<-stats::coef(res)[9]


for(i in (q+1):N_est){
VaR_boot[i]<-b0 + b1*abs(r_t_boot[i-1])  + b2*abs(r_t_boot[i-2])+ b3*abs(r_t_boot[i-3])+
b4*abs(r_t_boot[i-4])+b5*abs(r_t_boot[i-5])+b6*abs(r_t_boot[i-6])+ theta*tau_d[i]
+ z*X_boot[i-1]
r_t_boot[i]<-abs(VaR_boot[i])*resids_boot[i]
}
} else if(q==7){

b0<-stats::coef(res)[1]
b1<-stats::coef(res)[2]
b2<-stats::coef(res)[3]
b3<-stats::coef(res)[4]
b4<-stats::coef(res)[5]
b5<-stats::coef(res)[6]
b6<-stats::coef(res)[7]
b7<-stats::coef(res)[8]
z<-stats::coef(res)[9]
theta<-stats::coef(res)[10]

for(i in (q+1):N_est){
VaR_boot[i]<-b0 + b1*abs(r_t_boot[i-1]) + b2*abs(r_t_boot[i-2])+ b3*abs(r_t_boot[i-3])+
b4*abs(r_t_boot[i-4])+b5*abs(r_t_boot[i-5])+b6*abs(r_t_boot[i-6])+
b7*abs(r_t_boot[i-7])+ theta*tau_d[i] + z*X_boot[i-1]
r_t_boot[i]<-abs(VaR_boot[i])*resids_boot[i]
}
} else if(q==8){

b0<-stats::coef(res)[1]
b1<-stats::coef(res)[2]
b2<-stats::coef(res)[3]
b3<-stats::coef(res)[4]
b4<-stats::coef(res)[5]
b5<-stats::coef(res)[6]
b6<-stats::coef(res)[7]
b7<-stats::coef(res)[8]
b8<-stats::coef(res)[9]
z<-stats::coef(res)[10]
theta<-stats::coef(res)[11]

for(i in (q+1):N_est){
VaR_boot[i]<-b0 + b1*abs(r_t_boot[i-1]) + b2*abs(r_t_boot[i-2])+ b3*abs(r_t_boot[i-3])+
b4*abs(r_t_boot[i-4])+b5*abs(r_t_boot[i-5])+b6*abs(r_t_boot[i-6])+
b7*abs(r_t_boot[i-7])+
b8*abs(r_t_boot[i-8])+ theta*tau_d[i] + z*X_boot[i-1]
r_t_boot[i]<-abs(VaR_boot[i])*resids_boot[i]
} 
} else if(q==9){

b0<-stats::coef(res)[1]
b1<-stats::coef(res)[2]
b2<-stats::coef(res)[3]
b3<-stats::coef(res)[4]
b4<-stats::coef(res)[5]
b5<-stats::coef(res)[6]
b6<-stats::coef(res)[7]
b7<-stats::coef(res)[8]
b8<-stats::coef(res)[9]
b9<-stats::coef(res)[10]
z<-stats::coef(res)[11]
theta<-stats::coef(res)[12]

for(i in (q+1):N_est){
VaR_boot[i]<-b0 + b1*abs(r_t_boot[i-1]) + b2*abs(r_t_boot[i-2])+ b3*abs(r_t_boot[i-3])+
b4*abs(r_t_boot[i-4])+b5*abs(r_t_boot[i-5])+b6*abs(r_t_boot[i-6])+b7*abs(r_t_boot[i-7])+
b8*abs(r_t_boot[i-8])+b9*abs(r_t_boot[i-9])+ theta*tau_d[i] + z*X_boot[i-1]
r_t_boot[i]<-abs(VaR_boot[i])*resids_boot[i]
}

} else if(q==10){

b0<-stats::coef(res)[1]
b1<-stats::coef(res)[2]
b2<-stats::coef(res)[3]
b3<-stats::coef(res)[4]
b4<-stats::coef(res)[5]
b5<-stats::coef(res)[6]
b6<-stats::coef(res)[7]
b7<-stats::coef(res)[8]
b8<-stats::coef(res)[9]
b9<-stats::coef(res)[10]
b10<-stats::coef(res)[11]
z<-stats::coef(res)[12]
theta<-stats::coef(res)[13]


for(i in (q+1):N_est){
VaR_boot[i]<-b0 + b1*abs(r_t_boot[i-1])  + b2*abs(r_t_boot[i-2])+ b3*abs(r_t_boot[i-3])+
b4*abs(r_t_boot[i-4])+b5*abs(r_t_boot[i-5])+b6*abs(r_t_boot[i-6])+b7*abs(r_t_boot[i-7])+
b8*abs(r_t_boot[i-8])+b9*abs(r_t_boot[i-9])+b10*abs(r_t_boot[i-10])+ theta*tau_d[i]
+ z*X_boot[i-1]
r_t_boot[i]<-abs(VaR_boot[i])*resids_boot[i]
}

} 

######################################### est

r_t_boot<-as.xts(r_t_boot,time(r_t_in_s))

X_t_boot<-cbind(
rep(1,N_est),
stats::lag(abs(r_t_boot),k=1:q),
tau_d,
stats::lag(X_boot,k=1)) 

X_t_boot<-replace(X_t_boot,is.na(X_t_boot),0)

colnames(X_t_boot)[1:(q+1)]<-paste("beta","_",0:q,sep="")
colnames(X_t_boot)[ncol(X_t_boot)-1]<-c("theta")
colnames(X_t_boot)[ncol(X_t_boot)]<-c("beta_X")

## Estimation

X_t_boot<-zoo::coredata(cbind(X_t_boot[,2:(q+3)]))
Y_boot<-zoo::coredata(r_t_boot)
MAT_boot<-data.frame(Y_boot,X_t_boot)
colnames(MAT_boot)[1]<-"Y_t"

if(Exp_Short=="No"){
res_boot<-rq(formula = Y_t ~ ., tau = tau, data=MAT_boot)
} else if (Exp_Short=="Yes") {

res_boot_start<-rq(formula = Y_t ~ ., tau = tau, data=MAT_boot)
start_val_boot<-c(stats::coef(res_boot_start),0)
names(start_val_boot)[(q+4)]<-"theta_gamma"

res_boot<- maxLik(lARCH_fun_es, 
start=start_val_boot, 
ret=Y_boot,
X_t=cbind(1,X_t_boot),
tau=tau,
method="BFGS")

}
# 
}

##### third step: compute the new estimates (based on bootstrapped returns)

mat_boot_results[b,]<-stats::coef(res_boot)

utils::setTxtProgressBar(pb, b)
}  # end B

Sys.sleep(1)
close(pb)

## compute the standard errors 

sd_est<-apply(mat_boot_results,2,stats::sd)


mat_coef[,1]<-round(stats::coef(res),6)
mat_coef[,2]<-round(sd_est,6)
mat_coef[,3]<-round(stats::coef(res)/sd_est,6)
mat_coef[,4]<-round(apply(rbind(stats::coef(res)/sd_est),1,function(x) 2*(1-stats::pnorm(abs(x)))),6)



###################################################### omega_2

if(model=="lARCHMIDAS"|model=="lARCHMIDASX"){
omega_2<-w2_star} else{
omega_2<-c("Not present")
}

######## in-sample and out-of-sample estimation and evaluation

if (missing(out_of_sample)){

res_f<-list(
model=model,
coef_mat=round(mat_coef,5),
obs=N,
period=range(stats::time(r_t_in_s)),
omega_2=omega_2,
VaR=as.xts(VaR,stats::time(r_t_in_s)),
ES=ES,
hit_in_s=sum(ifelse(r_t_in_s<=VaR,1,0))/N,
loss_in_s=QL(r_t_in_s,VaR,tau)
)

} else {

res_f<-list(
model=model,
coef_mat=round(mat_coef,5),
obs=N,
period=range(stats::time(r_t_in_s)),
omega_2=omega_2,
VaR=as.xts(VaR,stats::time(r_t_in_s)),
ES=ES,
hit_in_s=sum(ifelse(r_t_in_s<=VaR,1,0))/N,
loss_in_s=QL(r_t_in_s,VaR,tau),
VaR_oos=as.xts(VaR_oos,stats::time(r_t_oos)),
ES_oos=ES_oos
)
}

class(res_f)<-c("rqmidas")

return(res_f)
print.rqmidas(res_f)


}

##

###############################################################################################
# Computes the matrix power of a symmetric matrix using its spectral decomposition
# Input  : x -> square symmetric matrix
#          p -> exponent
# Output : matrix raised to the power p
###############################################################################################


"%^%" <- function(x, p) 
   with(eigen(x), vectors %*% (values^p* t(vectors)))


