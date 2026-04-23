
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

tab_7_dt_f <- function(x) {
 
  # Assign column names
  colnames(x) <- c(
    "BT_pass_0.025", "MCS_in_0.025", "BT_and_MCS_0.025",
    "BT_pass_0.01",  "MCS_in_0.01",  "BT_and_MCS_0.01"
  )
  
  # Compute best and second-best values by column
  max_vals <- apply(x, 2, max, na.rm = TRUE)
  
  # secondo valore distinto più alto per colonna
  second_vals <- apply(x, 2, function(x) {
    ux <- sort(unique(x), decreasing = TRUE)
    if (length(ux) >= 2) ux[2] else NA
  })
  
  tab_dt <- datatable(
    x,
    rownames = TRUE,
    class = "compact",
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
            
            // best = verde scuro
            if (!isNaN(val) && val === max_vals[i]) {
              cell.css('background-color', '#39e639');
            }
            
            // second-best = verde più chiaro
            else if (!isNaN(val) && val === second_vals[i]) {
              cell.css('background-color', '#d9f7d9');
            }
            
            // bold per 9
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



