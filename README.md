---
output:
  html_document: default
  pdf_document: default
---
# Replication Package

## Combining Value-at-Risk and Expected Shortfall Forecasts via the Model Confidence Set

Amendola A., Candila V., Naimoli A., and G. Storti (2026),  
*International Journal of Forecasting.*

---

## 📅 Assembly Date and Contact

**Package assembled:** 3 June 2026  
**Contact:** Vincenzo Candila — vcandila@unisa.it  
*(Please reach out for any questions about the code or data)*

---

## 🗂️ Repository Structure

```
combining_var_es_mcs_combinations/
├── README.md
├── README.html
├── 1_replicate_table_2.R                         # Replicates Table 2
├── 2_from_raw_files_to_intermediary_files.R      # Computes the VaR and ES forecasts for the 32 models used in the paper
├── 3_from_intermediary_files_to_final_files.R    # Computes the proposed forecast combinations and benchmarks
├── 4_replicate_figure_2_and_tables_5_to_7.R      # Replicates Figure 2 and Table 5, 6, and 7
├── functions/
    ├── main_functions.R                          # Core functions used by all scripts
    ├── AL_lambda.R                               # Asymmetric Laplace loss (lambda variant)
    ├── AL_msc.R                                  # Asymmetric Laplace loss (Minimum score variant)
    ├── beta_esti_ms.R                            # Beta coefficient estimation (Minimum score)
    ├── lambda_esti.R                             # Lambda weight estimation
    ├── rel_comb.R                                # Relative combination utility
    ├── msc_main.R                      # Implementation of the Minimum Score Combining (MSC) method of Taylor (2020)
    ├── rsc_main.R                      # Implementation of the Relative Score Combining (MSC) method of Taylor (2020)
    └── transbeta.R                               # Beta transformation for quantile/ES weights
└── data/
    ├── raw/
        ├── sp500_data.RData                      # S&P 500 daily returns and realized volatility measures
        ├── shanghai_comp_data.RData              # Shanghai Composite returns and realized volatility measures
        ├── bovespa_data.RData                    # BOVESPA daily returns and realized volatility measures
        ├── bsesn_data.RData                      # BSESN returns and realized volatility measures
        ├── eurostoxx50_data.RData                # EUROSTOXX50 daily returns and realized volatility measures
        ├── hsi_data.RData                        # Hang Seng returns and realized volatility measures
        ├── ixic_data.RData                       # NASDAQ returns and realized volatility measures
        ├── mxx_data.RData                        # MXX returns and realized volatility measures
        ├── nikkei_data.RData                     # NIKKEI returns and realized volatility measures
        └── Global_Policy_Uncertainty_Data.csv    # Global Economic Policy Uncertainty index (EPU)
    ├── intermediary/
        ├── sp500_tau_0.025_files.RData           # S&P 500 intermediary files needed to compute the proposed combinations
        └── shanghai_comp_tau_0.025_files.RData   # Shanghai Comp. intermediary files needed to compute the proposed comb.
    └── results/
        ├── sp500_final_results_precomputed.RData             # Results for S&P 500 (τ = 2.5%)
        ├── shanghai_comp_final_results_precomputed.RData     # Results for Shanghai Comp. (τ = 2.5%)
        ├── all_indices_tau_0.025.RData                       # Result tables for all nine indices (τ = 2.5%)
        └── all_indices_tau_0.01.RData                        # Result tables for all nine indices (τ = 1%)
```

---

## 💻 Computing Environment

- **Language:** R version 4.4.3
- **Operating system:** Windows 10 (code is platform-independent)

### Required R packages

| Package     | Version (used)     | Purpose                                            |
|-------------|--------------------|----------------------------------------------------|
| `xts`       | 0.14.1             | Time series management                             |
| `zoo`       | 1.8.13             | Time series infrastructure                         |
| `fBasics`   | 4041.97            | Summary statistics                                 |
| `rugarch`   | 1.5.3              | Model Confidence Set (MCS) procedure via `mcsTest` |
| `GAS`       | 0.3.4.1 (Currently, the GAS package has been removed from the CRAN. But previous versions are available at https://cran.r-project.org/src/contrib/Archive/GAS/)           | Functions for VaR backtesting (`BacktestVaR`) and  Fissler-Ziegel loss evaluation (`FZLoss`)   |   
| `np`        | 0.60.18            | Automatic block length selection via `b.star` for the bootstrap procedure  within the MCS     |
| `esback`    | 0.3.1              | Expected Shortfall backtesting                     |
| `Rsolnp`    | 1.16               | Non-linear optimization (combinations only)        |
| `DT`        | 0.33               | Interactive HTML tables                            |
| `htmltools` | 0.5.8.1             | HTML rendering utilities                           |
| `fGarch`    | 4033.92             |  Estimation of standardized Student-t distribution parameters via `stdFit`                           |
| `rumidas` | 0.1.3            | MIDAS matrix construction via `mv_into_mat`                          |
| `quantreg` | 6.1             | Quantile regression estimation via `rq`                        |



Install all packages at once with:

```r
install.packages(c("xts", "zoo", "fBasics", "rugarch",
                   "np", "esback", "Rsolnp", "DT", "htmltools", "fGarch", "rumidas", "quantreg"))
# GAS has been removed from CRAN and must be installed manually from the CRAN archive:
install.packages("https://cran.r-project.org/src/contrib/Archive/GAS/GAS_0.3.4.1.tar.gz",
                 repos = NULL,
                 type = "source")
```

No environment manager (e.g., `renv`) was used. To record a snapshot of your installed versions, run `sessionInfo()` in R after installing the packages.

---

## 🔄 Workflow

The replication workflow is organized into three main stages:

1. **Raw data**
   - Raw financial and macroeconomic data are stored in `data/raw/`.
   - These include daily returns and realized volatility measures for all indices considered in the paper, as well as the Global Economic Policy Uncertainty (EPU) index.
   - Table 2 is replicated using the raw data and the script `1_replicate_table_2.R`.

2. **Intermediary data**
   - The raw data were originally used to generate the intermediary `.RData` objects stored in `data/intermediary/`.
   - These intermediary files contain the VaR and ES forecasts produced by the 32 competing models.
   - This allows users to fully reproduce the forecast combinations, tables, figures, and empirical results reported in the paper without regenerating the intermediary objects from the raw data.
   - Alternatively, the intermediary files can be regenerated from the raw data by running the script `2_from_raw_files_to_intermediary_files.R`.

3. **Results**
   - Pre-computed results are stored in `data/results/`.
   - The script `3_from_intermediary_files_to_final_files.R` (if executed) uses the intermediary datasets to compute the six proposed forecast combinations and the four benchmark methods. The resulting objects are stored in `data/results/` and allow the replication of Figure 2 and Tables 5 and 6 through the script `3_replicate_figure_2_and_tables_5_to_7.R`.
   - Alternatively, the script `4_replicate_figure_2_and_tables_5_to_7.R` reproduces Figure 2 and Tables 5–7 directly using the pre-computed files included in `data/results/`.

## 🗃️ Raw Data

### Index data

- **Files:** `data/raw/sp500_data.RData`, `data/raw/shanghai_comp_data.RData`, `data/raw/bovespa_data.RData`, `data/raw/bsesn_comp_data.RData`, `data/raw/eurostoxx50_data.RData`, `data/raw/hsi_comp_data.RData`, `data/raw/ixic_comp_data.RData`, `data/raw/mxx_data.RData`, `data/raw/nikkei_comp_data.RData`
- **Source:** Oxford-Man Institute's Realized Library 
- **Period used:** January 2013 – June 2022
- **Variables included in each file:**
  - `r_t`: daily close-to-close log returns
  - `rb_ss`: realized bipower variation with subsampling
  - `rk`: realized kernel volatility
  - `rvol5`: realized volatility based on 5-minute intervals
- **License:** Publicly available for academic use; please consult the Oxford-Man Institute terms.

### Global Economic Policy Uncertainty (EPU)

- **File:** `data/raw/Global_Policy_Uncertainty_Data.csv`
- **Source:** https://www.policyuncertainty.com/global_monthly.html  
  Column used: `GEPU_current`. Downloaded on **May 26, 2023**.
- **Note:** This series may be subject to revisions over time; the version included in the package reflects the data as of the download date.
- **License:** Publicly available for academic use; please consult the original website for terms of use.

## 📦 Intermediary Data

- **Files:** `data/intermediary/sp500_tau_0.025_files.RData`, `data/intermediary/shanghai_comp_tau_0.025_files.RData`
- **Variables included in each file:**
  - `Tin`					 : length of the training period
  - `nstep`					 : number of forecasting steps (model re-estimations) and out-of-sample period length
  - `N_model`					 : number of competing models (nmod)
  - `tau`						 : coverage level used	
  - `VaR_training_data_mod`	 : array of VaR forecasts over the training period  - Tin x nmod x nstep
  - `ES_training_data_mod`	  : array of ES forecasts over the training period   - Tin x nmod x nstep
  - `r_t_in_s_matrix`	       : matrix of returns over the training period       - Tin x nstep
  - `VaR_oos`	               : matrix of out-of-sample VaR forecasts            - nstep x nmod
  - `ES_oos`	                : matrix of out-of-sample ES forecasts             - nstep x nmod
  - `r_t_oos_full`	          : matrix of out-of-sample returns (xts/zoo format) - nstep x 1
  - `r_t_oos`	               : matrix of out-of-sample returns                  - nstep x 1
  - `list_of_models`	        : character vector of labels for the nmod competing models

## 📊 Results Data

- **Files:** `data/results/sp500_final_results_precomputed.RData`, `data/results/shanghai_comp_final_results_precomputed.RData`, `data/results/all_indices_tau_0.025.RData`, `data/results/all_indices_tau_0.01.RData.RData`, 
- **Variables included in** `data/results/sp500_final_results_precomputed.RData` and `data/results/shanghai_comp_final_results_precomputed.RData`:
- `r_t_oos_full_plot`: daily log-returns used for the training-period plots.
- `r_t_oos_full`: daily log-returns for the out-of-sample period.
- `nstep`: total number of model re-estimations, corresponding to the length of the out-of-sample period.
- `lab_full`: model labels, including the 32 individual models, the four benchmarks (`EW-Comb`, `Median-Comb`, `RS-Comb`, `MS-Comb`), and the six proposed forecast combinations (`MCS-Comb`, `WL-MCS-Comb`, `MCS-RS-Comb`, `WL-MCS-RS-Comb`, `MCS-MS-Comb`, `WL-MCS-MS-Comb`).
- `tau`: coverage level.
- `Backtesting_pvalues`: p-values of the six backtests reported in Table 4 of the paper, by model and forecasting step.
- `db_loss_oos`: matrix of Fissler-Ziegel losses for the out-of-sample period for the 42 models and forecast combinations.
- `VaR_oos_ev` and `ES_oos_ev`: matrices of VaR and ES forecasts for the out-of-sample period for the 42 models and forecast combinations.
- `MCS_est_included`: list of models included at each forecasting step in the Set of Superior Models (SSM) from the Model Confidence Set (MCS), using the unweighted Fissler-Ziegel loss.
- `MCS_est_included_w`: list of models included at each forecasting step in the Set of Superior Models (SSM) from the Model Confidence Set (MCS), using the weighted Fissler-Ziegel loss.
- **Variables included in** `data/results/all_indices_tau_0.025.RData` and `data/results/all_indices_tau_0.01.RData`:
- `list_of_tabs`: list of tables (one for each index) containing models in rows and backtesting p-values (`UC`, `CC`, `DQ`, `BD-1`, `BD-2`, `BD-3`) together with the Fissler-Ziegel loss in columns. LaTeX markers (e.g., `\cellcolor{...}`) indicate whether models simultaneously pass all backtests at the 5% significance level and belong to the Set of Superior Models (SSM) from the Model Confidence Set (MCS) at the 25% significance level.
- `lab_full`: model labels, including the 32 individual models, the four benchmarks (`EW-Comb`, `Median-Comb`, `RS-Comb`, `MS-Comb`), and the six proposed forecast combinations (`MCS-Comb`, `WL-MCS-Comb`, `MCS-RS-Comb`, `WL-MCS-RS-Comb`, `MCS-MS-Comb`, `WL-MCS-MS-Comb`).

---

## 📊 Code → Output Mapping

| Output     | Script  | Data required  |
|:-----------|:---------------|:---------------|
| Table 2 — Summary statistics | `code/1_replicate_table_2.R` | Files in `data/raw/`: `SP500_data.RData`, `Shanghai_comp_data.RData`, `Global_Policy_Uncertainty_Data.csv` |
| Figure 2(a) — Inclusion plot | `code/4_replicate_figure_2_and_tables_5_to_7.R` | Files in `data/results/`:Either `shanghai_comp_final_results_precomputed.RData` or `shanghai_comp_final_results_generated.RData` generated by `code/3_from_intermediary_files_to_final_files.R` |
| Figure 2(b) — Inclusion plot | `code/4_replicate_figure_2_and_tables_5_to_7.R` | Files in `data/results/`:Either `shanghai_comp_final_results_precomputed.RData` or `shanghai_comp_final_results_generated.RData` generated by `code/3_from_intermediary_files_to_final_files.R` |
| Figure 2(c) — Inclusion plot | `code/4_replicate_figure_2_and_tables_5_to_7.R` | Files in `data/results/`:Either `shanghai_comp_final_results_precomputed.RData` or `shanghai_comp_final_results_generated.RData` generated by `code/3_from_intermediary_files_to_final_files.R` |
| Table 5 — S&P 500 out-of-sample evaluation | `code/4_replicate_figure_2_and_tables_5_to_7.R` | Files in `data/results/`: Either `sp500_final_results_precomputed.RData` or `sp500_final_results_generated.RData` generated by `code/3_from_intermediary_files_to_final_files.R` |
| Table 6 — Shanghai Comp. out-of-sample evaluation | `code/4_replicate_figure_2_and_tables_5_to_7.R` | Files in `data/results/`: Either `shanghai_comp_final_results_precomputed.RData` or `shanghai_comp_final_results_generated.RData` generated by `code/3_from_intermediary_files_to_final_files.R` |
| Table 7 — Results across all indices and risk levels | `code/4_replicate_figure_2_and_tables_5_to_7.R` | Files in `data/results/`: `all_indices_tau_0.01.RData` and `all_indices_tau_0.025.RData` |

**Note on Table 5 and Table 6:** To reproduce exactly the p-values of the Bayer & Dimitriadis (2020) tests reported in the paper, the code should be run under R version 4.4.3.
When using different R versions or computing environments, small numerical differences may arise; these affect only the p-values of these tests and do not alter the main conclusions.


---

## ⏱️ Hardware and Expected Runtime

| Script         | Hardware used               | Approximate runtime |
|:---------------|:----------------------------|:--------------------|
| `1_replicate_table_2.R` | Standard laptop (Intel i5-1135G7 @ 2.4 GHz, 16 GB RAM) | ~few seconds |
| `2_from_raw_files_to_intermediary_files.R` | Dedicated workstation (Intel(R) Core(TM) i9-14900KF (3.20 GHz), 64 GB RAM) | ~40 hours per index |
| `3_from_intermediary_files_to_final_files.R` | Dedicated workstation (Intel(R) Core(TM) i9-14900KF (3.20 GHz), 64 GB RAM) | ~16 hours per index|
| `4_replicate_figure_2_and_tables_5_to_7.R` | Standard laptop (Intel i5-1135G7 @ 2.4 GHz, 16 GB RAM) | ~3-4 minutes|


---

## ⚙️ Special Setup Requirements

No special infrastructure is required. The code runs on a standard desktop or laptop computer. No GPU, cluster, or parallel computing setup is needed.

Make sure the **working directory** in R is set to the root of the replication folder (where this README is located) before running either script, so that all relative paths (`data/`, `functions/`) resolve correctly:

```r
setwd("path/to/combining_var_es_mcs_replication/")
```

---

## 📜 Licence

The code in this repository is shared for academic reproducibility purposes. If you use it, please cite the paper:

> Amendola A., Candila V., Naimoli A., and G. Storti (2026), "Combining Value-at-Risk and Expected Shortfall forecasts via the Model Confidence Set", *International Journal of Forecasting*.
