# Replication Package

## Combining Value-at-Risk and Expected Shortfall Forecasts via the Model Confidence Set

Amendola A., Candila V., Naimoli A., and G. Storti (2026),  
*International Journal of Forecasting.*

---

## 📅 Assembly Date and Contact

**Package assembled:** April 2026  
**Contact:** Vincenzo Candila — vcandila@unisa.it  
*(Please reach out for any questions about the code or data)*

---

## 🗂️ Repository Structure

```
combining_var_es_mcs_combinations/
├── README.md
├── README.html
├── replicate_tables_and_figures.R            # Replicates all tables and figures in the paper
├── generate_proposed_combinations.R          # Generates the six proposed combined predictors
├── functions/
│   ├── main_functions.R                      # Core functions used by both scripts
│   ├── AL_lambda.R                           # Asymmetric Laplace loss (lambda variant)
│   ├── AL_msc.R                              # Asymmetric Laplace loss (MSC variant)
│   ├── beta_esti_ms.R                        # Beta coefficient estimation (minimum score)
│   ├── lambda_esti.R                         # Lambda weight estimation
│   ├── rel_comb.R                            # Relative combination utility
│   └── transbeta.R                           # Beta transformation for quantile/ES weights
└── data/
    ├── SP500_data.RData                      # S&P 500 daily returns and realized volatility measures
    ├── Shanghai_comp_data.RData              # Shanghai Composite returns and realized volatility measures
    ├── Global_Policy_Uncertainty_Data.csv    # Global Economic Policy Uncertainty index (EPU)
    ├── res_sp500_tau_0.025.RData             # Pre-computed model results for S&P 500 (τ = 2.5%)
    ├── res_shanghai_comp_tau_0.025.RData     # Pre-computed model results for Shanghai Comp. (τ = 2.5%)
    ├── all_indices_tau_0.025.RData           # Result tables for all nine indices (τ = 2.5%)
    ├── all_indices_tau_0.01.RData            # Result tables for all nine indices (τ = 1%)
    └── shanghai_comp_tau_0.025_comb_example.RData  # Example dataset to be used with generate_proposed_combinations.R
```

---

## 💻 Computing Environment

- **Language:** R version 4.4.3
- **Operating system:** Windows 10 (code is platform-independent)

### Required R packages

| Package     | Version (used) | Purpose                                      |
|-------------|----------------|----------------------------------------------|
| `xts`       | ≥ 0.13         | Time series management                       |
| `zoo`       | ≥ 1.8          | Time series infrastructure                   |
| `fBasics`   | ≥ 4032         | Summary statistics                           |
| `rugarch`   | ≥ 1.5          | GARCH model estimation                       |
| `GAS`       | ≥ 0.3          | Generalized Autoregressive Score models      |
| `np`        | ≥ 0.60         | Non-parametric kernel methods                |
| `esback`    | ≥ 0.1          | Expected Shortfall backtesting               |
| `Rsolnp`    | ≥ 1.16         | Non-linear optimization (combinations only)  |
| `DT`        | ≥ 0.30         | Interactive HTML tables                      |
| `htmltools` | ≥ 0.5          | HTML rendering utilities                     |

Install all packages at once with:

```r
install.packages(c("xts", "zoo", "fBasics", "rugarch", "GAS",
                   "np", "esback", "Rsolnp", "DT", "htmltools"))
```

No environment manager (e.g., `renv`) was used. To record a snapshot of your installed versions, run `sessionInfo()` in R after installing the packages.

---

## 🗃️ Data

### SP500 and Shanghai Composite data

- **Files:** `data/SP500_data.RData`, `data/Shanghai_comp_data.RData`
- **Source:** Oxford-Man Institute's Realized Library 
- **Period used:** January 2013 – June 2022
- **Variables included in each file:**
  - `r_t`: daily close-to-close log returns
  - `rb_ss`: realized bipower variation with subsampling
  - `rk`: realized kernel volatility
  - `rvol5`: realized volatility based on 5-minute intervals
- **License:** Publicly available for academic use; please consult the Oxford-Man Institute terms.

### Global Economic Policy Uncertainty (EPU)

- **File:** `data/Global_Policy_Uncertainty_Data.csv`
- **Source:** https://www.policyuncertainty.com/global_monthly.html  
  Column used: `GEPU_current`. Downloaded on **May 26, 2023**.
- **Note:** This series may be subject to revisions over time; the version included in the package reflects the data as of the download date.

### Estimation results

| File | Description |
|------|-------------|
| `data/res_sp500_tau_0.025.RData` | VaR/ES results for S&P 500 (τ = 2.5%) |
| `data/res_shanghai_comp_tau_0.025.RData` | VaR/ES results for Shanghai Composite (τ = 2.5%) |
| `data/all_indices_tau_0.025.RData` | Result tables for all nine indices (τ = 2.5%) |
| `data/all_indices_tau_0.01.RData` | Result tables for all nine indices (τ = 1%) |
| `data/shanghai_comp_tau_0.025_comb_example.RData` | Forecast arrays for the combination example |


---

## 📊 Code → Output Mapping

### `replicate_tables_and_figures.R`

Run this script to replicate all tables and figures in the paper. Results are rendered as HTML tables and inline plots.

| Output     | Script section | Data required  |
|:-----------|:---------------|:---------------|
| Table 2 — Summary statistics | `#### Replication of Table 2` | `SP500_data.RData`, `Shanghai_comp_data.RData`, `Global_Policy_Uncertainty_Data.csv` |
| Figure 2(a) — Inclusion plot | `#### Replication of Figure 2(a)` | `res_shanghai_comp_tau_0.025.RData` |
| Figure 2(b) — Inclusion plot | `#### Replication of Figure 2(b)` | `res_shanghai_comp_tau_0.025.RData` |
| Figure 2(c) — Inclusion plot | `#### Replication of Figure 2(c)` | `res_shanghai_comp_tau_0.025.RData` |
| Table 5 — S&P 500 out-of-sample evaluation | `#### Replication of Table 5` | `res_sp500_tau_0.025.RData` |
| Table 6 — Shanghai Comp. out-of-sample evaluation | `#### Replication of Table 6` | `res_shanghai_comp_tau_0.025.RData` |
| Table 7 — Results across all indices and risk levels | `#### Replication of Table 7` | `all_indices_tau_0.025.RData`, `all_indices_tau_0.01.RData` |

**Note on Table 5 and Table 6:** To reproduce exactly the p-values of the Bayer & Dimitriadis (2020) tests reported in the paper, the code should be run under R version 4.4.3.
When using different R versions or computing environments, small numerical differences may arise; these affect only the p-values of these tests and do not alter the main conclusions.


### `generate_proposed_combinations.R`

Run this script to generate the six proposed forecast combinations. For ease of implementation and faster execution, it relies on a **simplified (reduced) version** of the original exercise, using the Shanghai Composite example dataset.

| Output | Description |
|--------|-------------|
| Saved `.RData` file (name set by `filename` variable) | Contains the six combined VaR and ES forecast series: `MCS_Comb`, `WL_MCS_Comb`, `MCS_RSC_Comb`, `WL_MCS_RSC_Comb`, `MCS_MSC_Comb`, `WL_MCS_MSC_Comb` |

Progress is printed to the console at each forecasting step via `message(tt)`.

---

## ⏱️ Hardware and Expected Runtime

| Script         | Hardware used | Approximate runtime |
|:---------------|:--------------|:--------------------|
| `replicate_tables_and_figures.R` | Standard laptop, 16 GB RAM | ~1–2 minutes |
| `generate_proposed_combinations.R` | Standard laptop, 16 GB RAM | ~3–4 minutes|


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
