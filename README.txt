This document illustrates how to replicate the tables and figures, as well as the proposed combined predictors, presented in:

Amendola A., Candila V., Naimoli A., and G. Storti (2026), Combining Value-at-Risk and Expected Shortfall forecasts via the Model Confidence Set, International Journal of Forecasting.

In particular, the shared folder includes:

- code/1_replicate_table_2.R 				: R script to reproduce Table and 2;
- code/2_from_raw_files_to_intermediary_files.R 	: R script to compute all the VaR and ES for the 32 models, starting from the raw data;
- code/3_from_intermediary_files_to_final_files.R 	: R script to compute the six proposed forecast combinations and the four benchmark methods from the intermediary datasets;
- code/4_replicate_figure_2_and_tables_5_to_7.R 	: R script to reproduce Figure 2 and Tables 5–7 using the pre-computed result files included in data/results/ or the files obtained from '2_from_intermediary_files_to_final_files.R';
- functions/ 						: directory containing all the functions required by the three scripts above;
- data/raw/ 						: directory containing the raw datasets;
- data/intermediary/ 					: directory containing the intermediary datasets generated from the raw data and used to compute the proposed forecast combinations and benchmark methods;
- data/results/ 					: directory containing the final pre-computed result files used to reproduce Figure 2 and Tables 5–7.

For any comments or requests, please contact Vincenzo Candila at vcandila at unisa.it.
