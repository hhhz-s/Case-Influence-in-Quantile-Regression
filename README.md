# Case-Influence-in-Quantile-Regression

## Title: 
Cross validation for penalized quantile regression with a case-weight adjusted solution path

## Introduction:
This archive contains the R codes to reproduce figures and tables for the King County Housing Data Analysis (Section 5) in the [paper](https://arxiv.org/abs/1902.07770)

## List of files: 
Data: 
- input.RData: contains the King County house sales data, preprocessed by data_preprocess.R
Code: 
- data_preprocess.R: code to pre-process the data (mention Bo's work)
- helper.R: customized helper functions for quantile regression
- loo_influ_y_x_single.R: code to reproduce Fig 10
- loo_analysis_single.R: code to reproduce Fig 11
- loo_influ_residual_leverage_multi.R: code to reproduce Fig 12
- case_influence_graph_multi.R: code to reproduce Fig 13
- EDA_sample.R: EDA for a subset of the King County house sales data, used to reproduce Fig 14
- case_influence_graph_single.R: code to reproduce Fig 15
- loo_analysis_multi: code to reproduce Tables 6 to 8
