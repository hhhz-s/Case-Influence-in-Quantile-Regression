# Case-Influence-in-Quantile-Regression

## Title: 
Cross validation for penalized quantile regression with a case-weight adjusted solution path

## Introduction:
This archive contains the R codes to reproduce figures and tables for the King County Housing Data Analysis (Section 5) in the [paper(v2)](https://arxiv.org/abs/1902.07770)

## List of files: 
### Data: 
- input.RData: contains the [King County house sales data](https://www.kaggle.com/datasets/harlfoxem/housesalesprediction) preprocessed by data_preprocess.R
- data_preprocess.R: code to pre-process the data by Bo, Luan
- EDA_sample.R: EDA for a subset of the King County house sales data in Fig 14
  
### Single variable QR: 
- helper.R: customized helper functions for quantile regression
- loo_influ_y_x_single.R: code for Fig 10
- loo_analysis_single.R: code for Fig 11
- case_influence_graph_single.R: code for Fig 15

### Multi-variable QR:
- helper.R: customized helper functions for quantile regression
- loo_influ_residual_leverage_multi.R: code for Fig 12
- case_influence_graph_multi.R: code for Fig 13
- loo_analysis_multi: code for Tables 6 to 8

Specify the working directory and install the required packages in R to reproduce the figures. 
