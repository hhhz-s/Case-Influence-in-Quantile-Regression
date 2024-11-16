rm(list=ls())

## all related packages are already listed in helper.R
## pkgs <- c("CVXR","ggplot2","gridExtra","quantreg","NonsmoothPath","Caseinflu")

# specify the working dir
setwd("path/to/github/repo") 

# source helper file
source("single_variable_qr/helper.R")

# Load data
load("data/Input.RData")

set.seed(1)
X = as.matrix(X)
n_total = nrow(X)
n = 200
subsample = sample(n_total, 2*n)
train_idx = subsample[1:n]
test_idx = subsample[(n+1):(2*n)]

X <- X[, !colnames(X) %in% c("waterfront1", "condition_average","zone_1","zone_2")]

bedrooms_3 <- ifelse(X[, "bedrooms_12"] == 0, 1, 0)
X <- cbind(X, bedrooms_3)
X <- X[, !colnames(X) %in% c("bedrooms_35", "bedrooms_12")]

view_recode <- ifelse(X[, "view"] == 0, 0, 1)
X <- cbind(X, view_recode)
X <- X[, !colnames(X) %in% "view"]

X_train = X[train_idx,]
p = ncol(X_train)
X_tilde <- cbind(1,X_train)
Y_train = Y[train_idx]
sigma <- sigma_func(X_train, Y_train)

a = -rep(1, n)
B = -X_train
c = Y_train

X_test = X[test_idx,]
Y_test = Y[test_idx]


# lambda value from the cross validation result in loo_influ_residual_leverage_multi.R
lam_list <- c(14.3845, 0.01623777, 0.02636651)
tau_values <- c(0.1, 0.5, 0.95)
beta <- CVXR::Variable(p)
beta_0 = CVXR::Variable()
xi = CVXR::Variable(n)
eta = CVXR::Variable(n)
lambda = Parameter(pos = TRUE)

index_list <- c(0,0,0)

for (i in seq_along(tau_values)) {
  tau <- tau_values[i]
  alpha_0 = tau
  alpha_1 = 1 - tau
  
  obj = CVXR::sum_entries(alpha_0 * xi + alpha_1 * eta) + lambda / 2 * CVXR::power(CVXR::norm2(beta), 2)
  constraints = list(xi >= 0, eta >= 0, beta_0 * a + B %*% beta + c + eta >= 0, - beta_0 * a - B %*% beta - c + xi >= 0)
  prob <- CVXR::Problem(Minimize(obj), constraints)
  
  value(lambda) = lam_list[i]
  lam = lam_list[i]
  result <- CVXR::solve(prob, warm_start = TRUE)
  theta_full = result$getDualValue(constraints[[3]]) - result$getDualValue(constraints[[4]])
  r_train = Y_train - (X_train %*% result$getValue(beta) + result$getValue(beta_0))
  fitted_value_full <- X_train %*% result$getValue(beta) + result$getValue(beta_0)
  
  loo_influences <- sapply(1:n, calculate_loo_influence, X_train, y_train, tau=tau,lam=lam,
                           fitted_value_full = fitted_value_full)
  
  top_index <- which.max(loo_influences)
  index_list[i] <- top_index
  print(tau)
  print(top_index)
} 

quantile_regression_with_deletion <- function(X_train, Y_train, tau = tau, index) {
  full_model <- rq(Y_train ~ ., data = as.data.frame(X_train), tau = tau)
  full_summary <- summary(full_model, se = 'boot')
  
  full_coef <- round(coef(full_model),2)
  full_p_values <- round(full_summary$coefficients[,4],2)
  
  full_model_df <- data.frame(Variable = names(full_coef),
                              Coefficients = full_coef,
                              P_values = full_p_values)
  print('Full Data Solution:')
  print(full_model_df)
 
  X_train_subset <- X_train[-index, ]
  Y_train_subset <- Y_train[-index]
  loo_model <- rq(Y_train_subset ~ ., data = as.data.frame(X_train_subset), tau = tau)
  loo_model_summary <- summary(loo_model, se = 'boot')
    
  loo_coef <- round(coef(loo_model),2)
  loo_p_values <- round(loo_model_summary$coefficients[, 4],2)

  loo_model_df <- data.frame(Variable=names(loo_coef),
                             coefficients = loo_coef,
                             P_values = loo_p_values)
  print('LOO Solution:')
  print(loo_model_df)
}


for (i in 1:3) {
  print(tau_values[i])
  quantile_regression_with_deletion(X_train, Y_train, tau = tau_values[i],index = index_list[i])
  print('-------------------------')
}