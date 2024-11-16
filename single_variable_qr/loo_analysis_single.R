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

############## single-predictor
X_single = X[,'sqrtSqft_living15', drop=FALSE]
n = 200
p = 1
subsample = sample(n_total, 2*n)
train_idx = subsample[1:n]

X_train = X_single[train_idx, , drop = FALSE]
X_tilde <- cbind(1,X_train)
Y_train = Y[train_idx]

a = -rep(1, n)
B = -X_train
c = Y_train
lam = n*0.01

tau_values <- c(0.1, 0.5, 0.95)
beta <- CVXR::Variable(p)
beta_0 = CVXR::Variable()
xi = CVXR::Variable(n)
eta = CVXR::Variable(n)
lambda = Parameter(pos = TRUE)

plot_list <- list()
x_seq <- seq(min(X_train), max(X_train), length.out = 100)

for (i in seq_along(tau_values)) {
  tau <- tau_values[i]
  alpha_0 = tau
  alpha_1 = 1 - tau
  
  obj = CVXR::sum_entries(alpha_0 * xi + alpha_1 * eta) + lambda / 2 * CVXR::power(CVXR::norm2(beta), 2)
  constraints = list(xi >= 0, eta >= 0, beta_0 * a + B %*% beta + c + eta >= 0, - beta_0 * a - B %*% beta - c + xi >= 0)
  prob <- CVXR::Problem(Minimize(obj), constraints)
  
  value(lambda) = n*0.01
  result <- CVXR::solve(prob, warm_start = TRUE)
  theta_full = result$getDualValue(constraints[[3]]) - result$getDualValue(constraints[[4]])
  r_train = Y_train - (X_train %*% result$getValue(beta) + result$getValue(beta_0))
  fitted_value_full <- X_train %*% result$getValue(beta) + result$getValue(beta_0)
  
  eps = 1e-3
  Elbow_idx <- which(abs(r_train)<=eps)
  Left_Elbow_idx <- which(r_train<=-eps)
  Right_Elbow_idx <- which(r_train>=eps)
  
  # Calculate LOO influence for all indices
  loo_influences <- sapply(1:n, calculate_loo_influence, X_train, y_train, tau, lam, fitted_value_full)
  
  top_5_index <- order(loo_influences, decreasing = TRUE)[1:5]
  solution_full <- list(beta = result$getValue(beta), beta_0 = result$getValue(beta_0))
  betas <- list(full = solution_full)
  
  for (j in 1:length(top_5_index)) {
    loo_solution <- calculate_loo_solution(top_5_index[j],X_train,y_train,tau,lam)
    betas[[paste0("loo_", j)]] <- loo_solution
  }

plot_data <- data.frame(
  x = rep(x_seq, times = length(betas)),
  y = unlist(lapply(betas, function(beta) beta$beta_0 + beta$beta * x_seq)),
  solution = rep(names(betas), each = length(x_seq))
)

colors <- c("full" = "black", "loo_1" = "#1F77B4",
            "loo_2" = "#FF7F0E", "loo_3" = "#2CA02C", 
            "loo_4" = "#D62728", "loo_5" = "#9467BD")

if(i==1){
P <- ggplot(plot_data, aes(x = x, y = y, color = solution, linetype = solution)) +
  geom_line(size = 1, alpha = .7) +
  scale_color_manual(values = colors) +
  scale_linetype_manual(values = c("full" = "solid", "loo_1" = "dashed", "loo_2" = "dashed", "loo_3" = "dashed", "loo_4" = "dashed", "loo_5" = "dashed")) +
  labs(x = "Sqft_living15", y = "Predicted house price", color = "Solution", 
       title = bquote(tau == .(sprintf("%.2f", tau_values[i])))) +
  guides(linetype = FALSE) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14), 
    legend.position = "none",
    axis.title.x = element_text(size = 14), 
    axis.title.y = element_text(size = 14)) + 
  ylim(400, 1500) 
} else{
  P <- ggplot(plot_data, aes(x = x, y = y, color = solution, linetype = solution)) +
    geom_line(size = 1, alpha = .7) +
    scale_color_manual(values = colors) +
    scale_linetype_manual(values = c("full" = "solid", "loo_1" = "dashed", "loo_2" = "dashed", "loo_3" = "dashed", "loo_4" = "dashed", "loo_5" = "dashed")) +
    labs(x = "Sqft_living15", color = "Solution",
         title = bquote(tau == .(sprintf("%.2f", tau_values[i])))) +
    guides(linetype = FALSE) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14), 
      legend.position = "none",
      axis.title.x = element_text(size = 14), 
      axis.title.y = element_blank()) + 
    ylim(400, 1500) 
}

plot_list[[i]] <- P
} 

grid.arrange(grobs = plot_list, nrow = 1, ncol = 3)
# 12 * 4 inch