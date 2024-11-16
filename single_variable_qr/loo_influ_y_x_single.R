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
X_single = X[,'sqrtSqft_living15', drop=FALSE]
n_total = nrow(X)
p = ncol(X_single)
n = 200 

subsample = sample(n_total, 2*n)
train_idx = subsample[1:n]
test_idx = subsample[(n+1):(2*n)]

X_train = X_single[train_idx, , drop = FALSE]
X_tilde <- cbind(1,X_train)
Y_train = Y[train_idx]
sigma <- sigma_func(X_train, Y_train)

X_test = X_single[test_idx, , drop = FALSE]
Y_test = Y[test_idx]

a = -rep(1, n)
B = -X_train
c = Y_train

tau_values <- c(0.10, 0.50, 0.95)

plots <- list()

for (i in 1:length(tau_values)) {
  tau <- tau_values[i]
  tau_suffix <- formatC(tau * 100, width = 2, format = "d", flag = "0")
  alpha_0 = tau
  alpha_1 = 1 - tau
  
  #################################### Select lambda ####################################
  
  # Compute full-data solution
  beta <- CVXR::Variable(p)
  beta_0 = CVXR::Variable()
  xi = CVXR::Variable(n)
  eta = CVXR::Variable(n)
  lambda = Parameter(pos = TRUE)
  obj = CVXR::sum_entries(alpha_0 * xi + alpha_1 * eta) + lambda / 2 * CVXR::power(CVXR::norm2(beta), 2)
  constraints = list(xi >= 0, eta >= 0, beta_0 * a + B %*% beta + c + eta >= 0, - beta_0 * a - B %*% beta - c + xi >= 0)
  prob <- CVXR::Problem(Minimize(obj), constraints)
  
  lam_final = 0.01*n
  lam=0.01*n
  
  value(lambda) = lam_final
  result <- CVXR::solve(prob, warm_start = TRUE)
  theta_full = result$getDualValue(constraints[[3]]) - result$getDualValue(constraints[[4]])
  r_train = Y_train - (X_train %*% result$getValue(beta) + result$getValue(beta_0))
  fitted_value_full <- X_train %*% result$getValue(beta) + result$getValue(beta_0)
  lam_final = lam_final / n
  
  #################################### Compute global influence ####################################
  influence_measure = Compute_CaseInflu_nonsmooth(X_train,
                                                  Y_train, lam_final, class = "quantile", tau = tau, influence_measure = "FMD")
  global_influ = influence_measure$global_influence
  
  # Local influence 
  local_influ = Compute_LocalInflu_nonsmooth(X_train,
                                             Y_train, lam_final, class = "quantile", tau = tau)
  
  # local_influ <- sapply(local_influ, function(t) ifelse(t < 1e-5, 0, t))
  
  # LOO influence 
  loo_influ = sapply(1:n, calculate_loo_influence, X_train, y_train, tau=tau,lam=lam,
                     fitted_value_full = fitted_value_full)
  
  
  eps = 1e-3
  Elbow_idx <- which(abs(r_train)<=eps)
  X_elbow_tilde <- X_tilde[Elbow_idx, , drop=FALSE]
  
  generalized_leverage_func <- function(X_tilde, X_elbow_tilde){
    m <- nrow(X_elbow_tilde)
    P_E <- t(X_elbow_tilde) %*% (solve(X_elbow_tilde %*% t(X_elbow_tilde))) %*% X_elbow_tilde
    I_minus_P_E <- diag(dim(P_E)[1]) - P_E
    one_vec <- matrix(1, nrow = m, ncol = 1) 
    Q <- t(one_vec) %*% solve(X_elbow_tilde %*% t(X_elbow_tilde)) %*% X_elbow_tilde
    Q_Qt <- Q %*% t(Q) 
    results_part1 <- numeric(n)
    results_part2 <- numeric(n)
    
    for (i in 1:n) {
      x_i_tilde <- X_tilde[i, , drop = FALSE]
      results_part1[i] <- as.numeric((x_i_tilde) %*% I_minus_P_E %*% t(x_i_tilde))
      Q_i <- Q %*% t(x_i_tilde)
      results_part2[i] <- as.numeric((Q_i - 1) %*% (Q_i - 1) / Q_Qt)
    }
    combined_results <- results_part1 + results_part2
    return(combined_results)
  }
  
  
  generalized_leverage <- generalized_leverage_func(X_tilde = X_tilde, X_elbow_tilde = X_elbow_tilde)
  
  df <- data.frame(X = X_train[,1], y = Y_train, 
                   global_influence = global_influ,
                   global_influence_rescale = global_influ*n/((p+1)*sigma^2),
                   local_influence = local_influ,
                   local_influence_rescale = local_influ*n/((p+1)*sigma^2),
                   loo_influence = loo_influ,
                   loo_influence_rescale = loo_influ*n/((p+1)*sigma^2))
  
  
# top; right; bottom; left
  if(i==1){
  scatter_plot_loo_shade_rescale <- ggplot(df, aes(x = X, y = y)) +
    geom_point(aes(alpha = loo_influence_rescale), color = "blue", size = 2) +  # Change the size of the points
    geom_abline(intercept = result$getValue(beta_0), slope = result$getValue(beta), linetype = "solid") +
    scale_alpha_continuous(range = c(0.1, 1)) + 
    labs(x = 'Sqft_living15', y = 'House price', alpha = 'LOO', title = bquote(tau == .(sprintf("%.2f", tau_values[i])))) + 
    theme_minimal() + 
    theme(
      plot.title = element_text(hjust = 0.5, size = 14), 
      axis.title.x = element_text(size = 14),  
      axis.title.y = element_text(size = 14), 
      legend.title = element_text(size = 10),  
      legend.text = element_text(size = 8),    
      legend.position = "bottom",               
      legend.key.height = unit(0.3, "cm"),     
      legend.key.width = unit(0.3, "cm"),        
      legend.margin = margin(0, 0, 0, 0),
      legend.box = "vertical"
    )
  } else{
    scatter_plot_loo_shade_rescale <- ggplot(df, aes(x = X, y = y)) +
      geom_point(aes(alpha = loo_influence_rescale), color = "blue", size = 2) +  # Change the size of the points
      geom_abline(intercept = result$getValue(beta_0), slope = result$getValue(beta), linetype = "solid") +
      scale_alpha_continuous(range = c(0.1, 1)) + 
      labs(x = 'Sqft_living15', alpha = 'LOO',title = bquote(tau == .(sprintf("%.2f", tau_values[i])))) + 
      theme_minimal() + 
      theme(
        plot.title = element_text(hjust = 0.5, size = 14), 
        axis.title.x = element_text(size = 14),  
        axis.title.y = element_blank(),   
        legend.title = element_text(size = 10),  
        legend.text = element_text(size = 8),    
        legend.position = "bottom",               
        legend.key.height = unit(0.3, "cm"),     
        legend.key.width = unit(0.3, "cm"),        
        legend.margin = margin(0, 0, 0, 0),
        legend.box = "vertical"
      )
  }
  plots[[i]] <- scatter_plot_loo_shade_rescale
}

# save as pdf file 4*12
combined_plots_loo <- grid.arrange(plots[[1]], 
                                   plots[[2]],plots[[3]],nrow = 1)