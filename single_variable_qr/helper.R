# required packages 
pkgs <- c(
  "CVXR",
  "ggplot2",
  "gridExtra",
  "quantreg"
)

for(pkg in pkgs) {
  if(!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

# Install NonsmoothPath from github repo
if(!require("NonsmoothPath", character.only = TRUE)) {
  devtools::install_github("qiuyu1995/NonsmoothPath")
  library(NonsmoothPath)
}

# Install Caseinflu from github repo
if(!require("Caseinflu", character.only = TRUE)) {
  devtools::install_github("qiuyu1995/Caseinflu")
  library(Caseinflu)
}


quantile_obj <- function(W, X, Y, beta, beta_0, lam, tau) {
  r = Y - (X %*% beta + beta_0)
  return(sum(W * (0.5*abs(r) + (tau - 0.5)*r)) + lam / 2 * sum(beta ^ 2) )
}

check_loss <- function(r, tau) {
  return (tau * max(r, 0) + (1 - tau) * max(-r, 0))
}


# Linear interpolation to compute the solution 
# at any w between (0, 1)
find_w_index <- function(w, w_vec){
  if (w > 1 || w < 0)
    return -1
  
  for (i in 1:length(w_vec))
    if (w >= w_vec[i])
      return(i)
}

solution_at_w <- function(w, w_vec, beta_0_vec, beta_mat) {
  # # Case when 0 < w < 1e-6
  # if (w >= 0 && w < 1e-6) {
  #   n_bkpt = length(w_vec)
  #   return(list(beta_0_vec[n_bkpt], beta_mat[, n_bkpt]))
  # }
  
  index = find_w_index(w, w_vec)
  if (index == -1) {
    print("w is between 0 and 1")
    return(0)
  }
  
  # Case when w = 1
  if (index == 1)
    return(list(beta_0_vec[1], beta_mat[, 1]))
  
  # Case when w between w_vec[1] and w_vec[2] and |E_init| = 1
  if (index == 2 && all(beta_mat[,1] == beta_mat[,2])) {
    return(list(beta_0_vec[1], beta_mat[, 1]))
  }
  
  w_m1 = w_vec[index]
  w_m = w_vec[index-1]
  beta_0_m1 = beta_0_vec[index]
  beta_0_m = beta_0_vec[index-1]
  beta_m1 = beta_mat[,index]
  beta_m = beta_mat[,index-1]
  
  beta_0_w = (w - w_m) / (w_m1 - w_m) * beta_0_m1 + (w_m1 - w) / (w_m1 - w_m) * beta_0_m
  beta_w = (w - w_m) / (w_m1 - w_m) * beta_m1 + (w_m1 - w) / (w_m1 - w_m) * beta_m
  return(list(beta_0_w, beta_w))
}

calculate_loo_solution <- function(index,X_train, y_train,tau=tau,lam = lam) {
  quantile_path <- nonsmooth_path(X_train, Y_train, lam, index, class = "quantile", tau = tau)
  w_vec <- quantile_path$W_vec
  beta_0_vec <- quantile_path$Beta_0
  beta_mat <- quantile_path$Beta
  solutions <- solution_at_w(0, w_vec, beta_0_vec, beta_mat)
  beta_0 <- solutions[[1]]
  beta <- solutions[[2]]
  return(list(beta = beta, beta_0 = beta_0))
}

calculate_loo_influence <- function(idx, X_train, y_train, tau = tau, lam=lam, fitted_value_full) {
  loo_solution <- calculate_loo_solution(idx, X_train,y_train, tau = tau,lam=lam)
  fitted_value_loo <- loo_solution$beta_0 + X_train%*%loo_solution$beta
  
  # Calculate the squared difference (influence)
  loo_influence <- mean((fitted_value_full - fitted_value_loo)^2)
  
  return(loo_influence)
}

case_weight_cooks_dist <- function(X, Y, beta, beta_0, beta_w, beta_0_w) {
  r_w = X %*% beta + beta_0 - X %*% beta_w - beta_0_w
  return(mean(r_w^2))
}

sigma_func <- function(X,Y){
  median_reg <- rq(Y ~ X, tau = 0.5)
  
  residuals <- Y - fitted(median_reg)
  sigma <- mad(residuals)
  return(sigma)
}

