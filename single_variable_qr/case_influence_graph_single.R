rm(list=ls())

## related packages are already listed and imported in helper.R
## pkgs <- c("CVXR","ggplot2","gridExtra","quantreg","NonsmoothPath","Caseinflu")
library(grid)

# specify the working dir
setwd("path/to/github/repo") 

# source helper file
source("single_variable_qr/helper.R")

# Load data
load("data/Input.RData")

set.seed(1)
X = as.matrix(X)
n_total = nrow(X)

############## single predictor
X_single = X[,'sqrtSqft_living15', drop=FALSE]
n = 200
p = 1
subsample = sample(n_total, 2*n)
train_idx = subsample[1:n]
test_idx = subsample[(n+1):(2*n)]

X_train = X_single[train_idx, , drop = FALSE]
X_tilde <- cbind(1,X_train)
Y_train = Y[train_idx]

sigma <- sigma_func(X_train, Y_train)

a = -rep(1, n)
B = -X_train
c = Y_train

plot_list_4 <- list()

tau_values <- c(0.10, 0.50, 0.95)
lam = n*0.01

beta <- CVXR::Variable(p)
beta_0 = CVXR::Variable()
xi = CVXR::Variable(n)
eta = CVXR::Variable(n)
lambda = Parameter(pos = TRUE)

y_ranges <- list(c(0, 1.0), c(0, 0.01), c(0, 2.0))

case_w_inf <- function(X_train, Y_train, lam, indices, tau) {
  Number_grid <- 100
  grid_w <- seq(from = 0, to = 1, length.out = Number_grid + 1)
  
  all_distances <- data.frame(grid_w = rep(grid_w, length(indices)), 
                              distance = rep(0, length(grid_w) * length(indices)), 
                              index = rep(indices, each = length(grid_w)))
  
  for (j in seq_along(indices)) {
    index <- indices[j]
    quantile_path <- nonsmooth_path(X_train, Y_train, lam, index, class = "quantile", tau = tau)
    w_vec <- quantile_path$W_vec
    beta_0_vec <- quantile_path$Beta_0
    beta_mat <- quantile_path$Beta
    beta_0_init <- beta_0_vec[1]
    beta_init <- beta_mat[, 1]
    
    for (k in 1:(Number_grid + 1)) {
      W <- rep(1, n)
      W[index] <- grid_w[k]
      solutions <- solution_at_w(grid_w[k], w_vec, beta_0_vec, beta_mat)
      beta_0_w <- solutions[[1]]
      beta_w <- solutions[[2]]
      distance <- case_weight_cooks_dist(X_train, Y_train, beta_init, beta_0_init, beta_w, beta_0_w)
      
      all_distances$distance[all_distances$index == index & all_distances$grid_w == grid_w[k]] <- distance
    }
  }
  all_distances$distance[all_distances$distance < 1e-6] <- 0
  return(all_distances)
}


for (j in 1:length(tau_values)) {
  all_distances_list <- list()
  tau <- tau_values[j]
  alpha_0 = tau
  alpha_1 = 1 - tau
  
  obj = CVXR::sum_entries(alpha_0 * xi + alpha_1 * eta) + lambda / 2 * CVXR::power(CVXR::norm2(beta), 2)
  constraints = list(xi >= 0, eta >= 0, beta_0 * a + B %*% beta + c + eta >= 0, - beta_0 * a - B %*% beta - c + xi >= 0)
  prob <- CVXR::Problem(Minimize(obj), constraints)
  
  value(lambda) = n*0.01
  result <- CVXR::solve(prob, warm_start = TRUE)
  theta_full = result$getDualValue(constraints[[3]]) - result$getDualValue(constraints[[4]])
  r_train = Y_train - (X_train %*% result$getValue(beta) + result$getValue(beta_0))
  
  rescale_factor = n / ((p + 1) * sigma^2)
  
  eps = 1e-3
  Elbow_idx <- which(abs(r_train)<=eps)
  Left_Elbow_idx <- which(r_train<=-eps)
  Right_Elbow_idx <- which(r_train>=eps)
  
  index_sets <- list(Left_Elbow_idx,Elbow_idx,Right_Elbow_idx)
  
  common_y_range <- y_ranges[[j]]
  y_breaks <- seq(common_y_range[1], common_y_range[2],length.out = 5)
  for (i in 1:3) {
    indices <- index_sets[[i]]
    dist_df <- case_w_inf(X_train, Y_train, lam, indices, tau)
    dist_df$distance_rescale <- dist_df$distance * rescale_factor
    print(range(dist_df$distance_rescale))
    all_distances_list[[i]] <- dist_df
    
    P <- ggplot(dist_df, aes(x = grid_w, y = distance_rescale, color = as.factor(index))) +
      geom_line() +
      scale_y_continuous(limits = common_y_range, breaks = y_breaks, labels = sprintf("%.4f", y_breaks)) + 
      theme_minimal() +
      theme(
        legend.position = "none",
        axis.title = element_blank()
      )
  
    if (i > 1) {
      P <- P + theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())
    }
    if (tau != 0.95) {
      P <- P + theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
    }
    plot_list_4 <- c(plot_list_4, list(P))
  }
}

grid_arrange_shared_legend_text <- function(plots, ncol, nrow, x_label, y_label) {
  grid.newpage()
  
  pushViewport(viewport(layout = grid.layout(nrow + 2, ncol + 2,
                                             widths = unit.c(unit(3, "lines"), 
                                                             unit(1.25, "null"), 
                                                             unit(1, "null"), unit(1, "null"), 
                                                             unit(3, "lines")), 
                                             heights = unit.c(unit(3, "lines"), rep(unit(1, "null"), nrow), unit(3, "lines")))))  # Space for column labels
  

  column_labels <- c("Left Elbow", "Elbow", "Right Elbow")
  for (col in 1:ncol) {
    grid.text(column_labels[col], vp = viewport(layout.pos.row = 1, layout.pos.col = col + 1), 
              gp = gpar(fontsize = 12, col = "black"))
  }
  
  row_labels <- c("Quantile: 0.10", "Quantile: 0.50", "Quantile: 0.95")
  for (row in 1:nrow) {
    grid.text(row_labels[row], vp = viewport(layout.pos.row = row + 1, layout.pos.col = ncol + 2), 
              gp = gpar(fontsize = 12, col = "black"), rot = 270)
  }
  
  for (i in 1:length(plots)) {
    row <- ceiling(i / ncol)
    col <- (i - 1) %% ncol + 1
    print(plots[[i]], vp = viewport(layout.pos.row = row + 1, layout.pos.col = col + 1))
  }
  
  grid.text(x_label, y = unit(1.5, "lines"), just = "center", gp = gpar(fontsize = 14))  
  grid.text(y_label, x = unit(1.5, "lines"), just = "center", rot = 90, gp = gpar(fontsize = 14)) 
}

grid_arrange_shared_legend_text(plot_list_4, ncol = 3, nrow = 3, 
                                x_label = expression("Case Weight " * omega), 
                                y_label = expression("Scaled Cook's Distance"))

# single predictor plot size: 6.1 * 6