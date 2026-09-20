library(irlba)
library(here)

set.seed(1)

setwd(here("analysis", "theta_estimator", "code"))

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Set defaults if arguments not provided
n <- if (length(args) > 0) as.numeric(args[1]) else 1000
d <- if (length(args) > 1) as.numeric(args[2]) else 1000
B <- if (length(args) > 2) as.numeric(args[3]) else 25


theta2_grid <- c(0.1, 0.5, 0.9)
theta1_grid <- seq(1.025, 1.5, by = 0.025)


run_iteration <- function(theta1, theta2, n, d) {

  u1 <- rnorm(n)
  u1 <- u1 / sqrt(sum(u1^2))

  u2 <- rnorm(n)
  u2 <- u2 / sqrt(sum(u2^2))

  v <- rnorm(d)
  v <- v / sqrt(sum(v^2))

  X1 <- theta1 * u1 %*% t(v) +
    matrix(rnorm(n * d, sd = 1 / sqrt(d)), n, d)

  X2 <- theta2 * u2 %*% t(v) +
    matrix(rnorm(n * d, sd = 1 / sqrt(d)), n, d)

  my.svd <- irlba(X1, nv = 1)

  sigma12 <- my.svd$d[1]^2

  c1 <- n / d
  c2 <- n / d

  # Estimate theta1
  if (sigma12 < (1 + sqrt(c1))^2) {

    theta1.hat <- c1

  } else {

    theta1.hat <- sqrt(
      sigma12 - (1 + c1) +
        sqrt((sigma12 - (1 + c1))^2 - 4 * c1)
    ) / sqrt(2)
  }

  # Estimate beta1
  beta1_sq <- (theta1.hat^4 - c1) /
    (theta1.hat^2 * (theta1.hat^2 + 1))

  beta1.hat <- if (is.finite(beta1_sq) && beta1_sq > 0) sqrt(beta1_sq) else NA_real_

  # Estimate theta2
  my.norm <- sum((X2 %*% my.svd$v[, 1])^2)

  if (!is.finite(beta1.hat) || my.norm <= c2) {

    theta2.hat <- 0

  } else {

    theta2.hat <- beta1.hat^(-1) * sqrt(my.norm - c2)

    if (!is.finite(theta2.hat)) {
      theta2.hat <- 0
    }
  }

  # Singular vector recovery
  X <- rbind(X1, X2)
  v.hat <- irlba(X, nv = 1)$v[, 1]
  unweighted_stack_svd <- sum(v.hat * v)^2

  w.hat <- sqrt(c(
    theta1.hat^2 / (theta1.hat^2 + 1),
    theta2.hat^2 / (theta2.hat^2 + 1)
  ))
  X <- rbind(w.hat[1] * X1, w.hat[2] * X2)
  v.hat <- irlba(X, nv = 1)$v[, 1]
  estimated_weight_stack_svd <- sum(v.hat * v)^2

  w.opt <- sqrt(c(
    theta1^2 / (theta1^2 + 1),
    theta2^2 / (theta2^2 + 1)
  ))
  X <- rbind(w.opt[1] * X1, w.opt[2] * X2)
  v.hat <- irlba(X, nv = 1)$v[, 1]
  optimal_weight_stack_svd <- sum(v.hat * v)^2

  return(c(
    theta2_hat = theta2.hat,
    unweighted_stack_svd = unweighted_stack_svd,
    estimated_weight_stack_svd = estimated_weight_stack_svd,
    optimal_weight_stack_svd = optimal_weight_stack_svd
  ))
}

# Run simulations
results <- expand.grid(
  theta1 = theta1_grid,
  theta2 = theta2_grid,
  rep = 1:B
)

results$theta2_hat <- NA_real_
results$unweighted_stack_svd <- NA_real_
results$estimated_weight_stack_svd <- NA_real_
results$optimal_weight_stack_svd <- NA_real_

for (i in seq_len(nrow(results))) {
  if (i %% 50 == 0 || i == 1 || i == nrow(results)) {
    message("Running simulation ", i, " of ", nrow(results))
  }

  sim_results <- run_iteration(
    theta1 = results$theta1[i],
    theta2 = results$theta2[i],
    n = n,
    d = d
  )

  results$theta2_hat[i] <- sim_results["theta2_hat"]
  results$unweighted_stack_svd[i] <- sim_results["unweighted_stack_svd"]
  results$estimated_weight_stack_svd[i] <- sim_results["estimated_weight_stack_svd"]
  results$optimal_weight_stack_svd[i] <- sim_results["optimal_weight_stack_svd"]
}



# Compute MSE
results$sq_error <- (results$theta2_hat - results$theta2)^2

# Save RDS
saveRDS(results, file = "../data/theta2_estimation_results.RDS")
