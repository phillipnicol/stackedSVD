library(irlba)
library(ggplot2)
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


estimate_theta2 <- function(theta1, theta2, n, d) {

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

  # Numerical safety
  if (!is.finite(beta1_sq) || beta1_sq <= 0) {
    return(0)
  }

  beta1.hat <- sqrt(beta1_sq)

  # Estimate theta2
  my.norm <- sum((X2 %*% my.svd$v[, 1])^2)

  if (my.norm <= c2) {

    theta2.hat <- 0

  } else {

    theta2.hat <- beta1.hat^(-1) * sqrt(my.norm - c2)

    if (!is.finite(theta2.hat)) {
      theta2.hat <- 0
    }
  }

  return(theta2.hat)
}

# Run simulations
results <- expand.grid(
  theta1 = theta1_grid,
  theta2 = theta2_grid,
  rep = 1:B
)

results$theta2_hat <- NA_real_

for (i in seq_len(nrow(results))) {

  results$theta2_hat[i] <- estimate_theta2(
    theta1 = results$theta1[i],
    theta2 = results$theta2[i],
    n = n,
    d = d
  )
}



# Compute MSE
results$sq_error <- (results$theta2_hat - results$theta2)^2

#Save RDS
saveRDS(results, file = "../data/theta2_estimation_results.RDS")

mse_results <- aggregate(
  sq_error ~ theta1 + theta2,
  data = results,
  FUN = mean
)

names(mse_results)[3] <- "mse"

library(dplyr)

# Summary statistics for error bars
plot_data <- results %>%
  group_by(theta1, theta2) %>%
  summarise(
    mse = mean(sq_error),
    se = sd(sq_error) / sqrt(n()),
    .groups = "drop"
  )

ggplot(plot_data,
       aes(x = theta1,
           y = mse,
           color = factor(theta2),
           fill = factor(theta2))) +

  geom_line(linewidth = 1.2) +

  geom_point(size = 2.5) +

  geom_errorbar(
    aes(ymin = mse - 1.96 * se,
        ymax = mse + 1.96 * se),
    width = 0.015,
    alpha = 0.6
  ) +

  scale_color_manual(
    values = c(
      "#0072B2",  # blue
      "#D55E00",  # vermillion
      "#009E73"   # green
    ),
    labels = c("0.1", "0.5", "0.9")
  ) +

  scale_fill_manual(
    values = c(
      "#0072B2",
      "#D55E00",
      "#009E73"
    ),
    labels = c("0.1", "0.5", "0.9")
  ) +

  theme_bw(base_size = 14) +

  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    ),
    legend.position = "top",
    legend.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  ) +

  labs(
    x = expression(theta[1]),
    y = expression(MSE(hat(theta)[2])),
    color = expression(theta[2]),
    fill = expression(theta[2]),
    title = expression(
      paste("MSE of ", hat(theta)[2], " vs ", theta[1])
    )
  )
