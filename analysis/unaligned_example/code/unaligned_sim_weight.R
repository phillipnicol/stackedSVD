
library(irlba)
n <- 3000
d <- 3000
c <- 1
c_stack <- 2


stackedSVD <- function(X.list,
                       rank = "auto",
                       max.rank = 50) {
  X.stacked <- do.call(rbind, X.list)
  c.tilde <- nrow(X.stacked)/ncol(X.stacked)
  my.svd <- irlba::irlba(X.stacked, nv=max.rank, nu=max.rank)
  rank <- ifelse(rank == "auto", sum(my.svd$d > 1 + sqrt(c.tilde)),
                 as.numeric(rank))
  V.hat <- my.svd$v[,1:rank]
  return(V.hat)
}

SVDstacked_new <- function(X.list,
                       rank = "auto",
                       max.rank = 50) {
  V.list <- lapply(X.list, function(X) {
    my.svd <- irlba::irlba(X, nv=max.rank, nu=max.rank)
    if(max.rank == 1) {
      my.svd$v <- as.matrix(my.svd$v)
    }
    ci <- nrow(X)/ncol(X)
    rank <- ifelse(rank == "auto", sum(my.svd$d > 1 + sqrt(ci)),
                   as.numeric(rank))
    return(t(my.svd$v[,1:max.rank]))
  })
  V.stacked <- do.call(rbind, V.list)
  my.svd <- svd(V.stacked)
  rank.keep <- ifelse(rank == "auto",
                      sum(my.svd$d > 1 + 10^{-10}),
                      as.numeric(rank))
  return(my.svd$v[,1:rank.keep])
}

#Generate random orthonormal V of dimension d x 2
V <- qr.Q(qr(matrix(rnorm(d*2), d, 2)))

#Generate u1

u1 <- rnorm(n)
u2 <- rnorm(n)

u1 <- u1 / sqrt(sum(u1^2))
u2 <- u2 / sqrt(sum(u2^2))

theta1 <- 1.25
theta2 <- 3.25
#theta <- 1.75
psi <- pi/16
beta <- sqrt((theta^4 - c)/(theta^4 + theta^2))

V1 <- V[,1]
X1 <- theta1 * u1 %*% t(V1) + matrix(rnorm(n*d, sd=1/sqrt(d)), n, d)

V2 <- sin(psi) * V[,1] + cos(psi) * V[,2]
X2 <- theta2 * u2 %*% t(V2) + matrix(rnorm(n*d, sd=1/sqrt(d)), n, d)

# Run simulations for multiple psi values
psi_values <- seq(pi/16, pi/2, by=pi/16)
n_reps <- 10

# Initialize results data frame
results <- data.frame(
  psi = numeric(),
  rep = numeric(),
  stackedSVD_result = numeric(),
  SVDstacked_result = numeric(),
  stackedSVD_theoretical = numeric(),
  SVDstacked_theoretical = numeric()
)

set.seed(123) # For reproducibility

for (psi_val in psi_values) {
  print(psi_val)
  psi <- psi_val
  beta <- sqrt((theta^4 - c)/(theta^4 + theta^2))

  for (rep in 1:n_reps) {
    # Generate random orthonormal V of dimension d x 2
    V <- qr.Q(qr(matrix(rnorm(d*2), d, 2)))

    # Generate u1 and u2
    u1 <- rnorm(n)
    u2 <- rnorm(n)

    u1 <- u1 / sqrt(sum(u1^2))
    u2 <- u2 / sqrt(sum(u2^2))

    V1 <- V[,1]
    X1 <- theta1 * u1 %*% t(V1) + matrix(rnorm(n*d, sd=1/sqrt(d)), n, d)

    V2 <- sin(psi) * V[,1] + cos(psi) * V[,2]
    X2 <- theta2 * u2 %*% t(V2) + matrix(rnorm(n*d, sd=1/sqrt(d)), n, d)

    # StackedSVD method
    V.hat.stacksvd <- stackedSVD(list(X1, X2), rank=2, max.rank = 2)
    stackedsvd_result <- sum((t(V.hat.stacksvd) %*% V)^2)

    # Theoretical prediction for StackedSVD
    theta_adj <- c(theta*sqrt(1+sin(psi)), theta*sqrt(1-sin(psi)))
    stackedsvd_theoretical <- ifelse(theta_adj^4 > c_stack, (theta_adj^4 - c_stack)/(theta_adj^4 + theta_adj^2), 0) |> sum()

    # SVDstacked method
    V.hat.svdstack <- SVDstacked_new(list(X1, X2), rank=2, max.rank = 1)
    svdstack_result <- sum((t(V.hat.svdstack) %*% V)^2)

    # Theoretical prediction for SVDstacked
    svdstack_theoretical <- beta^2*(1+sin(psi))/(1+sin(psi)*beta^2) + beta^2*(1-sin(psi))/(1-sin(psi)*beta^2)


    #Weighted Stack SVD
    v1 <- irlba::irlba(X1, nv=1, nu=1)$v
    v2 <- irlba::irlba(X2, nv=1, nu=1)$v
    V.tilde <- t(cbind(v1, v2))

    #Here we will use the oracle

    #Estimate beta
    b1.hat <- sqrt((theta1^4 - c)/(theta1^4 + theta1^2))
    b2.hat <- sqrt((theta2^4 - c)/(theta2^4 + theta2^2))
    Ab.hat <- matrix(c(1, b1.hat*b2.hat*sin(psi), b1.hat*b2.hat*sin(psi), 1), nrow=2, ncol=2)

    #Ab inv
    Ab.hat.sqrt.inv <- eigen(Ab.hat)$vectors %*% diag(eigen(Ab.hat)$values^{-1/2}) %*% t(eigen(Ab.hat)$vectors)

    W.hat <- t(eigen(Ab.hat.sqrt.inv %*% (Ab.hat - diag(c(b1.hat, b2.hat))) %*% Ab.hat.sqrt.inv)$vectors) %*% Ab.hat.sqrt.inv

    Vw.hat <- svd(W.hat %*% V.tilde)$v
    weighted_svdstack_result <- sum((t(Vw.hat) %*% V)^2)




    # Store results
    results <- rbind(results, data.frame(
      psi = psi,
      rep = rep,
      stackedSVD_result = stackedsvd_result,
      SVDstacked_result = svdstack_result,
      stackedSVD_theoretical = stackedsvd_theoretical,
      SVDstacked_theoretical = svdstack_theoretical,
      weighted_svdstack_result = weighted_svdstack_result
    ))
  }
}

# Save results to file
write.csv(results, file = "../data/unaligned_simulation_results_theta_1.75_weight.csv", row.names = FALSE)

results <- read.csv("../data/unaligned_simulation_results_theta_1.75_weight.csv")
# Print summary
print(results)

# Install/load ggplot2 if needed
library(ggplot2)
library(dplyr)

# Calculate summary statistics
summary_stats <- results %>%
  group_by(psi) %>%
  summarise(
    stackedSVD_mean = mean(stackedSVD_result),
    stackedSVD_sd = sd(stackedSVD_result),
    SVDstacked_mean = mean(SVDstacked_result),
    SVDstacked_sd = sd(SVDstacked_result),
    stackedSVD_theoretical = first(stackedSVD_theoretical),
    SVDstacked_theoretical = first(SVDstacked_theoretical),
    SVDstacked_weighted_mean = mean(weighted_svdstack_result),
    SVDstacked_weighted_sd = sd(weighted_svdstack_result),
    .groups = 'drop'
  ) %>%
  mutate(psi_deg = psi * 180 / pi)  # For labeling

# Reshape data for plotting empirical results
empirical_data <- summary_stats %>%
  pivot_longer(
    cols = c(stackedSVD_mean, SVDstacked_mean),
    names_to = "method",
    values_to = "mean_result"
  ) %>%
  pivot_longer(
    cols = c(stackedSVD_sd, SVDstacked_sd),
    names_to = "method_sd",
    values_to = "sd_result"
  ) %>%
  filter(
    (method == "stackedSVD_mean" & method_sd == "stackedSVD_sd") |
    (method == "SVDstacked_mean" & method_sd == "SVDstacked_sd") |
    (method == "SVDstacked_weighted_mean" & method_sd == "SVDstacked_weighted_sd")
  ) %>%
  select(psi, psi_deg, method, mean_result, sd_result) %>%
  mutate(method = ifelse(method == "stackedSVD_mean", "Stacked SVD", "SVD Stacked"))

# Reshape data for theoretical predictions
theoretical_data <- summary_stats %>%
  pivot_longer(
    cols = c(stackedSVD_theoretical, SVDstacked_theoretical),
    names_to = "method",
    values_to = "theoretical_value"
  ) %>%
  select(psi, psi_deg, method, theoretical_value) %>%
  mutate(method = ifelse(method == "stackedSVD_theoretical", "Stacked SVD", "SVD Stacked"))

# Create publication-quality plot
p <- ggplot() +
  # Empirical results with error bars
  geom_errorbar(
    data = empirical_data,
    aes(x = psi, ymin = mean_result - sd_result, ymax = mean_result + sd_result, color = method),
    width = 0.08,
    size = 0.6,
    alpha = 0.7
  ) +
  geom_point(
    data = empirical_data,
    aes(x = psi, y = mean_result, color = method),
    size = 3,
    alpha = 0.9
  ) +
  # Theoretical predictions as dashed lines
  geom_line(
    data = theoretical_data,
    aes(x = psi, y = theoretical_value, linetype = method),
    size = 1,
    color = "black",
    alpha = 0.6
  ) +
  scale_linetype_manual(
    name = "Theoretical",
    values = c("Stacked SVD" = "dashed", "SVD Stacked" = "dashed"),
    guide = "none"
  ) +
  scale_color_manual(
    name = "Method",
    values = c("Stacked SVD" = "#1b9e77", "SVD Stacked" = "#d95f02")
  ) +
  labs(
    x = expression(paste("Angle (", psi, ") [radians]")),
    y = "Squared Frobenius Norm",
    title = "Unaligned Simulation: Method Performance Comparison"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 11),
    legend.position = "right",
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10),
    panel.grid.major = element_line(color = "gray90", size = 0.3),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5)
  )

# Display and save plot
print(p)
ggsave("../plots/unaligned_simulation_comparison_theta_1.75.pdf", p, width = 10, height = 6, dpi = 300)
