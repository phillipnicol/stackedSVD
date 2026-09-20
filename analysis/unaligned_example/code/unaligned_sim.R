args <- commandArgs(trailingOnly = FALSE)
file_arg <- "--file="
script_path <- sub(file_arg, "", args[startsWith(args, file_arg)])
script_dir <- if (length(script_path) > 0) dirname(normalizePath(script_path)) else getwd()

data_dir <- file.path(script_dir, "..", "data")
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

library(irlba)

n <- 3000
d <- 3000
c <- 1
c_stack <- 2

theta_values <- c(0.9, 1.15, 1.5, 3)
sin_psi_values <- seq(0, 1, length.out = 32)
psi_values <- asin(sin_psi_values)
n_reps <- 10

stackedSVD <- function(X.list,
                       rank = "auto",
                       max.rank = 50) {
  X.stacked <- do.call(rbind, X.list)
  c.tilde <- nrow(X.stacked) / ncol(X.stacked)
  my.svd <- irlba::irlba(X.stacked, nv = max.rank, nu = max.rank)
  rank <- ifelse(
    rank == "auto",
    sum(my.svd$d > 1 + sqrt(c.tilde)),
    as.numeric(rank)
  )
  my.svd$v[, 1:rank, drop = FALSE]
}

SVDstacked_new <- function(X.list,
                           rank = "auto",
                           max.rank = 50) {
  V.list <- lapply(X.list, function(X) {
    my.svd <- irlba::irlba(X, nv = max.rank, nu = max.rank)
    if (max.rank == 1) {
      my.svd$v <- as.matrix(my.svd$v)
    }
    ci <- nrow(X) / ncol(X)
    rank <- ifelse(
      rank == "auto",
      sum(my.svd$d > 1 + sqrt(ci)),
      as.numeric(rank)
    )
    t(my.svd$v[, 1:max.rank, drop = FALSE])
  })
  V.stacked <- do.call(rbind, V.list)
  my.svd <- svd(V.stacked)
  rank.keep <- ifelse(
    rank == "auto",
    sum(my.svd$d > 1 + 10^-10),
    as.numeric(rank)
  )
  my.svd$v[, 1:rank.keep, drop = FALSE]
}

stackedsvd_theory <- function(theta, psi) {
  theta_adj <- c(theta * sqrt(1 + sin(psi)), theta * sqrt(1 - sin(psi)))
  sum(ifelse(
    theta_adj^4 > c_stack,
    (theta_adj^4 - c_stack) / (theta_adj^4 + theta_adj^2),
    0
  ))
}

svdstack_theory <- function(theta, psi) {
  beta_sq <- ifelse(theta^4 > c, (theta^4 - c) / (theta^4 + theta^2), 0)
  beta_sq * (1 + sin(psi)) / (1 + sin(psi) * beta_sq) +
    beta_sq * (1 - sin(psi)) / (1 - sin(psi) * beta_sq)
}

run_theta_simulation <- function(theta) {
  result_list <- vector("list", length(psi_values) * n_reps)
  result_i <- 1

  for (psi in psi_values) {
    message("theta = ", theta, ", psi = ", signif(psi, 4))

    for (rep in seq_len(n_reps)) {
      V <- qr.Q(qr(matrix(rnorm(d * 2), d, 2)))

      u1 <- rnorm(n)
      u2 <- rnorm(n)
      u1 <- u1 / sqrt(sum(u1^2))
      u2 <- u2 / sqrt(sum(u2^2))

      V1 <- V[, 1]
      X1 <- theta * u1 %*% t(V1) +
        matrix(rnorm(n * d, sd = 1 / sqrt(d)), n, d)

      V2 <- sin(psi) * V[, 1] + cos(psi) * V[, 2]
      X2 <- theta * u2 %*% t(V2) +
        matrix(rnorm(n * d, sd = 1 / sqrt(d)), n, d)

      V.hat.stacksvd <- stackedSVD(list(X1, X2), rank = 2, max.rank = 2)
      stackedsvd_result <- sum((t(V.hat.stacksvd) %*% V)^2)

      V.hat.svdstack <- SVDstacked_new(list(X1, X2), rank = 2, max.rank = 1)
      svdstack_result <- sum((t(V.hat.svdstack) %*% V)^2)

      result_list[[result_i]] <- data.frame(
        theta = theta,
        psi = psi,
        rep = rep,
        stackedSVD_result = stackedsvd_result,
        SVDstacked_result = svdstack_result,
        stackedSVD_theoretical = stackedsvd_theory(theta, psi),
        SVDstacked_theoretical = svdstack_theory(theta, psi)
      )
      result_i <- result_i + 1
    }
  }

  results <- do.call(rbind, result_list)
  output_file <- file.path(
    data_dir,
    paste0("unaligned_simulation_results_theta_", theta, ".csv")
  )
  write.csv(results, file = output_file, row.names = FALSE)
  results
}

set.seed(123)
all_results <- do.call(rbind, lapply(theta_values, run_theta_simulation))

write.csv(
  all_results,
  file = file.path(data_dir, "unaligned_simulation_results_all_thetas.csv"),
  row.names = FALSE
)

print(all_results)
