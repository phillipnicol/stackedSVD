source("../../../R/oracle.R")
source("../../../R/theory_pred.R")

runSim <- function(theta.mean, theta.sd, d = 1000, n_iter = 100, n_repeats = 10) {
  M <- 5

  res <- matrix(0, nrow = n_iter, ncol = 4)
  Theta.save <- matrix(0, nrow = n_iter, ncol = M)
  C.save <- matrix(0, nrow = n_iter, ncol = M)

  for (i in 1:n_iter) {
    X.list <- list()

    while (TRUE) {
      c <- rexp(n = M, rate = 1) + 0.1
      theta <- c * exp(rnorm(M, mean = log(theta.mean), sd = theta.sd))
      if (sort(theta^4 / c)[M - 1] < 1) {
        next
      } else {
        break
      }
    }

    # Save theta and c
    Theta.save[i, ] <- theta
    C.save[i, ] <- c

    results_matrix <- matrix(0, nrow = n_repeats, ncol = 4)

    for (r in 1:n_repeats) {
      v <- rnorm(n = d, sd = 1 / sqrt(d))
      v <- v / sqrt(sum(v^2))

      X.list <- list()
      for (m in 1:M) {
        n <- round(d * c[m])
        u <- rnorm(n = n, sd = 1 / sqrt(n))
        X.list[[m]] <- theta[m] * outer(u, v) + matrix(rnorm(n * d, sd = 1 / sqrt(d)), nrow = n, ncol = d)
      }

      results_matrix[r, 1] <- sum(stackedSVDOracle(X.list, theta.true = theta) * v)^2
      results_matrix[r, 2] <- sum(SVDstackedOracle(X.list, theta.true = theta) * v)^2
      results_matrix[r, 3] <- sum(stackedSVDOracle(X.list, theta.true = theta, weight = TRUE) * v)^2
      results_matrix[r, 4] <- sum(SVDstackedOracle(X.list, theta.true = theta, weight = TRUE) * v)^2
    }

    res[i, ] <- colMeans(results_matrix) - getTheoryPred(theta.true = theta, c = c)
    print(res[i, ])
  }

  colnames(res) <- c("Stack SVD", "SVD Stack", "Weighted Stack SVD", "Weighted SVD Stack")

  library(reshape2)
  df <- melt(res)

  return(df)
}


theta.mean <- c(1.1, 1.5, 2)
theta.sd <- c(0.1, 0.25, 1)

params <- expand.grid(theta.mean, theta.sd)

for(i in 1:nrow(params)) {
  df <- runSim(params[i,1], params[i,2])
  saveRDS(df,paste0("../data/ctheta_sim_", params[i,1], "_", params[i,2], ".RDS"))
}






