

#Test for Tavor
#results_matrix <- matrix(0, nrow = n_repeats, ncol = 4)

M <- 2

n <- 1000
d <- 1000

c <- c(1,1)
theta <- c(0.95, 0.8)

v <- rnorm(n = d, sd = 1)
v <- v / sqrt(sum(v^2))

X.list <- list()
for (m in 1:M) {
  n <- round(d * c[m])
  u <- rnorm(n = n, sd = 1 / sqrt(n))
  X.list[[m]] <- theta[m] * outer(u, v) + matrix(rnorm(n * d, sd = 1 / sqrt(d)), nrow = n, ncol = d)
}

test.me <- rep(1,5)
test.me.weight <- rep(1,5)
for(i in 1:5) {
  v <- rnorm(n = d, sd = 1)
  v <- v / sqrt(sum(v^2))

  X.list <- list()
  for (m in 1:M) {
    n <- round(d * c[m])
    u <- rnorm(n = n, sd = 1 / sqrt(n))
    X.list[[m]] <- theta[m] * outer(u, v) + matrix(rnorm(n * d, sd = 1 / sqrt(d)), nrow = n, ncol = d)
  }

  v.hat <- stackedSVDOracle(X.list, theta.true = theta)
  test.me[i] <- sum(v.hat*v)^2

  v.hat <- stackedSVDOracle(X.list,theta.true=theta,weight=TRUE)
  test.me.weight[i] <- sum(v.hat*v)^2
  print(test.me[i])
  print(test.me.weight[i])
}

getTheoryPred(theta.true = theta, c = c)
