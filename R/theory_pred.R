


compute_x_star <- function(theta.true, c.vec = rep(1, length(theta.true))) {
  # Objective function: the difference from 1
  f <- function(x) {
    sum((theta.true^2) * (1 - x) / (c.vec * theta.true^(-2) + x)) - 1
  }

  # Use uniroot to find the root in (0, 1)
  result <- uniroot(f, interval = c(1e-8, 1 - 1e-8), tol = 1e-10)

  return(result$root)
}

### Get theoretical values

getTheoryPred <- function(theta.true, c) {

  #Stack-SVD

  theta2 <- sqrt(sum(theta.true^2))

  stack.svd <- ifelse(theta2^4/sum(c) > 1, (theta2^4 - sum(c))/(theta2^4 + theta2^2), 0)

  #SVD-Stack

  beta.true <- ifelse(theta.true > c^{1/4},
                      sqrt(1 - (c+theta.true^2)/(theta.true^4+theta.true^2)),
                      0)

  A <- beta.true %*% t(beta.true) + diag(1 - beta.true^2)

  svd.stack <- (sum(beta.true * eigen(A)$vectors[,1]))^2/eigen(A)$values[1]

  #Stack-SVD W

  stack.svd.w <- compute_x_star(theta.true, c.vec=c)

  #SVD-Stack W
  S <- sum(beta.true^2/(1 - beta.true^2))
  svd.stack.w <- S/(S+1)

  vec <- c(stack.svd, svd.stack, stack.svd.w, svd.stack.w)
  names(vec) <- c("Stack SVD", "SVD Stack", "Weighted Stack SVD", "Weighted SVD Stack")
  return(vec)
}
