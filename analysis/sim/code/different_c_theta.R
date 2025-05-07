


M <- 5

c <- rep(1, M)

res <- matrix(0, nrow=100, ncol=4)

d <- 1000

for(i in 1:100) {
  X.list <- list()
  theta <- 1 + rexp(n=5, rate=10)
  v <- rnorm(n=d,sd=1/sqrt(d))
  v <- v/sqrt(sum(v^2))
  for(m in 1:5) {
    u <- rnorm(n=d,sd=1/sqrt(d))
    X.list[[m]] <- theta[m]*outer(u,v) + matrix(rnorm(d^2, sd=1/sqrt(d)), nrow=d,ncol=d)
  }

  res[i,1] <- sum(stackedSVDOracle(X.list,theta.true=theta)*v)^2
  res[i,2] <- sum(SVDstackedOracle(X.list,theta.true=theta)*v)^2
  res[i,3] <- sum(stackedSVDOracle(X.list,theta.true=theta, weight=TRUE)*v)^2
  res[i,4] <- sum(SVDstackedOracle(X.list,theta.true=theta, weight=TRUE)*v)^2

  res[i,] <- res[i,]-getTheoryPred(theta.true=theta, c=c)
  print(res[i,])
}

colnames(res) <- c("Stack SVD", "SVD Stack", "Weighted Stack SVD", "Weighted SVD Stack")

library(reshape2)

df <- melt(res)


p <- df |> ggplot(aes(x=Var2, y = value)) +
  geom_boxplot() +
  geom_hline(yintercept = 0, color="blue",linetype="dashed") +
  theme_bw()

