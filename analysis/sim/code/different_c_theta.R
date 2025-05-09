


M <- 5

res <- matrix(0, nrow=100, ncol=4)

d <- 1000

Theta.save <- matrix(0, nrow=100, ncol=M)
C.save <- matrix(0, nrow=100, ncol=M)

for(i in 1:100) {
  X.list <- list()

  while(TRUE) {
    c <- rexp(n=M, rate=1) + 0.1
    theta <- c*exp(rnorm(M,mean=log(1.1),sd=0.1))
    if(sort(theta^4/c)[M-1] < 1) {
      next
    } else{
      break
    }
  }
  v <- rnorm(n=d,sd=1/sqrt(d))
  v <- v/sqrt(sum(v^2))
  for(m in 1:M) {
    n <- round(d*c[m])
    u <- rnorm(n=n,sd=1/sqrt(n))
    #u <- u/sqrt(sum(u^2))
    X.list[[m]] <- theta[m]*outer(u,v) + matrix(rnorm(n*d, sd=1/sqrt(d)), nrow=n,ncol=d)
  }

  res[i,1] <- sum(stackedSVDOracle(X.list,theta.true=theta)*v)^2
  res[i,2] <- sum(SVDstackedOracle(X.list,theta.true=theta)*v)^2
  res[i,3] <- sum(stackedSVDOracle(X.list,theta.true=theta, weight=TRUE)*v)^2
  res[i,4] <- sum(SVDstackedOracle(X.list,theta.true=theta, weight=TRUE)*v)^2

  Theta.save[i,] <- theta
  C.save[i,] <- c

  res[i,] <- res[i,]-getTheoryPred(theta.true=theta, c=c)
  print(res[i,])
}

colnames(res) <- c("Stack SVD", "SVD Stack", "Weighted Stack SVD", "Weighted SVD Stack")

library(reshape2)

df <- melt(res)

library(tidyverse)

p <- df |> ggplot(aes(x=Var2, y = value)) +
  geom_boxplot() +
  geom_hline(yintercept = 0, color="blue",linetype="dashed") +
  theme_bw()








M <- 5

res <- matrix(0, nrow=100, ncol=4)

d <- 1000

for(i in 1:100) {
  X.list <- list()
  #theta <- 1 + rexp(n=5, rate=1)
  c <- rexp(n=M, rate=1)
  theta <- c*exp(rnorm(M,mean=0,sd=0.1))
  v <- rnorm(n=d,sd=1/sqrt(d))
  v <- v/sqrt(sum(v^2))
  for(m in 1:M) {
    n <- round(d*c[m])
    u <- rnorm(n=n,sd=1/sqrt(n))
    #u <- u/sqrt(sum(u^2))
    X.list[[m]] <- theta[m]*outer(u,v) + matrix(rnorm(n*d, sd=1/sqrt(d)), nrow=n,ncol=d)
  }

  res[i,1] <- sum(stackedSVD(X.list, rank = 1)*v)^2
  res[i,2] <- sum(SVDstacked(X.list,rank=1)*v)^2
  res[i,3] <- sum(stackedSVDWeighted(X.list,rank=1)*v)^2
  res[i,4] <- sum(SVDstackedWeighted(X.list,rank=1)*v)^2

  res[i,] <- res[i,]-getTheoryPred(theta.true=theta, c=c)
  print(res[i,])
}

colnames(res) <- c("Stack SVD", "SVD Stack", "Weighted Stack SVD", "Weighted SVD Stack")

library(reshape2)

df <- melt(res)

library(tidyverse)

p <- df |> ggplot(aes(x=Var2, y = value)) +
  geom_boxplot() +
  geom_hline(yintercept = 0, color="blue",linetype="dashed") +
  theme_bw()



