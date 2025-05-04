
#0 and 25 rank = 5

#IncreasingMtest

M.try <- c(3:10)

Res <- matrix(0, nrow=5, ncol=length(M.try))

for(K in 1:length(M.try)) {

  #Divide into 10 blocks


  Ytrain <- readRDS(file="../data/Ytrain.RDS")
  Y <- Ytrain

  #M <- M.try[K]

  ixs <- sample(1:10, size=nrow(Y), replace=TRUE)

  X <- list()
  for(i in 1:10) {
    #X[[i]] <- generate_Xi(Y[ixs==i,], n_ambient=ifelse(i == 1,
    #                                                   0,
    #                                                   25))

    X[[i]] <- generate_Xi(Y[ixs==i,], n_ambient = i*5)

    X[[i]] <- 2*sqrt(X[[i]]/d)
    X[[i]] <- scale(X[[i]], center=TRUE, scale=FALSE)
  }


  #X0 <- sweep(Y, MARGIN=1,STATS=rowSums(Y)/nmed, FUN="/")
  #X0 <- 2*sqrt(Y/d)
  #X0 <- scale(X0, center=TRUE, scale=FALSE)
  #svd.0 <- irlba::irlba(X0, nu=10, nv=10)
  #V.true <- svd.0$v

  V.true <- svd.test$v[,1:20]

  #4
  X.list <- list()
  for(m in 1:M.try[K]) {
    X.list[[m]] <- X[[m]]
  }

  my.rank <- 5

  V1 <- stackedSVD(X.list, rank = my.rank)
  V2 <- SVDstacked(X.list, rank = my.rank)
  V3 <- SVDstackedWeighted(X.list, rank = my.rank)
  V4 <- stackedSVDWeighted(X.list, rank = my.rank)


  #if(ncol(V4) != 10) {
  #  stop("ERROR")
  #}

  V.true <- V.true[,1:my.rank]

  cor.1 <- sqrt(sum((V.true %*% t(V.true) - V1 %*% t(V1))^2))
  cor.2 <- sqrt(sum((V.true %*% t(V.true) - V2 %*% t(V2))^2))
  cor.3 <- sqrt(sum((V.true %*% t(V.true) - V3 %*% t(V3))^2))
  cor.4 <- sqrt(sum((V.true %*% t(V.true) - V4 %*% solve(t(V4) %*% V4) %*% t(V4))^2))

  V5 <- irlba::irlba(X.list[[1]], nu=10 ,nv=10)$v

  cor.5 <-  sqrt(sum((V.true %*% t(V.true) - V5 %*% t(V5))^2))

  Res[,K] <- c(cor.1, cor.2, cor.3, cor.4, cor.5)

  print(Res[,K])
}


rownames(Res) <- c("Stack SVD", "SVD Stack", "SVD Stack (W)", "Stack SVD (W)", "Low noise matrix")
colnames(Res) <- M.try


library(tidyverse)
library(reshape2)

df <- Res |> melt()

colnames(df) <- c("Method", "M", "Metric")

p <- df |> ggplot(aes(x=M, y = Metric, color=Method)) +
  geom_point() + geom_line() + theme_bw() +
  xlab("M") + ylab("Norm metric")























#Increasing nnoise

kappa <- seq(10, 100, by=10)

Res <- matrix(0, nrow=5, ncol=length(kappa))

for(K in 1:length(M.try)) {

  M <- 5
  #M <- M.try[K]

  ixs <- sample(1:M, size=nrow(Y), replace=TRUE)

  X <- list()
  for(i in 1:M) {
    X[[i]] <- generate_Xi(Y[ixs==i,], n_ambient=i*kappa[K])

    X[[i]] <- 2*sqrt(X[[i]]/d)
    X[[i]] <- scale(X[[i]], center=TRUE, scale=FALSE)
  }


  #X0 <- sweep(Y, MARGIN=1,STATS=rowSums(Y)/nmed, FUN="/")
  X0 <- 2*sqrt(Y/d)
  X0 <- scale(X0, center=TRUE, scale=FALSE)
  svd.0 <- irlba::irlba(X0, nu=10, nv=10)
  V.true <- svd.0$v

  X.list <- X

  V1 <- stackedSVD(X.list, rank = 10)
  V2 <- SVDstacked(X.list, rank = 10)
  V3 <- SVDstackedWeighted(X.list, rank = 10)
  V4 <- stackedSVDWeighted(X.list, rank = 10)

  if(ncol(V4) != 10) {
    stop("ERROR")
  }

  cor.1 <- sqrt(sum((V.true %*% t(V.true) - V1 %*% t(V1))^2))
  cor.2 <- sqrt(sum((V.true %*% t(V.true) - V2 %*% t(V2))^2))
  cor.3 <- sqrt(sum((V.true %*% t(V.true) - V3 %*% t(V3))^2))
  cor.4 <- sqrt(sum((V.true %*% t(V.true) - V4 %*% solve(t(V4) %*% V4) %*% t(V4))^2))

  V5 <- irlba::irlba(X.list[[1]], nu=10 ,nv=10)$v

  cor.5 <-  sqrt(sum((V.true %*% t(V.true) - V5 %*% t(V5))^2))

  Res[,K] <- c(cor.1, cor.2, cor.3, cor.4, cor.5)

  print(Res[,K])
}


rownames(Res) <- c("Stack SVD", "SVD Stack", "SVD Stack (W)", "Stack SVD (W)", "Low noise matrix")
colnames(Res) <- M.try


library(tidyverse)
library(reshape2)

df <- Res |> melt()

colnames(df) <- c("Method", "M", "Metric")

p <- df |> ggplot(aes(x=M, y = Metric, color=Method)) +
  geom_point() + geom_line() + theme_bw() +
  xlab("M") + ylab("Norm metric")






