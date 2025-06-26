
set.seed(1)

generate_Xi <- function(Y, n_ambient) {
  n <- nrow(Y); d <- ncol(Y)

  p <- colSums(Y); p <- p/sum(p)

  #Lambda <- matrix(n_ambient*p, nrow=n, ncol = d, byrow=TRUE)
  Lambda <- matrix(n_ambient, nrow=n, ncol=d)

  return(Y + matrix(rpois(n=n*d, lambda=Lambda), nrow=n, ncol=d))
}


# 0 and 25 rank = 5
# Increasing M test
set.seed(1)
M.try <- 3:10
n.reps <- 10  # Number of repetitions
Res_all <- array(0, dim = c(5, length(M.try), n.reps))

for(rep in 1:n.reps) {
  cat("Repetition:", rep, "\n")

  for(K in 1:length(M.try)) {
    Ytrain <- readRDS(file="../data/Ytrain.RDS")
    Y <- Ytrain

    d <- ncol(Y)

    ixs <- sample(1:10, size=nrow(Y), replace=TRUE)

    X <- list()
    for(i in 1:10) {
      X[[i]] <- generate_Xi(Y[ixs==i,], n_ambient=ifelse(i == 1, 10, 50))
      X[[i]] <- 2 * sqrt(X[[i]] / d)
      X[[i]] <- scale(X[[i]], center=TRUE, scale=FALSE)
    }

    V.true <- svd.test$v[,1:20]

    X.list <- list()
    for(m in 1:M.try[K]) {
      X.list[[m]] <- X[[m]]
    }

    my.rank <- 5

    V1 <- stackedSVD(X.list, rank = my.rank)
    V2 <- SVDstacked(X.list, rank = my.rank)
    V3 <- SVDstackedWeighted(X.list, rank = my.rank)
    V4 <- stackedSVDWeighted(X.list, rank = my.rank)
    V5 <- irlba::irlba(X.list[[1]], nu=my.rank ,nv=my.rank)$v

    V.true <- V.true[,1:my.rank]

    cor.1 <- sqrt(sum((V.true %*% t(V.true) - V1 %*% t(V1))^2))
    cor.2 <- sqrt(sum((V.true %*% t(V.true) - V2 %*% t(V2))^2))
    cor.3 <- sqrt(sum((V.true %*% t(V.true) - V3 %*% t(V3))^2))
    cor.4 <- sqrt(sum((V.true %*% t(V.true) - V4 %*% solve(t(V4) %*% V4) %*% t(V4))^2))
    cor.5 <- sqrt(sum((V.true %*% t(V.true) - V5 %*% t(V5))^2))

    Res_all[,K,rep] <- c(cor.1, cor.2, cor.3, cor.4, cor.5)
    print(c(cor.1, cor.2, cor.3, cor.4, cor.5))
  }
}

# Compute the mean over the repetitions
Res_mean <- apply(Res_all, c(1,2), mean)

Res_mean <- t(Res_mean)

colnames(Res_mean) <- c("Unweighted Stack-SVD", "Unweighted SVD-Stack", "Weighted SVD-Stack", "Weighted Stack-SVD", "Best single matrix")
rownames(Res_mean) <- M.try

Res_mean <- as.data.frame(Res_mean)
Res_mean$M <- M.try

# Plotting
library(tidyverse)
library(reshape2)

df <- Res_mean |> melt(id.vars="M")
#colnames(df) <- c("Method", "M", "Metric")

p <- df |> ggplot(aes(x=M, y=value, color=variable)) +
  geom_point() + geom_line() + theme_bw() +
  xlab("M") + ylab("Norm metric")

print(p)

saveRDS(df, file="../data/increasing_M_one_low_rest_medium.RDS")
