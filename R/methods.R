

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


SVDstacked <- function(X.list,
                       rank = "auto",
                       max.rank = 50) {

  V.list <- lapply(X.list, function(X) {
    my.svd <- irlba::irlba(X, nv=max.rank, nu=max.rank)

    ci <- nrow(X)/ncol(X)

    rank <- ifelse(rank == "auto", sum(my.svd$d > 1 + sqrt(ci)),
                   as.numeric(rank))

    return(t(my.svd$v[,1:rank]))
  })

  V.stacked <- do.call(rbind, V.list)

  my.svd <- svd(V.stacked)

  return(my.svd$v)
}


SVDstackedWeighted <- function(X.list,
                       rank = "auto",
                       max.rank = 50) {
  M <- length(X.list)

  #Get Theta
  Theta <- lapply(X.list, function(X) {
    my.svd <- irlba::irlba(X, nv=max.rank, nu=max.rank)

    ci <- nrow(X)/ncol(X)

    sapply(my.svd$d, function(x) ifelse(x > 1 + sqrt(ci), max(Re(polyroot(c(ci,0,ci+1-x^2,0,1)))), 0))
  })

  #Get Rank
  max.rank <- max(unlist(lapply(Theta, function(theta) sum(theta > 10^{-10}))))

  c <- unlist(lapply(X.list, function(X) nrow(X)/ncol(X)))

  #Get w

  Beta <- matrix(0, nrow=M, ncol=max.rank)
  for(r in 1:max.rank) {
    best.m <- which.max(unlist(sapply(Theta, function(theta) theta[r])))

    V.best <- irlba::irlba(X.list[[best.m]], nv=max.rank, nu=max.rank)$v

    best.beta <- 1 - (c[best.m] + Theta[[best.m]][r]^2)/(Theta[[best.m]][r]^2 + Theta[[best.m]][r]^4)

    for(m in 1:M) {
      if(m == best.m) {
        Beta[m, r] <- best.beta
      } else{
        Beta[m, r] <- min(1/best.beta * (sqrt(sum((X.list[[m]] %*% V.best[,r])^2) - c[m])), 1 - 10^{-10})
      }
    }
  }

  W <- 1/sqrt(1 - Beta^2)

  V.list <- list()
  for(m in 1:M) {
    my.svd <- irlba::irlba(X.list[[m]], nv=max.rank, nu=max.rank)

    ci <- nrow(X.list[[m]])/ncol(X.list[[m]])

    V.list[[m]] <- W[m,] * t(my.svd$v[,1:max.rank])
  }

  V.stacked <- do.call(rbind, V.list)

  my.svd <- svd(V.stacked)

  return(my.svd$v)
}



