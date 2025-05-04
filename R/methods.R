

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

  rank.keep <- ifelse(rank == "auto",
                      sum(my.svd$d > 1 + 10^{-10}),
                      as.numeric(rank))

  return(my.svd$v[,1:rank.keep])
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

    best.beta <- sqrt(1 - (c[best.m] + Theta[[best.m]][r]^2)/(Theta[[best.m]][r]^2 + Theta[[best.m]][r]^4))
    #print(best.beta)
    #cat("Best beta for rank", r, " is ", best.beta, "\n")

    for(m in 1:M) {
      if(m == best.m) {
        Beta[m, r] <- best.beta
      } else{
        theta_mr <- 1/best.beta * (sqrt(sum((X.list[[m]] %*% V.best[,r])^2) - c[m]))
        if(is.nan(theta_mr)) {
          Beta[m,r] <- 0
        } else if(theta_mr > c[m]^{1/4}) {
          Beta[m,r] <- sqrt(1 - (c[m] + theta_mr^2)/(theta_mr^2 + theta_mr^4))
        } else{
          Beta[m,r] <- 0
        }
      }
    }
  }

  W <- 1/sqrt(1 - Beta^2)

  V.list <- list()
  #cat("M: ", M, "\n")
  for(m in 1:M) {
    my.svd <- irlba::irlba(X.list[[m]], nv=max.rank, nu=max.rank)

    V.list[[m]] <- t(my.svd$v) * W[m, ]
    #print(dim(V.list[[m]]))
  }

  V.stacked <- do.call(rbind, V.list)

  #print(dim(V.stacked))

  my.svd <- svd(V.stacked)

  rank.keep <- ifelse(rank == "auto",
                      sum(my.svd$d > 1 + 10^{-10}),
                      as.numeric(rank))

  return(my.svd$v[,1:rank.keep])
}

stackedSVDWeighted <- function(X.list, rank, max.rank = 50) {
  M <- length(X.list)

  d <- ncol(X.list[[1]])

  #Get Theta
  Theta <- lapply(X.list, function(X) {
    my.svd <- irlba::irlba(X, nv=max.rank, nu=max.rank)

    ci <- nrow(X)/ncol(X)

    sapply(my.svd$d, function(x) ifelse(x > 1 + sqrt(ci), max(Re(polyroot(c(ci,0,ci+1-x^2,0,1)))), 0))
  })

  #Get Rank
  max.rank <- max(unlist(lapply(Theta, function(theta) sum(theta > 10^{-10}))))

  max.rank <- rank

  c <- unlist(lapply(X.list, function(X) nrow(X)/ncol(X)))

  #Get w

  V.hat <- matrix(0, nrow=d, ncol=max.rank)
  W <- matrix(0, nrow = M, ncol=max.rank)
  for(r in 1:max.rank) {
    best.m <- which.max(unlist(sapply(Theta, function(theta) theta[r])))

    V.best <- irlba::irlba(X.list[[best.m]], nv=max.rank, nu=max.rank)$v

    best.beta <- sqrt(1 - (c[best.m] + Theta[[best.m]][r]^2)/(Theta[[best.m]][r]^2 + Theta[[best.m]][r]^4))
    #print(best.beta)
    #cat("Best beta for rank", r, " is ", best.beta, "\n")

    Xw <- list()

    for(m in 1:M) {
      if(m == best.m) {
        #print(Theta[[best.m]][r])
        W[m,r] <- Theta[[best.m]][r]* sqrt((Theta[[best.m]][r]^2 + 1)/(Theta[[best.m]][r]^2 + c[m]))
      } else{
        theta_mr <- 1/best.beta * (sqrt(sum((X.list[[m]] %*% V.best[,r])^2) - c[m]))
        if(is.nan(theta_mr)) {
          W[m,r] <- 0
        } else{
          W[m,r] <- theta_mr * sqrt((theta_mr^2 + 1)/(theta_mr^2 + c[m]))
        }
      }

      if(r == 1) {
        print(W[,m])
      }
      Xw[[m]] <- W[m,r]*X.list[[m]]
    }

    Xbigw <- do.call(rbind, Xw)

    try(V.hat[,r] <- irlba::irlba(Xbigw, nv=max.rank)$v[,r])
  }

  return(V.hat)
}
