


stackedSVDOracle <- function(X.list, theta.true,
                             weight=FALSE) {
  M <- length(X.list)

  d <- ncol(X.list[[1]])

  c <- unlist(lapply(X.list, function(X) nrow(X)/ncol(X)))

  print(c)

  if(weight) {
    w <- theta.true/(sqrt(theta.true^2 + c))
  } else{
    w <- rep(1,M)
  }

  print(w)

  for(m in 1:M) {
    X.list[[m]] <- w[m] * X.list[[m]]
  }

  X.stacked <- do.call(rbind, X.list)

  my.svd <- irlba::irlba(X.stacked, nv=1, nu=1)

  return(my.svd$v)
}




SVDstackedOracle <- function(X.list,
                             theta.true,
                             weight=FALSE) {
  M <- length(X.list)

  d <- ncol(X.list[[1]])

  c <- unlist(lapply(X.list, function(X) nrow(X)/ncol(X)))

  beta.true <- lapply(1:M, function(m) {
    if(theta.true[m] > c[m]^{1/4}) {
      sqrt(1 - (c[m] + theta.true[m]^2)/(theta.true[m]^4 + theta.true[m]^2))
    } else{
      0
    }
  }) |> unlist()

  #Get w

  if(weight) {
    #w <- 1/sqrt(1-beta.true^2)
    w <- theta*sqrt((theta^2 + 1)/(theta^2 + c))
  } else{
    w <- rep(1,M)
  }

  V.stacked <- matrix(0, nrow=M, ncol=d)

  for(m in 1:M) {
    my.svd <- irlba::irlba(X.list[[m]], nv = 1, nu = 1)

    V.stacked[m,] <- w[m]*my.svd$v
  }


  my.svd <- svd(V.stacked)

  return(my.svd$v[,1])
}

