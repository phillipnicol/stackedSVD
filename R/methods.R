

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

