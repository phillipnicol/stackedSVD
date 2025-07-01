#' Stacked SVD
#'
#' Performs truncated SVD on the row-stacked version of a list of matrices.
#'
#' @param X.list A list of numeric matrices with matching column dimensions.
#' @param rank Desired rank or "auto" for automatic selection.
#' @param max.rank Maximum allowable rank for truncated SVD.
#'
#' @return A matrix of right singular vectors.
#' @export
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

#' SVD Stacked
#'
#' Computes the right singular vectors by individually applying SVD and stacking the results.
#'
#' @param X.list A list of numeric matrices with matching column dimensions.
#' @param rank Desired rank or "auto" for automatic selection.
#' @param max.rank Maximum allowable rank for truncated SVD.
#'
#' @return A matrix of aggregated right singular vectors.
#' @export
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

#' Weighted Stacked SVD
#'
#' Constructs right singular vectors by weighting and stacking based on estimated signal strength.
#'
#' @param X.list A list of numeric matrices with matching column dimensions.
#' @param rank Desired rank or "auto" for automatic selection.
#' @param max.rank Maximum allowable rank for truncated SVD.
#'
#' @return A matrix of weighted consensus right singular vectors.
#' @export
SVDstackedWeighted <- function(X.list,
                               rank = "auto",
                               max.rank = 50) {
  M <- length(X.list)

  Theta <- lapply(X.list, function(X) {
    my.svd <- irlba::irlba(X, nv=max.rank, nu=max.rank)
    ci <- nrow(X)/ncol(X)
    sapply(my.svd$d, function(x) ifelse(x > 1 + sqrt(ci), max(Re(polyroot(c(ci,0,ci+1-x^2,0,1)))), 0))
  })

  max.rank <- max(unlist(lapply(Theta, function(theta) sum(theta > 10^{-10}))))
  c <- unlist(lapply(X.list, function(X) nrow(X)/ncol(X)))

  Beta <- matrix(0, nrow=M, ncol=max.rank)
  Wmr <- Beta
  for(r in 1:max.rank) {
    best.m <- which.max(unlist(sapply(Theta, function(theta) theta[r])))
    V.best <- irlba::irlba(X.list[[best.m]], nv=max.rank, nu=max.rank)$v
    best.beta <- sqrt(1 - (c[best.m] + Theta[[best.m]][r]^2)/(Theta[[best.m]][r]^2 + Theta[[best.m]][r]^4))
    for(m in 1:M) {
      if(m == best.m) {
        Beta[m, r] <- best.beta
        theta_mr <- Theta[[best.m]][r]
        Wmr[m,r] <- theta_mr*sqrt((theta_mr^2 + 1)/(theta_mr^2 + c[m]))
      } else{
        theta_mr <- 1/best.beta * (sqrt(sum((X.list[[m]] %*% V.best[,r])^2) - c[m]))
        if(is.nan(theta_mr)) {
          Beta[m,r] <- 0
          Wmr[m,r] <- 0
        } else if(theta_mr > c[m]^{1/4}) {
          Beta[m,r] <- sqrt(1 - (c[m] + theta_mr^2)/(theta_mr^2 + theta_mr^4))
          Wmr[m,r] <- theta_mr*sqrt((theta_mr^2 + 1)/(theta_mr^2 + c[m]))
        } else{
          Beta[m,r] <- 0
          Wmr[m,r] <- theta_mr*sqrt((theta_mr^2 + 1)/(theta_mr^2 + c[m]))
        }
      }
    }
  }

  W <- Wmr
  V.list <- list()
  for(m in 1:M) {
    my.svd <- irlba::irlba(X.list[[m]], nv=max.rank, nu=max.rank)
    V.list[[m]] <- t(my.svd$v) * W[m, ]
  }

  V.stacked <- do.call(rbind, V.list)
  my.svd <- svd(V.stacked)
  rank.keep <- ifelse(rank == "auto",
                      sum(my.svd$d > 1 + 10^{-10}),
                      as.numeric(rank))
  return(my.svd$v[,1:rank.keep])
}

#' SVD Stacked with Weighting
#'
#' Computes right singular vectors using a weighted stacking procedure based on theoretical theta values.
#'
#' @param X.list A list of numeric matrices with matching column dimensions.
#' @param rank Desired number of components to keep.
#' @param max.rank Maximum allowable rank for truncated SVD.
#'
#' @return A matrix of consensus right singular vectors.
#' @export
stackedSVDWeighted <- function(X.list, rank, max.rank = 50) {
  M <- length(X.list)
  d <- ncol(X.list[[1]])

  Theta <- lapply(X.list, function(X) {
    my.svd <- irlba::irlba(X, nv=max.rank, nu=max.rank)
    ci <- nrow(X)/ncol(X)
    sapply(my.svd$d, function(x) ifelse(x > 1 + sqrt(ci), max(Re(polyroot(c(ci,0,ci+1-x^2,0,1)))), 0))
  })

  max.rank <- max(unlist(lapply(Theta, function(theta) sum(theta > 10^{-10}))))
  max.rank <- rank
  c <- unlist(lapply(X.list, function(X) nrow(X)/ncol(X)))

  V.hat <- matrix(0, nrow=d, ncol=max.rank)
  W <- matrix(0, nrow = M, ncol=max.rank)

  for(r in 1:max.rank) {
    best.m <- which.max(unlist(sapply(Theta, function(theta) theta[r])))
    V.best <- irlba::irlba(X.list[[best.m]], nv=max.rank, nu=max.rank)$v
    best.beta <- sqrt(1 - (c[best.m] + Theta[[best.m]][r]^2)/(Theta[[best.m]][r]^2 + Theta[[best.m]][r]^4))

    Xw <- list()
    for(m in 1:M) {
      if(m == best.m) {
        W[m,r] <- Theta[[best.m]][r]/(sqrt(Theta[[best.m]][r]^2 + c[m]^2))
      } else{
        theta_mr <- 1/best.beta * (sqrt(sum((X.list[[m]] %*% V.best[,r])^2) - c[m]))
        if(is.nan(theta_mr)) {
          W[m,r] <- 0
        } else{
          W[m,r] <- theta_mr/(sqrt(theta_mr^2 + c[m]^2))
        }
      }
      Xw[[m]] <- W[m,r]*X.list[[m]]
    }

    Xbigw <- do.call(rbind, Xw)
    try(V.hat[,r] <- irlba::irlba(Xbigw, nv=max.rank)$v[,r])
  }

  return(V.hat)
}
