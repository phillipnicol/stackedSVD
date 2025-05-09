source("../../../R/oracle.R")
source("../../../R/theory_pred.R")

runSim <- function(theta.mean, theta.sd) {

  M <- 5

  res <- matrix(0, nrow=100, ncol=4)

  d <- 1000

  Theta.save <- matrix(0, nrow=100, ncol=M)
  C.save <- matrix(0, nrow=100, ncol=M)

  for(i in 1:100) {
    X.list <- list()

    while(TRUE) {
      c <- rexp(n=M, rate=1) + 0.1
      theta <- c*exp(rnorm(M,mean=log(theta.mean),sd=theta.sd))
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

  return(df)
}


theta.mean <- c(1.1, 1.5, 2)
theta.sd <- c(0.1, 0.25, 1)

params <- expand.grid(theta.mean, theta.sd)

for(i in 1:nrow(params)) {
  df <- runSim(params[i,1], params[i,2])
  saveRDS(df,paste0("../data/ctheta_sim_", params[i,1], "_", params[i,2], ".RDS"))
}






