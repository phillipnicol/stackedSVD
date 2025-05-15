

library(DuoClustering2018)
library(ggpubr)
library(tidyverse)

set.seed(1)

sce <- DuoClustering2018::sce_filteredHVG10_Zhengmix8eq()
Y <- t(sce@assays@data$counts)

generate_Xi <- function(Y, n_ambient) {
  n <- nrow(Y); d <- ncol(Y)

  p <- colSums(Y); p <- p/sum(p)

  #Lambda <- matrix(n_ambient*p, nrow=n, ncol = d, byrow=TRUE)
  Lambda <- matrix(n_ambient, nrow=n, ncol=d)

  return(Y + matrix(rpois(n=n*d, lambda=Lambda), nrow=n, ncol=d))
}

#Y <- t(sce@assays@data$counts)
#Y <- sqrt(sweep(Y, MARGIN=1,STATS=rowSums(Y), FUN="/"))
#Y <- scale(Y,center=TRUE,scale=FALSE)

d <- ncol(Y)

X0 <- 2*sqrt(Y/d)
X0 <- scale(X0, center=TRUE, scale=FALSE)
svd.0 <- irlba::irlba(X0, nu=10, nv=10)
V.true <- svd.0$v

svd.plot <- irlba::irlba(X0, nu=25, nv=25)
topd <- svd.plot$d
c <- nrow(Y)/ncol(Y)
p <- ggplot(data=data.frame(x=1:25, y= topd), aes(x=x,y=y)) +
  geom_point() + theme_bw() +
  geom_hline(yintercept = 1 + sqrt(c), color="blue",linetype="dashed") +
  xlab("") + ylab("Singular value")

ggsave(p, filename="/Users/phillipnicol/Desktop/tenx_sv_plot.png",
       width=6.06, height=4.65, units="in")

X1 <- generate_Xi(Y, n_ambient=0)
X2 <- generate_Xi(Y, n_ambient=100)
X3 <- generate_Xi(Y, n_ambient=500)
Xbig <- rbind(X1, X2, X3)
nmed <- median(rowSums(Xbig))


X1 <- sweep(X1, MARGIN=1,STATS=rowSums(X1)/nmed, FUN="/")
X1 <- 2*sqrt(X1/d)
X1 <- scale(X1, center=TRUE, scale=FALSE)


#X2 <- generate_Xi(Y, n_ambient=200)
X2 <- sweep(X2, MARGIN=1,STATS=rowSums(X2)/nmed, FUN="/")
X2 <- 2*sqrt(X2/d)
X2 <- scale(X2, center=TRUE, scale=FALSE)

X3 <- sweep(X3, MARGIN=1,STATS=rowSums(X3)/nmed, FUN="/")
X3 <- 2*sqrt(X3/d)
X3 <- scale(X3, center=TRUE, scale=FALSE)

#X3 <- rbind(X3, X3, X3, X3, X3)

#n <- nrow(Y); d <- ncol(Y)

#ixs <- sample(c(1:3), size=nrow(Y), replace=TRUE, prob=c(0.1,0.05,0.85))

#X1 <- Y[ixs==1,]; X2 <- Y[ixs==2,]; X3 <- Y[ixs==3,]

## Scramble da X3
#X3 <- matrix(sample(X3, size = nrow(X3)*ncol(X3), replace=TRUE), nrow=nrow(X3),
#             ncol=ncol(X3))

#phenoid <- c(sce$phenoid[ixs == 1], sce$phenoid[ixs == 2], sce$phenoid[ixs == 3])

#svd.1 <- irlba::irlba(X1, nv=50)

#svd.2 <- irlba::irlba(X2, nv=50)

#sigma2 <- sqrt(mean((X2 - svd.2$u %*% diag(svd.2$d) %*% t(svd.2$v))^2))

#X2 <- X2/(sqrt(d)*sigma2)


#svd.3 <- irlba::irlba(X3, nv=50)

#sigma3 <- sqrt(mean((X3 - svd.3$u %*% diag(svd.3$d) %*% t(svd.3$v))^2))

#X3 <- X3/(sqrt(d)*sigma3)


svd.1 <- irlba::irlba(X1,nv=10)

svd.2 <- irlba::irlba(X2,nv=10)

svd.3 <- irlba::irlba(X3,nv=10)



##Stack - SVD

library(tidyverse)

phenoid <- rep(sce$phenoid, 3)

Xbig <- rbind(X1, X2, X3)
Xbig <- scale(Xbig, center=TRUE, scale=FALSE)

stack.svd <- irlba::irlba(Xbig, nv=10)

n <- nrow(X1)

df.stacksvd <- data.frame(x=stack.svd$u[,3],
                          y=stack.svd$u[,4],
                          color=sce$phenoid) |>
  ggplot(aes(x=x,y=y,color=color)) +
  geom_point(size=0.5) +
  theme_bw() +
  ggtitle("Stack-SVD")

V.1 <- stack.svd$v

cor.1 <- sqrt(sum((V.true %*% t(V.true) - V.1 %*% t(V.1))^2))

df.svd1 <- data.frame(x=svd.1$u[,1],
                      y=svd.1$u[,2],
                      color=sce$phenoid) |>
  ggplot(aes(x=x,y=y,color=color)) +
  geom_point(size=0.5) +
  theme_bw() +
  ggtitle("ambient = 0") +
  xlab("V1") + ylab("V2") + labs(color="Cell type")

df.svd2 <- data.frame(x=svd.2$u[,1],
                      y=svd.2$u[,2],
                      color=sce$phenoid) |>
  ggplot(aes(x=x,y=y,color=color)) +
  geom_point(size=0.5) +
  theme_bw() +
  ggtitle("ambient = 100") +
  xlab("V1") + ylab("V2") + labs(color="Cell type")

df.svd3 <- data.frame(x=svd.3$u[,1],
                      y=svd.3$u[,2],
                      color=sce$phenoid) |>
  ggplot(aes(x=x,y=y,color=color)) +
  geom_point(size=0.5) +
  theme_bw() +
  ggtitle("ambient = 1000") +
  xlab("V1") + ylab("V2") + labs(color="Cell type")


p.indiv <- ggarrange(df.svd1, df.svd2, df.svd3, nrow=1, common.legend=TRUE, legend = "top")




df <- readRDS("../data/increasing_M_one_low_rest_medium.RDS")

avg.low <- df |> filter(variable == "Low noise matrix") |> select(value)
avg.low <- mean(avg.low$value)

df[df$variable == "Low noise matrix",]$value <- avg.low

p <- df |> ggplot(aes(x=M, y=value, color=variable)) +
  geom_point() + geom_line() + theme_bw() +
  xlab("M") + ylab("Norm metric") + labs(color="Method") +
  ggtitle("Low noise scenario")





df1 <- readRDS("../data/increasing_M_one_low.RDS")

avg.low <- df1 |> filter(variable == "Low noise matrix") |> select(value)
avg.low <- mean(avg.low$value)

df1[df1$variable == "Low noise matrix",]$value <- avg.low

p1 <- df1 |> ggplot(aes(x=M, y=value, color=variable)) +
  geom_point() + geom_line() + theme_bw() +
  xlab("M") + ylab("Norm metric") + labs(color="Method") +
  ggtitle("High noise scenario")



library(ggpubr)

p.full <- ggarrange(p.indiv,
          ggarrange(p, p1, nrow=1, ncol=2, labels=c("b","c"), common.legend=TRUE,
                    legend="top"),
          nrow=2, labels=c("a",""))

ggsave(p.full, filename="../plots/ambient_sim.png",
       width=8.46, height=7.44, units="in")



