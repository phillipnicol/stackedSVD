n <- 1000
d <- 1000

theta1 <- 2
theta2 <- 0.1

u1 <- rnorm(n); u1 <- u1/sqrt(sum(u1^2))
u2 <- rnorm(n); u2 <- u2/sqrt(sum(u2^2))
v <- rnorm(d); v <- v/sqrt(sum(v^2))

X1 <- theta1 * u1 %*% t(v) + matrix(rnorm(n*d, sd=1/sqrt(d)), n, d)
X2 <- theta2 * u2 %*% t(v) + matrix(rnorm(n*d, sd=1/sqrt(d)), n, d)

my.svd <- irlba::irlba(X1,nv=1)
sigma12 <- my.svd$d[1]^2
c1 <- n/d
c2 <- n/d
if(my.svd$d[1]^2 < (1 + sqrt(c1))^2) {
  print("The signal is too weak to be detected.")
  theta1.hat <- c1
} else {
  theta1.hat <- sqrt(sigma12 - (1 + c1) + sqrt((sigma12 - (1+c1))^2 - 4*c1))/sqrt(2)
}
beta1.hat <- ((theta1.hat^4 - c1)/(theta1.hat^2 * (theta1.hat^2 + 1)))^{1/2}

my.norm <- sum((X2 %*% my.svd$v[,1])^2)
if(my.norm < c2) {
  theta2.hat <- 0
} else{
  theta2.hat <- beta1.hat^{-1}*sqrt(my.norm - c2)
}

print(theta2.hat)
