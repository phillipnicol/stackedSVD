n <- 5000
d <- 5000

theta1 <- 1.1
theta2 <- 0.5

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

theta2.hat - theta2

###Singular vector recovery
#Unweighted Stack SVD
X <- rbind(X1, X2)
my.svd <- irlba::irlba(X,nv=1)
v.hat <- my.svd$v 
sum(v.hat*v)^2

#Estimated weight Stack SVD
w.hat <- sqrt(c(theta1.hat^2/(theta1.hat^2 + 1), theta2.hat^2/(theta2.hat^2 + 1)))
X <- rbind(w.hat[1] * X1, w.hat[2] * X2)
my.svd <- irlba::irlba(X,nv=1)
v.hat <- my.svd$v 
sum(v.hat*v)^2

#Optimal weight Stack SVD
w.opt <- sqrt(c(theta1^2/(theta1^2 + 1), theta2^2/(theta2^2 + 1)))
X <- rbind(w.opt[1] * X1, w.opt[2] * X2)
my.svd <- irlba::irlba(X,nv=1)
v.hat <- my.svd$v 
sum(v.hat*v)^2
