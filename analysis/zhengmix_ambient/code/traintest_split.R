
set.seed(1)
Y.perm <- Y[sample(1:nrow(Y), size=nrow(Y), replace = TRUE),]

Ytrain <- Y.perm[1:3000,]
Ytest <- Y.perm[3001:nrow(Y),]

saveRDS(Ytrain, file="../data/Ytrain.RDS")

Xtrain <- scale(2*sqrt(Ytrain/d), center=TRUE, scale=FALSE)
Xtest <- scale(2*sqrt(Ytest/d), center=TRUE, scale=FALSE)

library(irlba)

svd.train <- irlba(Xtrain, nv=20)
svd.test <- irlba(Xtest, nv=20)

saveRDS(Xtrain, file="../data/Xtrain.RDS")
saveRDS(Xtest, file="../data/Xtest.RDS")
