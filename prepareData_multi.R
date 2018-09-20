rm(list=ls())
library(BGLR)

# Load data
data(wheat)
X <- wheat.X
Y <- wheat.Y

# Select environments. For instance, environments 2,4, and 5
Y <- Y[,c(2,3,4)]

# Genomic relationship matrix
M <- scale(X)
G <- tcrossprod(M)/ncol(X)

# Design matrix for individuals. It connects individuals with environments
GID <- factor(rep(rownames(Y),ncol(Y)),levels=rownames(Y))
Zg <- model.matrix(~GID-1)

# Design matrix for environments. Used in the multi-environment R-Norm model
envID <- factor(rep(colnames(Y),each=nrow(Y)),levels=colnames(Y))
ZE <- model.matrix(~envID-1)   

#  Covariance structure for effects
ZgGZgt <- Zg%*%G%*%t(Zg)    # Genetic effect  
ZEZEt <- tcrossprod(ZE)     # Environmental effect
GE <- ZgGZgt*ZEZEt          # GxE interaction term (R-Norm model)

# Eigen decomposition (to speed computational time)
eigen_G <- eigen(G)
eigen_G0 <- eigen(ZgGZgt)
eigen_GE <- eigen(GE)

# Interaction terms (MxE model)
MxE_eigen <- vector("list",ncol(Y))
for(env in 1:ncol(Y)){ 
    tmp <- rep(0,ncol(Y)) ; tmp[env] <- 1; G1 <- kronecker(diag(tmp),G)
    MxE_eigen[[env]] <- eigen(G1)
}

# Save prepared data
dir.create("../multiEnvironment")
save(Y,envID,eigen_G,eigen_G0,eigen_GE,MxE_eigen,file="../multiEnvironment/prepData_multi.RData")