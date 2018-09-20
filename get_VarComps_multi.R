rm(list=ls())
library(BGLR)
load("../multiEnvironment/prepData_multi.RData")
n <- nrow(Y)
nEnv <- ncol(Y)
y <- as.vector(Y)

set.seed(123)

# Matrix to store results. It will save variance components for each model
outVAR <- matrix(NA,ncol=4,nrow=1+2*nEnv)
dimnames(outVAR) <- list(c("Main",rep(paste0("Env ",colnames(Y)),2)),c("Single","Across","MxE","R-Norm"))

# Number of iterations and burn-in for Bayesian models
nIter <- 30000; burnIn <- 2000

#--------------------------------------------------------
# 1. Single environment (within-environment) model
#--------------------------------------------------------
ETA <- list(G=list(V=eigen_G$vectors,d=eigen_G$values,model='RKHS'))
for(env in 1:nEnv){
    fm <-BGLR(y=Y[,env],ETA=ETA,nIter=nIter,burnIn=burnIn)
    outVAR[env+1,1] <- fm$ETA[[1]]$varU
    outVAR[env+4,1] <- fm$varE
}

#--------------------------------------------------------
# 2. Across-environments model
#--------------------------------------------------------
ETA <- list(list(~envID-1,model="FIXED"))
ETA[[2]] <- list(V=eigen_G0$vectors,d=eigen_G0$values,model='RKHS')

# Model Fitting
fm <- BGLR(y=y,ETA=ETA,nIter=nIter,burnIn=burnIn)
outVAR[1,2] <- fm$ETA[[2]]$varU
outVAR[(1:nEnv)+4,2] <- fm$varE

#--------------------------------------------------------
# 3. MxE interaction model
#--------------------------------------------------------
ETA <- list(list(~envID-1,model="FIXED"))
ETA[[2]] <- list(V=eigen_G0$vectors,d=eigen_G0$values,model='RKHS')

# Adding interaction terms
for(env in 1:nEnv){
    eigen_G1 <- MxE_eigen[[env]]
    ETA[[(env+2)]] <- list(V=eigen_G1$vectors,d=eigen_G1$values,model='RKHS')
}

# Model Fitting
fm <- BGLR(y=y,ETA=ETA,nIter=nIter,burnIn=burnIn)
outVAR[1,3] <- fm$ETA[[2]]$varU
for(env in 1:nEnv) outVAR[env+1,3] <- fm$ETA[[env+2]]$varU
outVAR[(1:nEnv)+4,3] <- fm$varE

#--------------------------------------------------------
# 4. Reaction-Norm model
#--------------------------------------------------------
ETA <- list(list(~envID-1,model="FIXED"))
ETA[[2]] <- list(V=eigen_G0$vectors,d=eigen_G0$values,model='RKHS')
ETA[[3]] <- list(V=eigen_GE$vectors,d=eigen_GE$values,model="RKHS")

# Model Fitting
fm <- BGLR(y=y,ETA=ETA,nIter=nIter,burnIn=burnIn)
outVAR[1,4] <- fm$ETA[[2]]$varU
outVAR[(1:nEnv)+1,4] <- fm$ETA[[3]]$varU
outVAR[(1:nEnv)+4,4] <- fm$varE
outVAR

# Save results
write.table(outVAR,file="../multiEnvironment/varComps.csv",sep=",",row.names=F)