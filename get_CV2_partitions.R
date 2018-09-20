rm(list=ls())
#=========================================================
# User specifications
#=========================================================
# Number of replicates
m <- 100

# Percentage of the data assigned to Testing set
percTST <- 0.3
#=========================================================

# Load data
load("../multiEnvironment/prepData_multi.RData")
n <- nrow(Y)
nEnv <- ncol(Y)

# Creation of seed for repeated randomizations
set.seed(123)
seeds <- round(seq(1E3,1E6,length=m))

nTST <- round(percTST*n)
YNA <- vector("list",m)
nNA <- nEnv*nTST

for(k in 1:m)
{
    set.seed(seeds[k])
    YNA0 <- Y
    
    if(nNA<n){ indexNA <- sample(1:n,nNA,replace=FALSE) }
    if(nNA>=n){
        nRep <- floor(nNA/n)
        remain <- sample(1:n,nNA%%n,replace=FALSE)
        a0 <- sample(1:n,n,replace=FALSE)
        indexNA <- rep(a0,nRep)
        if(length(remain)>0){
            a1 <- floor(length(indexNA)/nTST)*nTST
            a2 <- nNA - a1 - length(remain)
            bb <- sample(a0[!a0%in%remain],a2,replace=FALSE)
            noInIndexNA <- c(rep(a0,nRep-1),a0[!a0%in%bb])
            indexNA <- c(noInIndexNA,bb,remain)
        }
    }
    indexEnv <- rep(1:nEnv,each=nTST)
    for(j in 1:nEnv) YNA0[indexNA[indexEnv==j],j] <- NA
    YNA[[k]] <- YNA0
}

# Save YNA matrix
save(YNA,file="../multiEnvironment/YNA_CV2_multiEnv.RData")