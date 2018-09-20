setwd("/mnt/home/lopezcru/GS")
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
load("multiEnvironment/prepData_multi.RData")
n <- nrow(Y)

# Creation of seed for repeated randomizations
set.seed(123)
seeds <- round(seq(1E3,1E6,length=m))

nTST <- round(percTST*n)
YNA <- vector("list",m)

for(k in 1:m)
{
    set.seed(seeds[k])
    indexTST <- sample(1:n,size=nTST,replace=FALSE)
    YNA0 <- Y
    YNA0[indexTST,] <- NA
    YNA[[k]] <- YNA0
}

# Save YNA matrix
save(YNA,file="multiEnvironment/YNA_CV1_multiEnv.RData")
