# Data preparation
```
X <- wheat.X
Y <- wheat.Y

n <- nrow(Y)
p <- ncol(X)

# Genomic relationship matrix
Z <- scale(X)
G <- tcrossprod(Z)/p
I <- diag(n)

# Select only environments 2,4, and 5 to work with
Y <- Y[,c(2,3,4)]

# Matrix to store results
out <- matrix(NA,nrow=ncol(Y),ncol=3)
out[,1] <- colnames(Y)
colnames(out) <- c("Environment","Within","MxE")
```

# MxE model
Model GxE interaction using a marker x environment (MxE) approach that benefits of positively correlated environments 
and can be used in two predictions problems (CV1 and CV2) that mimic 2 evaluation situations:
1. Cross Validation 1 (CV1). Represent a scheme of prediction of lines that have not been evaluated in any field
trials.
2. Cross Validation 2 (CV2). Represent a scheme of prediction of lines that have been evaluated in some but all target environments.

Using a GBLUP approach, the prediction power of the multi-environment MxE model is compared with that that ignores GxE (across-environment) and with the GBLUP model fitted within environment. 

Reference: *[Lopez-Cruz et. al, 2015](https://www.ncbi.nlm.nih.gov/pubmed/25660166)*

## Training-Testing random partitions.
Data is randomly splitted into training and testing using 70% of the data for training and 30%  for testing as depicted in the figures below.

<img src="https://github.com/MarcooLopez/Genomic-Selection-Demo/blob/master/CV1_2_scheme.png" width="450">

CV1 scheme is made in such a way that 30% of the lines are missing in all environments. CV2 scheme is created by having 30% of the entries missing in one environment but present in all the rest of environments.

This procedure of TRN-TST can be repeated many times to allow for estimation of standard errors (SE)

## Cross Validation 1 (CV1)
### TRN and TST sets creation
```
set.seed(123)
env <- c(1,2,3) # choose any set of environments from 1:ncol(Y)
nEnv <- length(env)
Y <- Y[,env]
n <- nrow(Y)
percTST <- 0.3
nTST <- round(percTST*n)
tst <- sample(1:n,size=nTST,replace=FALSE)
YNA <- Y
YNA[tst,]<-NA
```

### Running models
```
#====== Single environments models ===================================
YHatSE <- matrix(nrow=nrow(Y),ncol=ncol(Y),NA)
ETA <- list(G=list(K=G,model='RKHS'))
for(j in 1:nEnv){
    prefix <- paste(colnames(Y)[j],"_",sep="")
    fm <-BGLR(y=YNA[,i],ETA=ETA,nIter=12000,burnIn=2000,saveAt=prefix)
    YHatSE[,j] <- fm$yHat
}

#====== Across environment model (ignoring GxE) ======================
yNA <- as.vector(YNA)

# Fixed effect (env-intercepts)
envID <- rep(env,each=nrow(Y))
ETA <- list(list(~factor(envID)-1,model="FIXED"))

# Main effects of markers
G0 <- kronecker(matrix(nrow=nEnv,ncol=nEnv,1),G)
ETA[[2]] <- list(K=G0,model='RKHS')

# Model Fitting
prefix <- paste(c('Across',colnames(Y),''),collapse='_')
fm <- BGLR(y=yNA,ETA=ETA,nIter=12000,burnIn=2000,saveAt=prefix)
YHatAcross <- matrix(fm$yHat,ncol=nEnv)

#====== MxE Interaction Model =======================================
# Adding interaction terms
    for(j in 1:nEnv){
    tmp <- rep(0,nEnv) ; tmp[j] <- 1; G1 <- kronecker(diag(tmp),G)
    ETA[[(j+2)]] <- list(K=G1,model='RKHS')
}
# Model Fitting
prefix <- paste(c('MxE',colnames(Y),''),collapse='_')
fm <- BGLR(y=yNA,ETA=ETA,nIter=12000,burnIn=2000,saveAt=prefix)
YHatInt <- matrix(fm$yHat,ncol=nEnv)
```
## Results

|       |GBLUP  |B-GBLUP | BRR  | LASSO | Bayes B |
|-------|-------|--------|------|-------|-------|
|Fold 1  | 0.54  | 0.53  | 0.50 | 0.51 | 0.51 |
|Fold 2  | 0.48  | 0.49  | 0.46 | 0.45 | 0.45 |
|Fold 3  | 0.54  | 0.54  | 0.52 | 0.52 | 0.52 |
|Fold 4  | 0.52  | 0.52  | 0.53 | 0.53 | 0.53 |
|Fold 5  | 0.51  | 0.52  | 0.50 | 0.49 | 0.48 |
|Mean    | 0.52  | 0.52  | 0.50 | 0.50 | 0.50 |


#
* **[back](https://github.com/MarcooLopez/Genomic-Selection-Demo/blob/master/README.md)**
