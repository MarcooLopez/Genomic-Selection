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
Model GxE interaction using a marker x environment (MxE) approach that benefits of positively correlated environments. The bennefit is more remarkable depending of the prediction problem faced in breeding programs.

MxE descomposes marker effects into an effect that is common to all environments and an effect that is specific to each environment.

Using a GBLUP approach, the prediction power of the multi-environment MxE model is compared with that that ignores GxE (across-environment) and with the GBLUP model fitted within environment. 

Reference: *[Lopez-Cruz et. al, 2015](https://www.ncbi.nlm.nih.gov/pubmed/25660166)*

## Training-Testing random partitions.
The prediction power of the model will be assessed using the training-testing (TRN-TST) random partitions approach. 
Data is randomly splitted into training and testing sets. Model parameters are estimated in training set and model is tested in TST set.  Two main estimations problems are addressed using the MxE interaction model. 

1. Cross Validation 1 (CV1). Represent a scheme of prediction of lines that have not been evaluated in any field
trials.
2. Cross Validation 2 (CV2). Represent a scheme of prediction of lines that have been evaluated in some but all target environments.

<img src="https://github.com/MarcooLopez/Genomic-Selection-Demo/blob/master/CV1_2_scheme.png" width="450">

In our case, we will use 70% of the data for training set and the remaining 30% for the testing set.
For CV1, we will create a scheme in which 30% of the lines are missing in all environments. 
CV2 scheme is created by having 30% of the entries missing in one environment but present in all the rest of environments.

This procedure of TRN-TST can be repeated many times to allow for estimation of standard errors (SE).

### Cross Validation 1 (CV1)

Code below will generate a matrix YNA containing "NA" values for the entries corresponding to the TST set mimicing the CV1 prediction problem. 

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

### Cross Validation 2 (CV2)

Code below will generate a matrix YNA containing "NA" values for the entries corresponding to the TST set mimicing the CV2 prediction problem. 

```
set.seed(12345)

env <- c(4,5) # choose any set of environments from 1:ncol(Y)
nEnv <- length(env)
Y <- Y[,env]
n <- nrow(Y)

percTST<-0.3
nTST <- round(percTST*n)
nNA <- nEnv*nTST
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
YNA <- Y
for(j in 1:nEnv) YNA[indexNA[indexEnv==j],j] <- NA

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
### Results
Computing the within-environment correlation
```
COR <- matrix(nrow=length(env),ncol=3,NA)
colnames(COR) <- c('SingleEnv', 'AcrossEnv', 'MxE')
rownames(COR) <- colnames(Y)
for(j in 1:nEnv){
    tst <- which(is.na(YNA[,j]))
    COR[j,1] <- cor(Y[tst,j],YHatSE[tst,j])
    COR[j,2] <- cor(Y[tst,j],YHatAcross[tst,j])
    COR[j,3] <- cor(Y[tst,j],YHatInt[tst,j])
}
COR
```

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
