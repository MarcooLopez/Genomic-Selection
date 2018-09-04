
# Multi-environment models

### * Reaction Norm model
Is an extention of the G-BLUP model that incorporates GxE by introducing covariance structures as a funcion of the marker information.
The reaction norm model can model main and interaction effects of environmental covariates (EC) and markers.


### * MxE model
It models GxE interaction using a marker-by-environment (MxE) approach that benefits of positively correlated environments. MxE descomposes marker effects into an effect that is common to all environments and an effect that is specific to each environment.

Performance of the *Reaction Norm* and *MxE model* will be compared with that of the *across-environments* model that ignores GxE modeling and the *single-environment* model which is fitted within each environment.

## Model assessment
### Training-Testing random partitions.
The prediction power of the model will be assessed using the training-testing (TRN-TST) random partitions approach. 
Data is randomly splitted into training and testing sets. Model parameters are estimated in training set and model is tested in TST set.  Two main estimations problems are addressed using the multi-environments models. 

**1. Cross Validation 1 (CV1)**. Represent a scheme of prediction of lines that have not been evaluated in any field
trials.

**2. Cross Validation 2 (CV2)**. Represent a scheme of prediction of lines that have been evaluated in some but all target environments. Thus, prediction of non-evaluated lines benefits from borrowing of information from lines that were evaluated in other environments.

<img src="https://github.com/MarcooLopez/Genomic-Selection-Demo/blob/master/CV1_2_scheme.png" width="450">

In our case, we will use 70% of the data for training set and the remaining 30% for the testing set.
For CV1, we will create a scheme in which 30% of the lines are missing in all environments. 
CV2 scheme is created by having 30% of the entries missing in one environment but present in all the rest of environments.

This procedure of TRN-TST can be repeated many times to allow for estimation of standard errors (SE).

## Data preparation
### Load data, generate G-matrix
```
# Load libraries
library(BGLR)
library(rrBLUP)

# Load data
data(wheat)
X <- wheat.X
Y <- wheat.Y

n <- nrow(Y)
p <- ncol(X)

# Select environments. For instance, environments 2,4, and 5
Y <- Y[,c(2,3,4)]

y <- as.vector(Y)

# Genomic relationship matrix
M <- scale(X)
G <- tcrossprod(M)/p

# Design matrix for individuals. It connects individuals with environments
GID <- factor(rep(rownames(Y),ncol(Y)),levels=rownames(Y))
Zg <- model.matrix(~GID-1)

# Design matrix for environments. Used in the multi-environment R-Norm model
envID <- factor(rep(colnames(Y),each=nrow(Y)),levels=colnames(Y))
ZE <- model.matrix(~envID-1)   

#  Covariance structure for genetic and environmental effects
ZgGZgt <- Zg%*%G%*%t(Zg)  
ZEZEt <- tcrossprod(ZE)
```

## Running models

### 1. Variance components estimation 
Code below can be used after 'data preparation' part to fit all the models and to extract variance components.

```
set.seed(123)
nEnv <- ncol(Y)

# Eigen decomposition (to speed computational time) of main effects of markers, ZgGZgt'
eigen_G0 <- eigen(ZgGZgt)

# Matrix to store results. It will save variance components for each model
outVAR <- matrix(NA,ncol=4,nrow=1+2*ncol(Y))
dimnames(outVAR) <- list(c("Main",rep(paste0("Env ",colnames(Y)),2)),c("Single","Across","MxE","R-Norm"))

# Number of iterations and burn-in for Bayesian models
nIter <- 30000; burnIn <- 2000

#==============================================================
# 1. Single environment (within-environment) model, ignoring GxE effect
#==============================================================
ETA <- list(G=list(K=G,model='RKHS'))
for(env in 1:nEnv){
    fm <-BGLR(y=Y[,env],ETA=ETA,nIter=nIter,burnIn=burnIn)
    outVAR[env+1,1] <- fm$ETA[[1]]$varU
    outVAR[env+4,1] <- fm$varE
}

#==============================================================
# 2. Across-environments model. Factor 'environment' as fixed effect
#==============================================================
ETA <- list(list(~envID-1,model="FIXED"))
ETA[[2]] <- list(V=eigen_G0$vectors,d=eigen_G0$values,model='RKHS')

# Model Fitting
fm <- BGLR(y=y,ETA=ETA,nIter=nIter,burnIn=burnIn)
outVAR[1,2] <- fm$ETA[[2]]$varU
outVAR[(1:nEnv)+4,2] <- fm$varE

#==============================================================
# 3. MxE interaction model. Factor 'environment' as fixed effect
#==============================================================
ETA <- list(list(~envID-1,model="FIXED"))
ETA[[2]] <- list(V=eigen_G0$vectors,d=eigen_G0$values,model='RKHS')

# Adding interaction terms
for(env in 1:nEnv){
    tmp <- rep(0,nEnv) ; tmp[env] <- 1; G1 <- kronecker(diag(tmp),G)
    eigen_G1 <- eigen(G1)
    ETA[[(env+2)]] <- list(V=eigen_G1$vectors,d=eigen_G1$values,model='RKHS')
}

# Model Fitting
fm <- BGLR(y=y,ETA=ETA,nIter=nIter,burnIn=burnIn)
outVAR[1,3] <- fm$ETA[[2]]$varU
for(env in 1:nEnv) outVAR[env+1,3] <- fm$ETA[[env+2]]$varU
outVAR[(1:nEnv)+4,3] <- fm$varE

#==============================================================
# 4. Reaction-Norm model. Factor 'environment' as fixed effect
#==============================================================
ETA <- list(list(~envID-1,model="FIXED"))
ETA[[2]] <- list(V=eigen_G0$vectors,d=eigen_G0$values,model='RKHS')

# Adding GxE interaction term as Hadamart product
GE <- ZgGZgt*ZEZEt
eigen_GE <- eigen(GE)
ETA[[3]] <- list(V=eigen_GE$vectors,d=eigen_GE$values,model="RKHS")

# Model Fitting
fm <- BGLR(y=y,ETA=ETA,nIter=nIter,burnIn=burnIn)
outVAR[1,4] <- fm$ETA[[2]]$varU
outVAR[(1:nEnv)+1,4] <- fm$ETA[[3]]$varU
outVAR[(1:nEnv)+4,4] <- fm$varE

```

#### Results
<img src="https://github.com/MarcooLopez/Genomic-Selection/blob/master/varComp.png" width="360">

##
### 2. Replicates of partitions to obtain standard deviations of predictions
Using a GBLUP approach, the prediction power of the multi-environment models (MxE and Reaction Norm) will be compared with that that ignores GxE (across-environment) and with the GBLUP model fitted within environment. 

#### 2.1. Training-Testing partitions
After running the 'data preparation' part, it can be chosen either to perform CV1 or CV2 aproaches

* **Cross Validation 1 (CV1)**

Code below will generate a matrix YNA containing "NA" values for the entries corresponding to the TST set mimicing the CV1 prediction problem. It generates a 'list' with 'm' matrices containing the TRN-TST partitions

```
# Number of replicates
m <- 100

# Percentage of the data assigned to Testing set
percTST <- 0.3

# Creation of seed for repeated randomizations
set.seed(123)
seeds <- round(seq(1E3,1E6,length=m))

nEnv <- ncol(Y)
n <- nrow(Y)
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
```

* **Cross Validation 2 (CV2)**

Code below will generate a matrix YNA containing "NA" values for the entries corresponding to the TST set mimicing the CV2 prediction problem. It generates a 'list' with 'm' matrices containing the TRN-TST partitions

```
# Number of replicates
m <- 100

# Percentage of the data assigned to Testing set
percTST <- 0.3

# Creation of seed for repeated randomizations
set.seed(123)
seeds <- round(seq(1E3,1E6,length=m))

nEnv <- ncol(Y)
n <- nrow(Y)
nTST <- round(percTST*n)

nNA <- nEnv*nTST
YNA <- vector("list",m)
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
```

After running the code to generate partitions for either CV1 or CV2 scenarios, the following code can be run to fit the models repeatealy for all partitions.

```
# Models
models <- c("single","across","MxE","R-Norm")

# Choose one model. 1: single; 2:across; 3:MxE; 4:R-Norm
mod <- 2

# Matrix to store results. It will save the accuracy for each partition
outCOR <- matrix(NA,nrow=m,ncol=nEnv)
colnames(outCOR) <- colnames(Y)

# Eigen decomposition (to speed computational time) of main effects of markers, ZgGtZg'
eigen_G0 <- eigen(ZgGZgt)

# Number of iterations and burn-in for Bayesian models
nIter <- 1200
burnIn <- 200

model <- models[mod]
for(k in 1:10)   # Loop for the replicates
{
    YNA0 <- YNA[[k]]
    yNA <- as.vector(YNA0)
    
    #==============================================================
    # 1. Single environment (within-environment) model, ignoring GxE effect
    #==============================================================
    if(model=="single")
    {
        ETA <- list(G=list(K=G,model='RKHS'))
        for(env in 1:nEnv){
            fm <-BGLR(y=YNA0[,env],ETA=ETA,nIter=nIter,burnIn=burnIn)
            indexTST <- which(is.na(YNA0[,env]))
            outCOR[k,env] <- cor(Y[indexTST,env],fm$yHat[indexTST])
        }
    }

    #==============================================================
    # 2. Across-environments model. Factor 'environment' as fixed effect
    #==============================================================
    if(model=="across")
    {
        ETA <- list(list(~envID-1,model="FIXED"))
        ETA[[2]] <- list(V=eigen_G0$vectors,d=eigen_G0$values,model='RKHS')

        # Model Fitting
        fm <- BGLR(y=yNA,ETA=ETA,nIter=nIter,burnIn=burnIn)
        YHat <- matrix(fm$yHat,ncol=nEnv)
        for(env in 1:nEnv){
            indexTST <- which(is.na(YNA0[,env]))
            outCOR[k,env] <- cor(Y[indexTST,env],YHat[indexTST,env])
        }
    }
    
    #==============================================================
    # 3. MxE interaction model. Factor 'environment' as fixed effect
    #==============================================================
    if(model=="MxE")
    {
        ETA <- list(list(~envID-1,model="FIXED"))
        ETA[[2]] <- list(V=eigen_G0$vectors,d=eigen_G0$values,model='RKHS')

        # Adding interaction terms
        for(env in 1:nEnv){
            tmp <- rep(0,nEnv) ; tmp[env] <- 1; G1 <- kronecker(diag(tmp),G)
            eigen_G1 <- eigen(G1)
            ETA[[(env+2)]] <- list(V=eigen_G1$vectors,d=eigen_G1$values,model='RKHS')
        }

        # Model Fitting
        fm <- BGLR(y=yNA,ETA=ETA,nIter=nIter,burnIn=burnIn)
        YHat <- matrix(fm$yHat,ncol=nEnv)
        for(env in 1:nEnv){
            indexTST <- which(is.na(YNA0[,env]))
            outCOR[k,env] <- cor(Y[indexTST,env],YHat[indexTST,env])
        }
    }

    #==============================================================
    # 4. Reaction-Norm model. Factor 'environment' as fixed effect
    #==============================================================
    if(model=="R-Norm")
        ETA <- list(list(~envID-1,model="FIXED"))
        ETA[[2]] <- list(V=eigen_G0$vectors,d=eigen_G0$values,model='RKHS')

        # Adding GxE interaction term as Hadamart product
        GE <- ZgGZgt*ZEZEt
        eigen_GE <- eigen(GE)
        ETA[[3]] <- list(V=eigen_GE$vectors,d=eigen_GE$values,model="RKHS")

        # Model Fitting
        fm <- BGLR(y=yNA,ETA=ETA,nIter=nIter,burnIn=burnIn)
        YHat <- matrix(fm$yHat,ncol=nEnv)
        for(env in 1:nEnv){
            indexTST <- which(is.na(YNA0[,env]))
            outCOR[k,env] <- cor(Y[indexTST,env],YHat[indexTST,env])
        }
    }
}

```

### Results
Computing the within-environment correlation for all three models
```
COR <- matrix(nrow=length(env),ncol=4,NA)
colnames(COR) <- c('SingleEnv', 'AcrossEnv', 'MxE','RNorm')
rownames(COR) <- colnames(Y)
for(j in 1:nEnv){
    tst <- which(is.na(YNA[,j]))
    COR[j,1] <- cor(Y[tst,j],YHat1[tst,j])
    COR[j,2] <- cor(Y[tst,j],YHat2[tst,j])
    COR[j,3] <- cor(Y[tst,j],YHat3[tst,j])
    COR[j,4] <- cor(Y[tst,j],YHat4[tst,j])
}
COR
```

**CV1. One TRN-TST partition**

|       |Single-Env |Across-Env | MxE  | RNorm |
|-------|-------|--------|------|------|
|Env 2  | 0.41  | 0.33  | 0.38 | 0.38 |
|Env 4  | 0.29  | 0.32  | 0.32 | 0.32 |
|Env 5  | 0.47  | 0.49  | 0.51 | 0.51 |

##
**CV2. One TRN-TST partition**

|       |Single-Env |Across-Env | MxE  | RNorm |
|-------|-------|--------|------|-----|
|Env 2  | 0.41  | 0.62  | 0.64 | 0.63 |
|Env 4  | 0.38  | 0.61  | 0.59 | 0.59 |
|Env 5  | 0.43  | 0.41  | 0.46 | 0.45 |


#
* **[back](https://github.com/MarcooLopez/Genomic-Selection-Demo/blob/master/README.md)**

#
# References
* Jarquín, D., Crossa, J., Lacaze, X., Du Cheyron, P., Daucourt, J., Lorgeou, J., … de los Campos, G. (2014). **A reaction norm model for genomic selection using high-dimensional genomic and environmental data**. Theoretical and Applied Genetics, 127(3), 595–607. 
* Lopez-Cruz, M., Crossa, J., Bonnett, D., Dreisigacker, S., Poland, J., Jannink, J.-L., … de los Campos, G. (2015). **Increased Prediction Accuracy in Wheat Breeding Trials Using a Marker × Environment Interaction Genomic Selection Model**. G3: Genes, Genomes, Genetics, 5(4), 569–582. 



