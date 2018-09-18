
# Multi-environment models

* **Single-environment model**

This model is obtained by regressing the phenotype vector containing the *n* records available in the *k*<sup>th</sup> environment (*k=1,2,...,s* environments), ![](https://latex.codecogs.com/gif.latex?%5Ctextbf%7By%7D_k%3D%5By_%7B1k%7D%2C...%2Cy_%7Bnk%7D%5D%27), where *i=1,2,...,n* indexes lines (individuals), on *p* markers using a linear model in the form

<p align="center">
<img src="https://latex.codecogs.com/gif.latex?%5Ctextbf%7By%7D_k%3D%5Cmu_k%5Ctextbf%7B1%7D&plus;%5Ctextbf%7BX%7D_k%5Cboldsymbol%7B%5Cbeta%7D_k&plus;%5Cboldsymbol%7B%5Cvarepsilon%7D_k">
</p>

where ![](https://latex.codecogs.com/gif.latex?%5Cmu_k) is the intercept, ![](https://latex.codecogs.com/gif.latex?%5Ctextbf%7B1%7D%3D%5B1%2C...%2C1%5D%27) is a *n*-vector of ones, ![](https://latex.codecogs.com/gif.latex?%5Ctextbf%7BX%7D_k%3D%5C%7Bx_%7Bijk%7D%5C%7D) is the ![](https://latex.codecogs.com/gif.latex?n%5Ctimes%20p) matrix of centered and standardized markers available in the *k*<sup>th</sup> environment, ![](https://latex.codecogs.com/gif.latex?%5Cboldsymbol%7B%5Cbeta%7D_k%3D%5B%5Cbeta_%7B1k%7D%2C...%2C%5Cbeta_%7Bpk%7D%5D%27) is a *p*-vector of marker effects, and ![](https://latex.codecogs.com/gif.latex?%5Cboldsymbol%7B%5Cvarepsilon%7D_k%3D%5B%5Cvarepsilon_%7Bik%7D%2C...%2C%5Cvarepsilon_%7Bnk%7D%5D%27) is the *n*-vector of residuals.

This model can be represented as a G-BLUP model by setting ![](https://latex.codecogs.com/gif.latex?%5Ctextbf%7Bu%7D_k%3D%5Ctextbf%7BX%7D_k%5Cboldsymbol%7B%5Cbeta%7D_k), this is: 

<p align="center">
<img src="https://latex.codecogs.com/gif.latex?%5Ctextbf%7By%7D_k%3D%5Cmu_k%5Ctextbf%7B1%7D&plus;%5Ctextbf%7Bu%7D_k&plus;%5Cboldsymbol%7B%5Cvarepsilon%7D_k">
</p>

where ![](https://latex.codecogs.com/gif.latex?%5Ctextbf%7Bu%7D_k%3D%5Bu_%7B1k%7D%2C...%2Cu_%7Bnk%7D%5D%27) is a random effect assumed ![](https://latex.codecogs.com/gif.latex?%5Ctextbf%7Bu%7D_k%5Csim%20N%28%5Ctextbf%7B0%7D%2C%5Csigma%5E2_%7Bu_k%7D%5Ctextbf%7BG%7D_k%29) with the genomic relationship matrix estimated as ![](https://latex.codecogs.com/gif.latex?%5Ctextbf%7BG%7D_k%3D%5Ctextbf%7BX%7D_k%5Ctextbf%7BX%7D_k%27/p).

* **Across-environments model**

This model assumes that effects of markers are the same across environments, this is ![](https://latex.codecogs.com/gif.latex?%5Cboldsymbol%7B%5Cbeta%7D_1%3D%5Cboldsymbol%7B%5Cbeta%7D_2%3D...%3D%5Cboldsymbol%7B%5Cbeta%7D_s%3D%5Cboldsymbol%7B%5Cbeta%7D). The model above can be simultaneously fitted for all environments as (assume *s=3* environments)

<p align="center">
<img src="https://latex.codecogs.com/gif.latex?%5Cbegin%7Bbmatrix%7D%20%5Ctextbf%7By%7D_1%5C%5C%20%5Ctextbf%7By%7D_2%5C%5C%20%5Ctextbf%7By%7D_3%20%5Cend%7Bbmatrix%7D%3D%5Cbegin%7Bbmatrix%7D%20%5Cmu_1%5Ctextbf%7B1%7D%5C%5C%20%5Cmu_2%5Ctextbf%7B1%7D%5C%5C%20%5Cmu_3%5Ctextbf%7B1%7D%20%5Cend%7Bbmatrix%7D&plus;%5Cbegin%7Bbmatrix%7D%20%5Ctextbf%7BX%7D_1%5C%5C%20%5Ctextbf%7BX%7D_2%5C%5C%20%5Ctextbf%7BX%7D_3%20%5Cend%7Bbmatrix%7D%5Cboldsymbol%7B%5Cbeta%7D&plus;%20%5Cbegin%7Bbmatrix%7D%20%5Cboldsymbol%7B%5Cvarepsilon%7D_1%5C%5C%20%5Cboldsymbol%7B%5Cvarepsilon%7D_2%5C%5C%20%5Cboldsymbol%7B%5Cvarepsilon%7D_3%20%5Cend%7Bbmatrix%7D">
</p>

Similarly, it can be represented also as a G-BLUP model, by making ![](https://latex.codecogs.com/gif.latex?%5Ctextbf%7Bu%7D_k%3D%5Ctextbf%7BX%7D_k%5Cboldsymbol%7B%5Cbeta%7D), as: 

<p align="center">
<img src="https://latex.codecogs.com/gif.latex?%5Cbegin%7Bbmatrix%7D%20%5Ctextbf%7By%7D_1%5C%5C%20%5Ctextbf%7By%7D_2%5C%5C%20%5Ctextbf%7By%7D_3%20%5Cend%7Bbmatrix%7D%3D%5Cbegin%7Bbmatrix%7D%20%5Cmu_1%5Ctextbf%7B1%7D%5C%5C%20%5Cmu_2%5Ctextbf%7B1%7D%5C%5C%20%5Cmu_3%5Ctextbf%7B1%7D%20%5Cend%7Bbmatrix%7D&plus;%5Cbegin%7Bbmatrix%7D%20%5Ctextbf%7Bu%7D_1%5C%5C%20%5Ctextbf%7Bu%7D_2%5C%5C%20%5Ctextbf%7Bu%7D_3%20%5Cend%7Bbmatrix%7D&plus;%20%5Cbegin%7Bbmatrix%7D%20%5Cboldsymbol%7B%5Cvarepsilon%7D_1%5C%5C%20%5Cboldsymbol%7B%5Cvarepsilon%7D_2%5C%5C%20%5Cboldsymbol%7B%5Cvarepsilon%7D_3%20%5Cend%7Bbmatrix%7D">,
</p>

where the random effect ![](https://latex.codecogs.com/gif.latex?%5Ctextbf%7Bu%7D_0%3D%5B%5Ctextbf%7Bu%7D%27_1%2C%5Ctextbf%7Bu%7D%27_2%2C%5Ctextbf%7Bu%7D%27_3%5D%27) is assumed ![](https://latex.codecogs.com/gif.latex?%5Ctextbf%7Bu%7D_0%5Csim%20N%28%5Ctextbf%7B0%7D%2C%5Csigma%5E2_%7Bu_0%7D%5Ctextbf%7BG%7D_0%29) with ![](https://latex.codecogs.com/gif.latex?%5Ctextbf%7BG%7D_0) being the marker-derived genomic relationship calculated as

<p align="center">
<img src="https://latex.codecogs.com/gif.latex?%5Ctextbf%7BG%7D_0%3D%5Cbegin%7Bbmatrix%7D%20%5Ctextbf%7BX%7D_1%5Ctextbf%7BX%7D%27_1%5Cquad%20%5Ctextbf%7BX%7D_1%5Ctextbf%7BX%7D%27_2%5Cquad%20%5Ctextbf%7BX%7D_1%5Ctextbf%7BX%7D%27_3%5C%5C%20%5Ctextbf%7BX%7D_2%5Ctextbf%7BX%7D%27_1%5Cquad%20%5Ctextbf%7BX%7D_2%5Ctextbf%7BX%7D%27_2%5Cquad%20%5Ctextbf%7BX%7D_2%5Ctextbf%7BX%7D%27_3%5C%5C%20%5Ctextbf%7BX%7D_3%5Ctextbf%7BX%7D%27_1%5Cquad%20%5Ctextbf%7BX%7D_3%5Ctextbf%7BX%7D%27_2%5Cquad%20%5Ctextbf%7BX%7D_3%5Ctextbf%7BX%7D%27_3%5C%5C%20%5Cend%7Bbmatrix%7D/p">
</p>

* **MxE model**

It models GxE interaction using a marker-by-environment (MxE) approach that benefits of positively correlated environments. MxE model descomposes the effect of the *j*<sup>th</sup> marker on the *k*<sup>th</sup> environment, ![](https://latex.codecogs.com/gif.latex?%5Cbeta_%7Bjk%7D), into an effect that is common to all environments (*b*<sub>j0</sub>) and an effect that is specific to each environment (*b*<sub>jk</sub>), this is ![](https://latex.codecogs.com/gif.latex?%5Cbeta_%7Bjk%7D%3Db_%7Bj0%7D&plus;b_%7Bjk%7D). Thus, the multi-environmental model is

<p align="center">
<img src="https://latex.codecogs.com/gif.latex?%5Cbegin%7Bbmatrix%7D%20%5Ctextbf%7By%7D_1%5C%5C%20%5Ctextbf%7By%7D_2%5C%5C%20%5Ctextbf%7By%7D_3%20%5Cend%7Bbmatrix%7D%3D%5Cbegin%7Bbmatrix%7D%20%5Cmu_1%5Ctextbf%7B1%7D%5C%5C%20%5Cmu_2%5Ctextbf%7B1%7D%5C%5C%20%5Cmu_3%5Ctextbf%7B1%7D%20%5Cend%7Bbmatrix%7D&plus;%5Cbegin%7Bbmatrix%7D%20%5Ctextbf%7BX%7D_1%5C%5C%20%5Ctextbf%7BX%7D_2%5C%5C%20%5Ctextbf%7BX%7D_3%20%5Cend%7Bbmatrix%7D%5Ctextbf%7Bb%7D_0&plus;%20%5Cbegin%7Bbmatrix%7D%20%5Ctextbf%7BX%7D_1%5Cquad%5Ctextbf%7B0%7D%5Cquad%5Ctextbf%7B0%7D%5C%5C%20%5Ctextbf%7B0%7D%5Cquad%5Ctextbf%7BX%7D_2%5Cquad%5Ctextbf%7B0%7D%5C%5C%20%5Ctextbf%7B0%7D%5Cquad%5Ctextbf%7B0%7D%5Cquad%5Ctextbf%7BX%7D_3%20%5Cend%7Bbmatrix%7D%20%5Cbegin%7Bbmatrix%7D%20%5Ctextbf%7Bb%7D_1%5C%5C%20%5Ctextbf%7Bb%7D_2%5C%5C%20%5Ctextbf%7Bb%7D_3%20%5Cend%7Bbmatrix%7D%20&plus;%5Cbegin%7Bbmatrix%7D%20%5Cboldsymbol%7B%5Cvarepsilon%7D_1%5C%5C%20%5Cboldsymbol%7B%5Cvarepsilon%7D_2%5C%5C%20%5Cboldsymbol%7B%5Cvarepsilon%7D_3%20%5Cend%7Bbmatrix%7D">
</p>

Likewise, by making ![](https://latex.codecogs.com/gif.latex?%5Ctextbf%7Bu%7D_%7B0,k%7D%3D%5Ctextbf%7BX%7D_k%5Ctextbf%7Bb%7D_0) and ![](https://latex.codecogs.com/gif.latex?%5Ctextbf%7Bu%7D_%7B1,k%7D%3D%5Ctextbf%7BX%7D_k%5Ctextbf%7Bb%7D_k), the G-BLUP representation is

<p align="center">
<img src="https://latex.codecogs.com/gif.latex?%5Cbegin%7Bbmatrix%7D%20%5Ctextbf%7By%7D_1%5C%5C%20%5Ctextbf%7By%7D_2%5C%5C%20%5Ctextbf%7By%7D_3%20%5Cend%7Bbmatrix%7D%3D%5Cbegin%7Bbmatrix%7D%20%5Cmu_1%5Ctextbf%7B1%7D%5C%5C%20%5Cmu_2%5Ctextbf%7B1%7D%5C%5C%20%5Cmu_3%5Ctextbf%7B1%7D%20%5Cend%7Bbmatrix%7D&plus;%20%5Cbegin%7Bbmatrix%7D%20%5Ctextbf%7Bu%7D_%7B0%2C1%7D%5C%5C%20%5Ctextbf%7Bu%7D_%7B0%2C2%7D%5C%5C%20%5Ctextbf%7Bu%7D_%7B0%2C3%7D%20%5Cend%7Bbmatrix%7D&plus;%20%5Cbegin%7Bbmatrix%7D%20%5Ctextbf%7Bu%7D_%7B1%2C1%7D%5C%5C%20%5Ctextbf%7Bu%7D_%7B1%2C2%7D%5C%5C%20%5Ctextbf%7Bu%7D_%7B1%2C3%7D%20%5Cend%7Bbmatrix%7D&plus;%20%5Cbegin%7Bbmatrix%7D%20%5Cboldsymbol%7B%5Cvarepsilon%7D_1%5C%5C%20%5Cboldsymbol%7B%5Cvarepsilon%7D_2%5C%5C%20%5Cboldsymbol%7B%5Cvarepsilon%7D_3%20%5Cend%7Bbmatrix%7D">
</p>

* **Reaction Norm model**
Is an extention of the G-BLUP model that incorporates GxE by introducing covariance structures as a funcion of the marker information.
The reaction norm model can model main and interaction effects of environmental covariates (EC) and markers.


Performance of the *Reaction Norm* and *MxE model* will be compared with that of the *across-environments* model that ignores GxE modeling and the *single-environment* model which is fitted within each environment.

## Model assessment
### Training-Testing random partitions.
The prediction power of the model will be assessed using the training-testing (TRN-TST) random partitions approach. 
Data is randomly splitted into training and testing sets. Model parameters are estimated in training set and model is tested in TST set.  Two main estimations problems are addressed using the multi-environments models. 

*   **Cross Validation 1 (CV1)**. Represent a scheme of prediction of lines that have not been evaluated in any field
trials.

*   **Cross Validation 2 (CV2)**. Represent a scheme of prediction of lines that have been evaluated in some but all target environments. Thus, prediction of non-evaluated lines benefits from borrowing of information from lines that were evaluated in other environments.

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

# Select environments. For instance, environments 2,4, and 5
Y <- Y[,c(2,3,4)]

n <- nrow(Y)
p <- ncol(X)
nEnv <- ncol(Y)
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

#  Covariance structure for effects
ZgGZgt <- Zg%*%G%*%t(Zg)    # Genetic effect  
ZEZEt <- tcrossprod(ZE)     # Environmental effect
GE <- ZgGZgt*ZEZEt          # GxE interaction term (R-Norm model)

# Eigen decomposition (to speed computational time)
eigen_G <- eigen(G)
eigen_G0 <- eigen(ZgGZgt)
eigen_GE <- eigen(GE)

# Interaction terms (MxE model)
MxE_eigen <- vector("list",nEnv)
for(env in 1:nEnv){ 
    tmp <- rep(0,nEnv) ; tmp[env] <- 1; G1 <- kronecker(diag(tmp),G)
    MxE_eigen[[env]] <- eigen(G1)
}
```

## Running models

### 1. Variance components estimation 
Code below can be used after 'data preparation' part to fit all the models and to extract variance components.

```
set.seed(123)

# Matrix to store results. It will save variance components for each model
outVAR <- matrix(NA,ncol=4,nrow=1+2*ncol(Y))
dimnames(outVAR) <- list(c("Main",rep(paste0("Env ",colnames(Y)),2)),c("Single","Across","MxE","R-Norm"))

# Number of iterations and burn-in for Bayesian models
nIter <- 30000; burnIn <- 2000

#==============================================================
# 1. Single environment (within-environment) model, ignoring GxE effect
#==============================================================
ETA <- list(G=list(V=eigen_G$vectors,d=eigen_G$values,model='RKHS'))
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
    eigen_G1 <- MxE_eigen[[env]]
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
#==================================================
# User specifications
#==================================================
# Number of replicates
m <- 100

# Percentage of the data assigned to Testing set
percTST <- 0.3
#==================================================

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
save(YNA,file="YNA_CV1_multiEnv.RData")
```

* **Cross Validation 2 (CV2)**

Code below will generate a matrix YNA containing "NA" values for the entries corresponding to the TST set mimicing the CV2 prediction problem. It generates a 'list' with 'm' matrices containing the TRN-TST partitions

```
#==================================================
# User specifications
#==================================================
# Number of replicates
m <- 100

# Percentage of the data assigned to Testing set
percTST <- 0.3
#==================================================

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
save(YNA,file="YNA_CV2_multiEnv.RData")
```

After running the code to generate partitions for either CV1 or CV2 scenarios, the following code can be run to fit the models repeatealy for all partitions.

```
# Models
models <- c("Single","Across","MxE","R-Norm")

#==================================================
# User specifications
#==================================================
# Choose one model. 1: single; 2:across; 3:MxE; 4:R-Norm
mod <- 4

# Type of CV. 1:CV1; 2:CV2
CV <- 1
#==================================================

load(paste0("YNA_CV",CV,"_multiEnv.RData"))

m <- length(YNA)
model <- models[mod]

# List to store results. Each element will save the predictions for each partition
YHat <- list("list",m)

# Number of iterations and burn-in for Bayesian models
nIter <- 12000
burnIn <- 2000

for(k in 1:m)   # Loop for the replicates
{
    YNA0 <- YNA[[k]]
    yNA <- as.vector(YNA0)
    
    #==============================================================
    # 1. Single environment (within-environment) model, ignoring GxE effect
    #==============================================================
    if(model=="Single")
    {
        YHat0 <- matrix(NA,nrow=nrow(Y),ncol=ncol(Y))
        ETA <- list(G=list(V=eigen_G$vectors,d=eigen_G$values,model='RKHS'))
        for(env in 1:nEnv){
            fm <-BGLR(y=YNA0[,env],ETA=ETA,nIter=nIter,burnIn=burnIn)
            YHat0[,env] <- fm$yHat
        }
    }

    #==============================================================
    # 2. Across-environments model. Factor 'environment' as fixed effect
    #==============================================================
    if(model=="Across")
    {
        ETA <- list(list(~envID-1,model="FIXED"))
        ETA[[2]] <- list(V=eigen_G0$vectors,d=eigen_G0$values,model='RKHS')

        # Model Fitting
        fm <- BGLR(y=yNA,ETA=ETA,nIter=nIter,burnIn=burnIn)
        YHat0 <- matrix(fm$yHat,ncol=nEnv)
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
            eigen_G1 <- MxE_eigen[[env]]
            ETA[[(env+2)]] <- list(V=eigen_G1$vectors,d=eigen_G1$values,model='RKHS')
        }

        # Model Fitting
        fm <- BGLR(y=yNA,ETA=ETA,nIter=nIter,burnIn=burnIn)
        YHat0 <- matrix(fm$yHat,ncol=nEnv)
    }

    #==============================================================
    # 4. Reaction-Norm model. Factor 'environment' as fixed effect
    #==============================================================
    if(model=="R-Norm")
    {
        ETA <- list(list(~envID-1,model="FIXED"))
        ETA[[2]] <- list(V=eigen_G0$vectors,d=eigen_G0$values,model='RKHS')
        ETA[[3]] <- list(V=eigen_GE$vectors,d=eigen_GE$values,model="RKHS")

        # Model Fitting
        fm <- BGLR(y=yNA,ETA=ETA,nIter=nIter,burnIn=burnIn)
        YHat0 <- matrix(fm$yHat,ncol=nEnv)
    }
    YHat[[k]] <- YHat0
}

# Save results
save(YHat,file=paste0("outPRED_multiEnv_CV",CV,"_",model,".RData"))
```

#### 2.2 Retrieving results

The code below will retrieve results for all models fitted previously showing the within-environment correlation for all fitted models

```
library(ggplot2)
library(reshape)

CV <- 1
models <- c("Single","Across","MxE","R-Norm")

load(paste0("YNA_CV",CV,"_multiEnv.RData"))

# Calculate within-environment correlation
outCOR <- vector("list",length(models))
names(outCOR) <- models
for(mod in seq_along(models))
{
    filename <- paste0("outPRED_multiEnv_CV",CV,"_",models[mod],".RData")
    if(file.exists(filename)){
        load(filename,verbose=T)
        outcor <- c()
        for(k in 1:length(YHat))
        {
            YNA0 <- YNA[[k]]; YHat0 <- YHat[[k]]
            tmp <- rep(NA,ncol(YNA0))
            for(env in 1:ncol(YNA0)){
                indexTST <- which(is.na(YNA0[,env]))
                tmp[env] <- cor(Y[indexTST,env],YHat0[indexTST,env])
            }
            outcor <- rbind(outcor,tmp)
        }
        colnames(outcor) <- paste0("Env ",colnames(YNA0))
        rownames(outcor) <- NULL
        outcor <- data.frame(model=models[mod],outcor)
        outCOR[[mod]] <- outcor
    }    
}
outCOR <- outCOR[!sapply(outCOR,is.null)]

# Calculate means and SD's
t(do.call("rbind",lapply(outCOR,function(x)apply(x[,-1],2,mean))))
t(do.call("rbind",lapply(outCOR,function(x)apply(x[,-1],2,sd))))

toplot <- do.call("rbind",lapply(outCOR,function(x)melt(x,id="model")))
png(paste0("Accuacy_distn_CV",CV,".png"),height=350)
ggplot(toplot,aes(x=model,y=value,fill=variable)) + geom_boxplot()+labs(fill="Env",y="Accuracy")
dev.off()

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

##
Distribution of accuracy over 100 partitions by model. CV1

<img src="https://github.com/MarcooLopez/Genomic-Selection/blob/master/Accuacy_distn_CV1.png" width="350">


#
* **[back](https://github.com/MarcooLopez/Genomic-Selection-Demo/blob/master/README.md)**

#
# References
* Jarquín, D., Crossa, J., Lacaze, X., Du Cheyron, P., Daucourt, J., Lorgeou, J., … de los Campos, G. (2014). **A reaction norm model for genomic selection using high-dimensional genomic and environmental data**. Theoretical and Applied Genetics, 127(3), 595–607. 
* Lopez-Cruz, M., Crossa, J., Bonnett, D., Dreisigacker, S., Poland, J., Jannink, J.-L., … de los Campos, G. (2015). **Increased Prediction Accuracy in Wheat Breeding Trials Using a Marker × Environment Interaction Genomic Selection Model**. G3: Genes, Genomes, Genetics, 5(4), 569–582. 



