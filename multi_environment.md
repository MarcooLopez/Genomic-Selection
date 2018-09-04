
# Multi-environment models

### * Reaction Norm model
Is an extention of the G-BLUP model that incorporates GxE by introducing covariance structures as a funcion of the marker information.
The reaction norm model can model main and interaction effects of environmental covariates (EC) and markers.


### * MxE model
It models GxE interaction using a marker-by-environment (MxE) approach that benefits of positively correlated environments. MxE descomposes marker effects into an effect that is common to all environments and an effect that is specific to each environment.


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
X <- wheat.X
Y <- wheat.Y

n <- nrow(Y)
p <- ncol(X)

# Select environments. For instance, environments 2,4, and 5
Y <- Y[,c(2,3,4)]

# Genomic relationship matrix
Z <- scale(X)
G <- tcrossprod(Z)/p

# Z matrix for individuals. In this case it is a diagonal since there are no replicates
GID <- factor(rownames(Y),levels=rownames(Y))
Z <- model.matrix(~GID-1)
```

## Running models
### 1. Variance components estimation 
```
```

### 2. Replicates of partitions to obtain standard deviations of predictions
Using a GBLUP approach, the prediction power of the multi-environment models (MxE and Reaction Norm) will be compared with that that ignores GxE (across-environment) and with the GBLUP model fitted within environment. 

#### 2.1. Training-Testing partitions
User can choose either to perform CV1 or CV2 aproaches

* **Cross Validation 1 (CV1)**

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

* **Cross Validation 2 (CV2)**

Code below will generate a matrix YNA containing "NA" values for the entries corresponding to the TST set mimicing the CV2 prediction problem. 

```
set.seed(123)

env <- c(1,2,3) # choose any set of environments from 1:ncol(Y)
nEnv <- length(env)
Y <- Y[,env]
n <- nrow(Y)

percTST <- 0.3
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

#### Single environment models
**1. Within-environment model, ignoring GxE effect**

```
YHat1 <- matrix(nrow=nrow(Y),ncol=ncol(Y),NA)
ETA <- list(G=list(K=G,model='RKHS'))
for(j in 1:nEnv){
    prefix <- paste(colnames(Y)[j],"_",sep="")
    fm1 <-BGLR(y=YNA[,j],ETA=ETA,nIter=12000,burnIn=2000,saveAt=prefix)
    YHat1[,j] <- fm1$yHat
}
```

#### Multi environment models
The environment as fixed effect and main effects of markers will be common to all multi-environment models.
Eigen value decompostion (EVD) will be used instead of the whole matrix to make speed computational time.

```
yNA <- as.vector(YNA)

# Fixed effect (environment-intercepts)
envID <- factor(rep(colnames(Y),each=nrow(Y)),levels=colnames(Y))

# Main effects of markers
GID <- factor(rep(rownames(Y),ncol(Y)),levels=rownames(Y))
Zg <- model.matrix(~GID-1)
G0 <- Zg%*%G%*%t(Zg)
eigen_G0 <- eigen(G0)
```

**2. Accros-environment**

Including only main effect that is common to all environments but ignoring GxE interaction.

```
# Fixed effect and main effects of markers
ETA <- list(list(~envID-1,model="FIXED"))
ETA[[2]] <- list(V=eigen_G0$vectors,d=eigen_G0$values,model='RKHS')

# Model Fitting
fm2 <- BGLR(y=yNA,ETA=ETA,nIter=12000,burnIn=2000,saveAt="AcrossEnv_")
YHat2 <- matrix(fm2$yHat,ncol=nEnv)
```

**3. MxE Interaction**

Including main effect that is common to all environments and an environment-specific effect (interaction) that accounts for GxE.

```
# Fixed effect and main effects of markers
ETA <- list(list(~envID-1,model="FIXED"))
ETA[[2]] <- list(V=eigen_G0$vectors,d=eigen_G0$values,model='RKHS')

# Adding interaction terms
for(j in 1:nEnv){
    tmp <- rep(0,nEnv) ; tmp[j] <- 1; G1 <- kronecker(diag(tmp),G)
    eigen_G1 <- eigen(G1)
    ETA[[(j+2)]] <- list(V=eigen_G1$vectors,d=eigen_G1$values,model='RKHS')
}

# Model Fitting
fm3 <- BGLR(y=yNA,ETA=ETA,nIter=12000,burnIn=2000,saveAt="MxE_")
YHat3 <- matrix(fm3$yHat,ncol=nEnv)
```

**4. Reaction norm model**

Including main effect that is common to all environments and modeling GxE using a covariance structure.

```
# Fixed effect and main effects of markers
ETA <- list(list(~envID-1,model="FIXED"))
ETA[[2]] <- list(V=eigen_G0$vectors,d=eigen_G0$values,model='RKHS')

# Incidence matrix for the effects of environments
ZE <- model.matrix(~envID-1)   
ZEZEt <- tcrossprod(ZE)

# Adding GxE interaction term as Hadamart product
GE <- G0*ZEZEt
eigen_GE <- eigen(GE)
ETA[[3]] <- list(V=eigen_GE$vectors,d=eigen_GE$values,model="RKHS")

# Model Fitting
fm4 <- BGLR(y=yNA,ETA=ETA,nIter=12000,burnIn=2000,saveAt="RNorm_")
YHat4 <- matrix(fm4$yHat,ncol=nEnv)
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



