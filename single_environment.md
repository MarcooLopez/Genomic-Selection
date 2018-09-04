
# Model assessment
## Training-Testing random partitions
The prediction power of the model will be assessed using the training-testing (TRN-TST) random partitions approach. Data is randomly splitted into training and testing sets. Model parameters are estimated in training set and model is tested in TST set. 

# Data preparation
### Load data, generate G-matrix and create objects to store results
```
X <- wheat.X
Y <- wheat.Y

n <- nrow(Y)
p <- ncol(X)

# Select environment 2
y <- Y[,2]

# Genomic relationship matrix
Z <- scale(X)
G <- tcrossprod(Z)/p
I <- diag(n)
```

# Running models

## 1. One single instance of 5-folds partition
The following code runs a single 5-folds partition for each model

```
# Creation of folds
set.seed(123)
folds <- rep(1:5,ceiling(n/5))
folds <- folds[sample(1:length(folds))]
folds <- folds[1:n]

# Matrix to store results
out <- matrix(NA,ncol=5,nrow=6)
dimnames(out) <- list(c(paste0("fold_",1:5),"mean"),c("GBLUP","BGBLUP","BRR","LASSO","BayesB"))

# Number of iterations and burn-in for Bayesian models
nIter <- 2000
burnIn <- 500

for(i in 1:5)   # Loop for the 5 folds
{
    indexTST <- which(folds==i)
    yNA <- y
    yNA[indexTST] <- NA
    
    # G-BLUP model using rrBLUP package
    fm <- mixed.solve(y=yNA,Z=I,K=G)
    out[i,1] <- cor(fm$u[indexTST],y[indexTST])
    
    # G-BLUP (Bayesian) model using BGLR package. RKHS model with K=G
    fm <- BGLR(yNA,ETA=list(list(K=G,model="RKHS")),nIter=nIter,burnIn=burnIn)
    out[i,2] <- cor(fm$yHat[indexTST],y[indexTST])
    
    # Bayesian Ridge Regression using BGLR package.
    fm <- BGLR(yNA,ETA=list(list(X=X,model="BRR")),nIter=nIter,burnIn=burnIn)
    out[i,3] <- cor(fm$yHat[indexTST],y[indexTST])
    
    # Bayesian LASSO model using BGLR package.
    fm <- BGLR(yNA,ETA=list(list(X=X,model="BL")),nIter=nIter,burnIn=burnIn)
    out[i,4] <- cor(fm$yHat[indexTST],y[indexTST])
    
    # Bayes B model using BGLR package.
    fm <- BGLR(yNA,ETA=list(list(X=X,model="BayesB")),nIter=nIter,burnIn=burnIn)
    out[i,5] <- cor(fm$yHat[indexTST],y[indexTST])
}

# Take the average across folds for each model
out[6,] <- apply(out[1:5,],2,mean)
print(out)
```

### Results

|       |GBLUP  |B-GBLUP | BRR  | LASSO | Bayes B |
|-------|-------|--------|------|-------|-------|
|Fold 1  | 0.52  | 0.52  | 0.51 | 0.52 | 0.52 |
|Fold 2  | 0.53  | 0.53  | 0.54 | 0.53 | 0.53 |
|Fold 3  | 0.53  | 0.53  | 0.54 | 0.54 | 0.53 |
|Fold 4  | 0.49  | 0.49  | 0.49 | 0.48 | 0.49 |
|Fold 5  | 0.53  | 0.53  | 0.54 | 0.54 | 0.53 |
|Mean    | 0.52  | 0.52  | 0.52 | 0.52 | 0.52 |

#
## 2. Replicates of partitions to obtain standard deviations of predictions
The code below runs repeated partitions to obtain mean and standard deviations of accuracies for a single model.

All the models will be run using 'BGLR' package.

```
# Models
models <- c("GBLUP","BRR","LASSO","BayesB")

# Percentage of the data assigned to Testing set
percTST <- 0.3

# Number of replicates
m <- 100

# Creation of seed for repeated randomizations
set.seed(123)
seeds <- round(seq(1E3,1E6,length=m))

# Matrix to store results. It will save the corelation for each partition
outCOR <- matrix(NA,nrow=m,ncol=length(models))
colnames(outCOR) <- models

# Number of iterations and burn-in for Bayesian models
nIter <- 200
burnIn <- 50

nTST <- round(percTST*n)

for(k in 1:m)   # Loop for the replicates
{
    set.seed(seeds[k])
    indexTST <- sample(1:n,size=nTST,replace=FALSE)
    yNA <- y
    yNA[indexTST] <- NA

    for(mod in 1:length(models))   # Loop for the models
    {
        model <- models[mod]
        
        if(model=="GBLUP")  ETA <- list(list(K=G,model="RKHS"))
        if(model=="BRR")    ETA <- list(list(X=X,model="BRR"))
        if(model=="LASSO")  ETA <- list(list(X=X,model="BL"))
        if(model=="BayesB") ETA <- list(list(X=X,model="BayesB"))

        fm <- BGLR(yNA,ETA=ETA,nIter=nIter,burnIn=burnIn)
        outCOR[k,mod] <- cor(fm$yHat[indexTST],y[indexTST])
    }
}
```

### Results
```
rbind(mean=apply(outCOR,2,mean),sd=apply(outCOR,2,sd))
boxplot(outCOR,ylab="Accuracy",xlab="Model")
```

Boxplot of distribution of the accuracies by model

![](https://github.com/MarcooLopez/Genomic-Selection/blob/master/boxplot1.png)

#
* **[back](https://github.com/MarcooLopez/Genomic-Selection-Demo/blob/master/README.md)**
