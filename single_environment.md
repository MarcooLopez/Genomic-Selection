
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

## 1. Replicates of partitions to obtain standard deviations of predictions
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
nIter <- 12000
burnIn <- 2000

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
round(rbind(mean=apply(outCOR,2,mean),sd=apply(outCOR,2,sd)),4)
boxplot(outCOR,ylab="Accuracy",xlab="Model")
```

|       |GBLUP  | BRR  | LASSO | Bayes B |
|-------|-------|--------|------|-------|
|Mean | 0.475  | 0.476  | 0.474 | 0.465 |
|SD  | 0.049 | 0.050 | 0.049 | 0.048 |

Boxplot of distribution of the accuracies by model

![](https://github.com/MarcooLopez/Genomic-Selection/blob/master/boxplot1.png)

#
* **[back](https://github.com/MarcooLopez/Genomic-Selection-Demo/blob/master/README.md)**
