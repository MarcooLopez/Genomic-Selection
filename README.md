# Genomic Selection Demo

## Data
Data from CIMMYT’s Global Wheat Program. It contains information on 599 wheat lines whose grain
yield was evaluated in four environments (E1, low rainfall and
irrigated; E2, high rainfall; E3, low rainfall and high temperature;
and E4, low humidity and hot). 
Data is available for download in R-package 'BGLR' (Perez and de los Campos, 2014).

## R-packages installation
```
install.packages("BGLR")
install.packages("rrBLUP")
```

## Download data
```
library(BGLR)
data(wheat)
X <- wheat.X
Y <- wheat.Y
A <- wheat.A

# Visualize data
head(Y)
X[1:10,1:5]
```

# Data preparation
```
n <- nrow(Y)
p <- ncol(X)

# Select first environment
y <- Y[,1]

# Genomic relationship matrix
Z <- scale(X)
G <- tcrossprod(Z)/p
I <- diag(n)

# Matrix to store results
out <- matrix(NA,ncol=5,nrow=6)
dimnames(out) <- list(c(paste0("fold_",1:5),"mean"),c("GBLUP","BGBLUP","BRR","LASSO","BayesB"))
```

# Model assessment
## Cross validation with 5-folds.
Model is trained using data from 4 folds and tested in the remaining fold, this means that the model is
trained using 80% of the data and tested in the other 20%. This procedure is repeated for all the 5 folds. 

## Folds creation
```
set.seed(123)
folds <- rep(1:5,ceiling(n/5))
folds <- folds[sample(1:length(folds))]
folds <- folds[1:n]
```

# Running models
### Mixed Model G-BLUP. rrBLUP package
```
for(i in 1:5)
{
    indexTST <- which(folds==i)
    yNA <- y
    yNA[indexTST] <- NA
    fm <- mixed.solve(y=yNA,Z=I,K=G)
    out[i,1] <- cor(fm$u[indexTST],y[indexTST])
}
out[6,1] <- mean(out[1:5,1])
out[,1]
```

### Bayesian G-BLUP. BGLR package using RKHS model with K=G
```
for(i in 1:5)
{
    indexTST <- which(folds==i)
    yNA <- y
    yNA[indexTST] <- NA
    fm <- BGLR(yNA,ETA=list(list(K=G,model="RKHS")),nIter=4000,burnIn=1000)
    out[i,2] <- cor(fm$yHat[indexTST],y[indexTST])
}
out[6,2] <- mean(out[1:5,2])
out[,2]
```

### Bayesian Ridge Regression. BGLR package using keyword 'BRR'
```
for(i in 1:5)
{
    indexTST <- which(folds==i)
    yNA <- y
    yNA[indexTST] <- NA
    fm <- BGLR(yNA,ETA=list(list(X=X,model="BRR")),nIter=4000,burnIn=1000)
    out[i,3] <- cor(fm$yHat[indexTST],y[indexTST])
}
out[6,3] <- mean(out[1:5,3])
out[,3]
```

### Bayesian LASSO. BGLR package using keyword 'BL'
```
for(i in 1:5)
{
    indexTST <- which(folds==i)
    yNA <- y
    yNA[indexTST] <- NA
    fm <- BGLR(yNA,ETA=list(list(X=X,model="BL")),nIter=4000,burnIn=1000)
    out[i,4] <- cor(fm$yHat[indexTST],y[indexTST])
}
out[6,4] <- mean(out[1:5,4])
out[,4]
```

### Bayesian LASSO. BGLR package using keyword 'BayesB'
```
for(i in 1:5)
{
    indexTST <- which(folds==i)
    yNA <- y
    yNA[indexTST] <- NA
    fm <- BGLR(yNA,ETA=list(list(X=X,model="BayesB")),nIter=4000,burnIn=1000)
    out[i,5] <- cor(fm$yHat[indexTST],y[indexTST])
}
out[6,5] <- mean(out[1:5,5])
out[,5]
```



