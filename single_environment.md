# Data preparation
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

# Matrix to store results
out <- matrix(NA,ncol=5,nrow=6)
dimnames(out) <- list(c(paste0("fold_",1:5),"mean"),c("GBLUP","BGBLUP","BRR","LASSO","BayesB"))
```

# Model assessment
## Cross validation with 5-folds.
To mimic prediction of GEBV of new untested breeding material, data is randomly splitted into 5 sets where training set is comprised of any 4 folds and testing set will consist of the remaining fold. This means that the model is
trained using 80% of the data and tested in the other 20%. This procedure is repeated for all the 5 folds.

## Folds creation
```
set.seed(123)
folds <- rep(1:5,ceiling(n/5))
folds <- folds[sample(1:length(folds))]
folds <- folds[1:n]
```

# Running models
### 1. Mixed Model G-BLUP. rrBLUP package
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

### 2. Bayesian G-BLUP. BGLR package using RKHS model with K=G
```
for(i in 1:5)
{
    indexTST <- which(folds==i)
    yNA <- y
    yNA[indexTST] <- NA
    fm <- BGLR(yNA,ETA=list(list(K=G,model="RKHS")),nIter=12000,burnIn=2000)
    out[i,2] <- cor(fm$yHat[indexTST],y[indexTST])
}
out[6,2] <- mean(out[1:5,2])
out[,2]
```

### 3. Bayesian Ridge Regression. BGLR package using keyword 'BRR'
```
for(i in 1:5)
{
    indexTST <- which(folds==i)
    yNA <- y
    yNA[indexTST] <- NA
    fm <- BGLR(yNA,ETA=list(list(X=X,model="BRR")),nIter=12000,burnIn=2000)
    out[i,3] <- cor(fm$yHat[indexTST],y[indexTST])
}
out[6,3] <- mean(out[1:5,3])
out[,3]
```

### 4. Bayesian LASSO. BGLR package using keyword 'BL'
```
for(i in 1:5)
{
    indexTST <- which(folds==i)
    yNA <- y
    yNA[indexTST] <- NA
    fm <- BGLR(yNA,ETA=list(list(X=X,model="BL")),nIter=12000,burnIn=2000)
    out[i,4] <- cor(fm$yHat[indexTST],y[indexTST])
}
out[6,4] <- mean(out[1:5,4])
out[,4]
```

### 4. Bayesian B. BGLR package using keyword 'BayesB'
```
for(i in 1:5)
{
    indexTST <- which(folds==i)
    yNA <- y
    yNA[indexTST] <- NA
    fm <- BGLR(yNA,ETA=list(list(X=X,model="BayesB")),nIter=12000,burnIn=2000)
    out[i,5] <- cor(fm$yHat[indexTST],y[indexTST])
}
out[6,5] <- mean(out[1:5,5])
out[,5]
```
## Results

|       |GBLUP  |B-GBLUP | BRR  | LASSO | Bayes B |
|-------|-------|--------|------|-------|-------|
|Fold 1  | 0.52  | 0.52  | 0.51 | 0.52 | 0.52 |
|Fold 2  | 0.53  | 0.53  | 0.54 | 0.53 | 0.53 |
|Fold 3  | 0.53  | 0.53  | 0.54 | 0.54 | 0.53 |
|Fold 4  | 0.49  | 0.49  | 0.49 | 0.48 | 0.49 |
|Fold 5  | 0.53  | 0.53  | 0.54 | 0.54 | 0.53 |
|Mean    | 0.52  | 0.52  | 0.52 | 0.52 | 0.52 |

![](https://latex.codecogs.com/gif.latex?%5Cfrac%7B1%7D%7B2%7De%5E%7B-%28x-%5Cmu%29%5E2%7D)
<img src="https://latex.codecogs.com/gif.latex?%5Cfrac%7B1%7D%7B2%7De%5E%7B-%28x-%5Cmu%29%5E2%7D">

#
* **[back](https://github.com/MarcooLopez/Genomic-Selection-Demo/blob/master/README.md)**
