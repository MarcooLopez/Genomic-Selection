# Data preparation
```
n <- nrow(Y)
p <- ncol(X)

# Genomic relationship matrix
Z <- scale(X)
G <- tcrossprod(Z)/p
I <- diag(n)

# Select only environments 2,4, and 5 to work with
Y0 <- Y[,c(2,3,4)]

# Matrix to store results
out <- matrix(NA,nrow=ncol(Y0),ncol=3)
out[,1] <- colnames(Y0)
colnames(out) <- c("Environment","Within","MxE")
```

# MxE model
Model GxE interaction using a marker x environment (MxE) approach that benefits of positively correlated environments 
and can be used in two predictions problems (CV1 and CV2) that mimic 2 evaluation situations:
1. Cross Validation 1 (CV1). Represent a scheme of prediction of lines that have not been evaluated in any field
trials.
2. Cross Validation 2 (CV2). Represent a scheme of prediction of lines that have been evaluated in some but all target environments.

Reference: *[Lopez-Cruz et. al, 2015](https://www.ncbi.nlm.nih.gov/pubmed/25660166)*

# Model assessment
## Training-Testing random partitions.
Data is randomly splitted into training and testing using 70% of the data for training and 30%  for testing as depicted in the figures below.

![picture alt](https://github.com/MarcooLopez/Genomic-Selection-Demo/blob/master/CV1_2_scheme.png)

This procedure of TRN-TST can be repeated many times to allow for estimation of standard errors (SE)

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
    fm <- BGLR(yNA,ETA=list(list(K=G,model="RKHS")),nIter=4000,burnIn=1000)
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
    fm <- BGLR(yNA,ETA=list(list(X=X,model="BRR")),nIter=4000,burnIn=1000)
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
    fm <- BGLR(yNA,ETA=list(list(X=X,model="BL")),nIter=4000,burnIn=1000)
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
    fm <- BGLR(yNA,ETA=list(list(X=X,model="BayesB")),nIter=4000,burnIn=1000)
    out[i,5] <- cor(fm$yHat[indexTST],y[indexTST])
}
out[6,5] <- mean(out[1:5,5])
out[,5]
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
