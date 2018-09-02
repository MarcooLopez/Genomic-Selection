# Genomic Selection Demo

Genomic Selection uses genetic markers covering the whole genome and potentially explaining all the genetic variance. These markers are asumed to be in Linkage Disequilibrium (LD) with the QTL thus models including all markers can estimate breeding values as linear combinatons of these QTL's.

## Model
Response variable *y* for the *i*-th individual (*i=1,...,n*) is regressed on allelic content dictated by whole-genome markers *x* through the model

![](https://latex.codecogs.com/gif.latex?y_i%3D%5Cmu&plus;%5Csum_%7Bj%3D1%7D%5Epx_%7Bij%7D%5Cbeta_j&plus;%5Cmathbf%7B%5Cvarepsilon%7D_i) 

where ![](https://latex.codecogs.com/gif.latex?%5Cmu) is the population mean, ![](https://latex.codecogs.com/gif.latex?x_%7Bij%7D) is the genotype of the *i*-th individual at the *j*-th marker, ![](https://latex.codecogs.com/gif.latex?%5Cbeta_%7Bj%7D) is the corresponding marker effect, and ![](https://latex.codecogs.com/gif.latex?%5Cmathbold%7B%5Cvarepsilon%7D) are the residuals which are assumed to be distributed Normal with constant variance ![](https://latex.codecogs.com/gif.latex?%5Cmathbf%7B%5Cvarepsilon%7D%5Csim%20N%28%5Ctextbf%7B0%7D%2CI%5Csigma%5E2_%5Cvarepsilon%29).

Model above presents some estimation difficulties when *p* is much bigger than *n* so penalization ans regularization aproaches are used to overcome this problem. Penalization and regularization solutions can be seen as posterior solutions in the Bayesian context.

### 1. Bayesian Ridge Regression.
Is a penalization regression that assumes that the regression coefficients follow a Gaussian (Normal) prior distribution, this is ![](https://latex.codecogs.com/gif.latex?%5Cbeta_%7Bj%7D%5Csim%20N%280%2C%5Csigma%5E2_%5Cbeta%29).
This prior induces shrinkage of estimates toward zero.

### 2. Bayesian LASSO.
It assumes that the regression coefficients have a prior distribution double-exponential (or Laplace)


### 1. G-BLUP model
It is equivalent to the solution of the BRR model when in the model above we make the sustitution ![](https://latex.codecogs.com/gif.latex?u_i%3D%5Csum_%7Bj%3D1%7D%5Epx_%7Bij%7D%5Cbeta_j). It can be shown that the random vector ![](https://latex.codecogs.com/gif.latex?%5Ctextbf%7Bu%7D%3D%5Bu_1%2C...%2Cu_n%5D%27) follows a Normal distribution 



## Data
Data from CIMMYT’s Global Wheat Program. It contains information on 599 wheat lines whose grain
yield was evaluated in four environments (E1, low rainfall and
irrigated; E2, high rainfall; E3, low rainfall and high temperature;
and E4, low humidity and hot). 
Data is available for download in R-package 'BGLR' (Perez and de los Campos, 2014).

## R-packages installation
```
if(!"BGLR"%in%rownames(installed.packages()))  install.packages("BGLR")
if(!"rrBLUP"%in%rownames(installed.packages())) install.packages("rrBLUP")
library(BGLR)
library(rrBLUP)
```

## Download data
```
data(wheat)
X <- wheat.X
Y <- wheat.Y
A <- wheat.A

# Visualize data
head(Y)
X[1:10,1:5]
```

## Analyses
* **[Single-environment](https://github.com/MarcooLopez/Genomic-Selection-Demo/blob/master/single_environment.md)**
* **[Multi-environment](https://github.com/MarcooLopez/Genomic-Selection-Demo/blob/master/multi_environment.md)**


