# Genomic Selection Demo

Genomic Selection uses genetic markers covering the whole genome and potentially explaining all the genetic variance. These markers are asumed to be in Linkage Disequilibrium (LD) with the QTL thus models including all markers can estimate breeding values as linear combinatons of these QTL's.

## Model
Response variable *y* for the *i*-th individual (*i=1,...,n*) is regressed on allelic content dictated by whole-genome markers *x* through the model

![](https://latex.codecogs.com/gif.latex?y_i%3D%5Cmu&plus;%5Csum_%7Bj%3D1%7D%5Epx_%7Bij%7D%5Cbeta_j&plus;%5Cmathbf%7B%5Cvarepsilon%7D_i) 

where ![](https://latex.codecogs.com/gif.latex?%5Cmu) is the intercept, ![](https://latex.codecogs.com/gif.latex?x_%7Bij%7D) is the genotype of the *i*-th individual at the *j*-th marker, ![](https://latex.codecogs.com/gif.latex?%5Cbeta_%7Bj%7D) is the corresponding marker effect, and ![](https://latex.codecogs.com/gif.latex?%5Cmathbf%5Cvarepsilon%3D%5B%5Cvarepsilon_1%2C...%2C%5Cvarepsilon_n%5D%27) are the residuals which are assumed to be distributed Normal with constant variance ![](https://latex.codecogs.com/gif.latex?%5Cmathbf%7B%5Cvarepsilon%7D%5Csim%20N%28%5Ctextbf%7B0%7D%2CI%5Csigma%5E2_%5Cvarepsilon%29).

Model above presents some estimation difficulties when *p* is much bigger than *n* so penalization ans regularization aproaches are used to overcome this problem. Penalization and regularization solutions can be seen as posterior solutions in the Bayesian context.

### 1. Bayesian Ridge Regression (BRR).
Is a penalization regression that assumes that the regression coefficients follow a Gaussian (Normal) prior distribution, this is ![](https://latex.codecogs.com/gif.latex?%5Cbeta_%7Bj%7D%5Csim%20N%280%2C%5Csigma%5E2_%5Cbeta%29).
This prior induces shrinkage of estimates toward zero.

### 2. Bayesian LASSO.
It assumes that the regression coefficients have a prior distribution double-exponential (*DE*, or Laplace) with parameters ![](https://latex.codecogs.com/gif.latex?%5Clambda%5E2) and ![](https://latex.codecogs.com/gif.latex?%5Csigma%5E2_%5Cvarepsilon). This prior is a thick-tailed prior that can be represented as a infinite mixture of normal densities scaled by exponential (![](https://latex.codecogs.com/gif.latex?Exp)) densities, this is

![](https://latex.codecogs.com/gif.latex?%5Cbeta_j%5Csim%20DE%28%5Clambda%5E2%2C%5Csigma%5E2_%5Cvarepsilon%29%3D%5Cint%20N%28%5Cbeta_j%7C0%2C%5Csigma%5E2_%5Cvarepsilon%5Ctau%5E2_j%29Exp%28%5Ctau%5E2_j%7C%5Clambda%5E2/2%29%5Cpartial%20%5Csigma%5E2_%7B%5Cbeta_j%7D)

### 3. Bayes A
The regression effects are assumed another thick-tailed prior, a scaled *t* distribution with degree of freedom ![](https://latex.codecogs.com/gif.latex?df_%5Cbeta) and scale ![](https://latex.codecogs.com/gif.latex?S_%5Cbeta) parameters. Similar as for doble-exponential, the scaled *t* distribution is represented as mixture of normal densities scaled with a scaled-inverse Chi-squared (![](https://latex.codecogs.com/gif.latex?%5Cchi%5E%7B-1%7D)) density, this is

![](https://latex.codecogs.com/gif.latex?%5Cbeta_j%5Csim%20t%28df_%5Cbeta%2CS_%5Cbeta%29%3D%5Cint%20N%28%5Cbeta_j%7C0%2C%5Csigma%5E2_%7B%5Cbeta_j%7D%29%5Cchi%5E%7B-1%7D%28%5Csigma%5E2_%7B%5Cbeta_j%7D%7Cdf_%5Cbeta%2CS_%5Cbeta%29%5Cpartial%20%5Csigma%5E2_%7B%5Cbeta_j%7D)

### 4. Bayes B
Markers effects are asummed to be equal to zero with probability ![](https://latex.codecogs.com/gif.latex?%5Cpi) and with probability 1-![](https://latex.codecogs.com/gif.latex?%5Cpi) are assumed to follow a scaled *t* distribution as in Bayes A model.

### 4. G-BLUP model (RR-BLUP)
The response is modeled as ![](https://latex.codecogs.com/gif.latex?y_i%3D%5Cmu&plus;u_i&plus;%5Cvarepsilon_i) and its solution is equivalent to that of the BRR model arised when in the model above we make the sustitution 

![](https://latex.codecogs.com/gif.latex?u_i%3D%5Csum_%7Bj%3D1%7D%5Epx_%7Bij%7D%5Cbeta_j)

It can be shown that the random vector ![](https://latex.codecogs.com/gif.latex?%5Ctextbf%7Bu%7D%3D%5Bu_1%2C...%2Cu_n%5D%27) follows a Normal distribution ![](https://latex.codecogs.com/gif.latex?%5Ctextbf%7Bu%7D%5Csim%20N%28%5Ctextbf%7B0%7D%2C%5Ctextbf%7BG%7D%5Csigma%5E2_u%29), where ![](https://latex.codecogs.com/gif.latex?%5Ctextbf%7BG%7D%3D%5Ctextbf%7BX%7D%5Ctextbf%7BX%7D%27/p) with ***X*** is the matrix of centered and standardized marker genotypes and it is called genomic relationship matrix.



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


