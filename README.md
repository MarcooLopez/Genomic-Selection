# Genomic Selection Demo

Genomic Selection uses genetic markers covering the whole genome and potentially explaining all the genetic variance. These markers are asumed to be in Linkage Disequilibrium (LD) with the QTL thus models including all markers can estimate breeding values as linear combinatons of these QTL's.

## Model
Response variable 'y' is regressed on allelic content dictated by molecular markers 'X'. Breeding values are estimated using information from markers through the model

![](https://latex.codecogs.com/gif.latex?%5Ctextbf%7By%7D%3DX%5Cmathbf%7B%5Cbeta%7D&plus;Z%5Cmathbf%7Bu%7D&plus;%5Cmathbf%7B%5Cvarepsilon%7D)  

where 'X' is a full-rank desig matrix for the fixed effects ![](https://latex.codecogs.com/gif.latex?%5Cmathbf%7B%5Cbeta%7D), Z is the 
design matrix for the random effects ![](https://latex.codecogs.com/gif.latex?%5Cmathbf%7Bu%7D), and the residuals are normal with constant variance.

### 1. G-BLUP model
Random effects are assumed to follow a Normal distribution
![](https://latex.codecogs.com/gif.latex?%5Ctextbf%7Bu%7D%5Csim%20N%280%2C%5Csigma%5E2_u%5Ctextbf%7BK%7D%29),
where K is a positive semidefinite matrix representing the relationship among individuals.


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


