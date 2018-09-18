# Genomic Selection Demo

The standard genetic model assumes that phenotype is the sum of a genetic component and a non-genetic component (residual), ![](https://latex.codecogs.com/gif.latex?y_i%3Dg_i&plus;%5Cvarepsilon_i). Genomic Selection uses genetic markers covering the whole genome and potentially explaining all the genetic variance. These markers are asumed to be in Linkage Disequilibrium (LD) with the QTL thus models including all markers can estimate breeding values ![](https://latex.codecogs.com/gif.latex?g_i) as combinatons of these QTL's.

## Model
Response variable *y* for the *i*-th individual (*i=1,...,n*) is regressed on a function of *p* marker genotypes ![](https://latex.codecogs.com/gif.latex?%5Ctextbf%7Bx%7D_i%3D%5Bx_%7Bi1%7D%2C...%2Cx_%7Bip%7D%5D%27) that seeks to aproximate to the true genetic value of the individual, this is

![](https://latex.codecogs.com/gif.latex?y_i%3Df%28%5Ctextbf%7Bx%7D_i%29&plus;%5Cvarepsilon_i),

where function ![](https://latex.codecogs.com/gif.latex?f%28%5Ctextbf%7Bx%7D_i%29) can be a parametric or non-parametric and ![](https://latex.codecogs.com/gif.latex?%5Cmathbf%5Cvarepsilon%3D%5B%5Cvarepsilon_1%2C...%2C%5Cvarepsilon_n%5D%27) are the residuals which are usually assumed to be distributed Normal with constant variance ![](https://latex.codecogs.com/gif.latex?%5Cmathbf%7B%5Cvarepsilon%7D%5Csim%20N%28%5Ctextbf%7B0%7D%2CI%5Csigma%5E2_%5Cvarepsilon%29).

### Parametric regression
The genotypic value of an individual is estimated using a **linear model** in which a linear combination of the marker genotypes are used, that is

![](https://latex.codecogs.com/gif.latex?f%28%5Ctextbf%7Bx%7D_i%29%3D%5Cmu&plus;%5Csum_%7Bj%3D1%7D%5Epx_%7Bij%7D%5Cbeta_j) 

where ![](https://latex.codecogs.com/gif.latex?%5Cmu) is the intercept, ![](https://latex.codecogs.com/gif.latex?x_%7Bij%7D) is the genotype of the *i*-th individual at the *j*-th marker, ![](https://latex.codecogs.com/gif.latex?%5Cbeta_%7Bj%7D) is the corresponding marker effect.

Model above presents some estimation difficulties when *p* is much bigger than *n* so penalization ans regularization aproaches are used to overcome this problem. Penalization and regularization solutions can be seen as posterior solutions in the Bayesian context.

#### 1. Bayesian Ridge Regression (BRR).
Is a penalization regression that assumes that the regression coefficients follow a Gaussian (Normal) prior distribution, this is ![](https://latex.codecogs.com/gif.latex?%5Cbeta_%7Bj%7D%5Csim%20N%280%2C%5Csigma%5E2_%5Cbeta%29).
This prior induces shrinkage of estimates toward zero.

#### 2. Bayesian LASSO.
It assumes that the regression coefficients have a prior distribution double-exponential (*DE*, or Laplace) with parameters ![](https://latex.codecogs.com/gif.latex?%5Clambda%5E2) and ![](https://latex.codecogs.com/gif.latex?%5Csigma%5E2_%5Cvarepsilon). This prior is a thick-tailed prior that can be represented as a infinite mixture of normal densities scaled by exponential (![](https://latex.codecogs.com/gif.latex?Exp)) densities, this is

![](https://latex.codecogs.com/gif.latex?%5Cbeta_j%5Csim%20DE%28%5Clambda%5E2%2C%5Csigma%5E2_%5Cvarepsilon%29%3D%5Cint%20N%28%5Cbeta_j%7C0%2C%5Csigma%5E2_%5Cvarepsilon%5Ctau%5E2_j%29Exp%28%5Ctau%5E2_j%7C%5Clambda%5E2/2%29%5Cpartial%20%5Csigma%5E2_%7B%5Cbeta_j%7D)

#### 3. Bayes A
The regression effects are assumed another thick-tailed prior, a scaled *t* distribution with degree of freedom ![](https://latex.codecogs.com/gif.latex?df_%5Cbeta) and scale ![](https://latex.codecogs.com/gif.latex?S_%5Cbeta) parameters. Similar as for doble-exponential, the scaled *t* distribution is represented as mixture of normal densities scaled with a scaled-inverse Chi-squared (![](https://latex.codecogs.com/gif.latex?%5Cchi%5E%7B-1%7D)) density, this is

![](https://latex.codecogs.com/gif.latex?%5Cbeta_j%5Csim%20t%28df_%5Cbeta%2CS_%5Cbeta%29%3D%5Cint%20N%28%5Cbeta_j%7C0%2C%5Csigma%5E2_%7B%5Cbeta_j%7D%29%5Cchi%5E%7B-1%7D%28%5Csigma%5E2_%7B%5Cbeta_j%7D%7Cdf_%5Cbeta%2CS_%5Cbeta%29%5Cpartial%20%5Csigma%5E2_%7B%5Cbeta_j%7D)

#### 4. Bayes B
Markers effects are asummed to be equal to zero with probability ![](https://latex.codecogs.com/gif.latex?%5Cpi) and with probability 1-![](https://latex.codecogs.com/gif.latex?%5Cpi) are assumed to follow a scaled *t* distribution as in Bayes A model.

#### 5. Bayes C
Similar to Bayes B, markers effects are asummed to be equal to zero with probability ![](https://latex.codecogs.com/gif.latex?%5Cpi) and with probability 1-![](https://latex.codecogs.com/gif.latex?%5Cpi) are assumed to follow a Gaussian distribution as in BRR model.

#### 6. G-BLUP model (RR-BLUP)
The response is modeled as ![](https://latex.codecogs.com/gif.latex?y_i%3D%5Cmu&plus;u_i&plus;%5Cvarepsilon_i) and its solution is equivalent to that of the BRR model arised when in the model above we make the sustitution 

![](https://latex.codecogs.com/gif.latex?g_i%3D%5Csum_%7Bj%3D1%7D%5Epx_%7Bij%7D%5Cbeta_j)

It can be shown that the random vector ![](https://latex.codecogs.com/gif.latex?%5Ctextbf%7Bg%7D%3D%5Bg_1%2C...%2Cg_n%5D%27) follows a Normal distribution ![](https://latex.codecogs.com/gif.latex?%5Ctextbf%7Bg%7D%5Csim%20N%28%5Ctextbf%7B0%7D%2C%5Ctextbf%7BG%7D%5Csigma%5E2_g%29), where ![](https://latex.codecogs.com/gif.latex?%5Ctextbf%7BG%7D%3D%5Ctextbf%7BX%7D%5Ctextbf%7BX%7D%27/p) with ***X*** is the matrix of centered and standardized marker genotypes and it is called genomic relationship matrix.


### Semi-parametric regression

#### 7. RKHS regression.
The genomic function ![](https://latex.codecogs.com/gif.latex?f%28%5Ctextbf%7Bx%7D_i%29) is expressed as a linear combination of some positive semi-definite basis functions called Reproducing Kernels (RK), ![](https://latex.codecogs.com/gif.latex?K%28%5Ctextbf%7Bx%7D_i%2C%5Ctextbf%7Bx%7D_%7Bi%27%7D%29), as follows

![](https://latex.codecogs.com/gif.latex?f%28%5Ctextbf%7Bx%7D_i%29%3D%5Csum_%7Bi%27%3D1%7D%5EnK%28%5Ctextbf%7Bx%7D_i%2C%5Ctextbf%7Bx%7D_%7Bi%27%7D%29%5Calpha_%7Bi%27%7D)

This model can be rewritten as ![](https://latex.codecogs.com/gif.latex?%5Ctextbf%7By%7D%3D%5Ctextbf%7BK%7D%5Cmathbf%7B%5Calpha%7D&plus;%5Cmathbf%7B%5Cvarepsilon%7D) where ![](https://latex.codecogs.com/gif.latex?%5Ctextbf%7BK%7D%3D%5C%7BK%28%5Ctextbf%7Bx%7D_i%2C%5Ctextbf%7Bx%7D_%7Bi%27%7D%29%5C%7D) is a ![](https://latex.codecogs.com/gif.latex?n%5Ctimes%20n) matrix containing all the evaluations of the RK function at the point (*i*,*i'*) and ![](https://latex.codecogs.com/gif.latex?%5Cmathbf%7B%5Calpha%7D%3D%5B%5Calpha_1%2C...%2C%5Calpha_n%5D%27).

This problem can be solved in a Bayesian fashion by assuming a prior ![](https://latex.codecogs.com/gif.latex?%5Cmathbf%7B%5Calpha%7D%5Csim%20N%28%5Ctextbf%7B0%7D%2C%5Ctextbf%7BK%7D%5E%7B-1%7D%5Csigma%5E2_%5Calpha%29).

**Note:** The Ridge Regression (and consequently, G-BLUP) can be represented as a RKHS model by setting **K**=**G**.  


## Implementation of models
Models previously above described will be implemented in R software using R-packages 'BGLR' and 'rrBLUP'. Using public data, it will be shown how to run the models for the single-environment case and then how to perform a multi-environment analysis with the G-BLUP model using a marker-by-environment (MxE) and a Reaction Norm approaches that account for GxE interaction.

### Data
Data from CIMMYT’s Global Wheat Program. Lines were evaluated for grain yield (each entry corresponds to an average of two plot records) at four different environments; phenotypes (*wheat.Y* object) were centered and standardized to a unit variance within environment. Each of the lines were genotyped for 1279 diversity array technology (DArT) markers. At each marker two homozygous genotypes were possible and these were coded as 0/1. Marker genotypes are given in the object *wheat.X*. Finally a matrix *wheat.A* provides the pedigree relationships between lines computed from the pedigree records.
Data is available for download in the R-package 'BGLR'.

### R-packages installation
```
if(!"BGLR"%in%rownames(installed.packages()))  install.packages("BGLR")
if(!"rrBLUP"%in%rownames(installed.packages())) install.packages("rrBLUP")
library(BGLR)
library(rrBLUP)
```

### Download data
```
data(wheat)
X <- wheat.X
Y <- wheat.Y
A <- wheat.A

# Visualize data
head(Y)
X[1:10,1:5]
```

### Type of analyses
* **[Single-environment](https://github.com/MarcooLopez/Genomic-Selection-Demo/blob/master/single_environment.md)**
* **[Multi-environment](https://github.com/MarcooLopez/Genomic-Selection-Demo/blob/master/multi_environment.md)**

##

# References
* de los Campos, G., Gianola, D., Rosa, G. J. M., Weigel, K. A., & Crossa, J. (2010). **Semi-parametric genomic-enabled prediction of genetic values using reproducing kernel Hilbert spaces methods**. Genetics Research, 92(4), 295–308. 
* de los Campos, G., Hickey, J. M., Pong-Wong, R., Daetwyler, H. D., & Calus, M. P. L. (2013). **Whole-genome regression and prediction methods applied to plant and animal breeding**. Genetics, 193(2), 327–345.
* Endelman, J. B. (2011). **Ridge Regression and Other Kernels for Genomic Selection with R Package rrBLUP**. The Plant Genome Journal, 4(3), 250–255. 
* Habier, D., Fernando, R. L., Kizilkaya, K., & Garrick, D. J. (2011). **Extension of the bayesian alphabet for genomic selection**. BMC Bioinformatics, 12(186), 1-12.
* Jarquín, D., Crossa, J., Lacaze, X., Du Cheyron, P., Daucourt, J., Lorgeou, J., … de los Campos, G. (2014). **A reaction norm model for genomic selection using high-dimensional genomic and environmental data**. Theoretical and Applied Genetics, 127(3), 595–607.
* Lopez-Cruz, M., Crossa, J., Bonnett, D., Dreisigacker, S., Poland, J., Jannink, J.-L., … de los Campos, G. (2015). **Increased prediction accuracy in wheat breeding trials using a marker × environment interaction genomic selection model**. G3: Genes, Genomes, Genetics, 5(4), 569–582.
* Meuwissen, T. H. E., Hayes, B. J., & Goddard, M. E. (2001). **Prediction of total genetic value using genome-wide dense marker maps**. Genetics, 157(4), 1819–1829.
* Park, T., & Casella, G. (2008). **The Bayesian Lasso**. Journal of the American Statistical Association, 103(482), 681–686.
* Perez, P., & de los Campos, G. (2014). **Genome-wide regression and prediction with the BGLR statistical package**. Genetics, 198(2), 483–495. 
R Development Core Team. (2015). **R: A Language and Environment for Statistical Computing**. Vienna, Austria: R Foundation for Statistical Computing.


