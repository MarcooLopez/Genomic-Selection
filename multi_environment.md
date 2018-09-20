
# Multi-environment models

* **Single-environment model**

This model is obtained by regressing the phenotype vector containing the *n* records available in the *k*<sup>th</sup> environment (*k=1,2,...,s* environments), ![](https://latex.codecogs.com/gif.latex?%5Ctextbf%7By%7D_k%3D%5By_%7B1k%7D%2C...%2Cy_%7Bnk%7D%5D%27), where *i=1,2,...,n* indexes lines (individuals), on *p* markers using a linear model in the form

<p align="center">
<img src="https://latex.codecogs.com/gif.latex?%5Ctextbf%7By%7D_k%3D%5Cmu_k%5Ctextbf%7B1%7D&plus;%5Ctextbf%7BX%7D_k%5Cboldsymbol%7B%5Cbeta%7D_k&plus;%5Cboldsymbol%7B%5Cvarepsilon%7D_k">
</p>

where ![](https://latex.codecogs.com/gif.latex?%5Cmu_k) is the intercept, ![](https://latex.codecogs.com/gif.latex?%5Ctextbf%7B1%7D%3D%5B1%2C...%2C1%5D%27) is a *n*-vector of ones, ![](https://latex.codecogs.com/gif.latex?%5Ctextbf%7BX%7D_k%3D%5C%7Bx_%7Bijk%7D%5C%7D) is the ![](https://latex.codecogs.com/gif.latex?n%5Ctimes%20p) matrix of centered and standardized markers available in the *k*<sup>th</sup> environment, ![](https://latex.codecogs.com/gif.latex?%5Cboldsymbol%7B%5Cbeta%7D_k%3D%5B%5Cbeta_%7B1k%7D%2C...%2C%5Cbeta_%7Bpk%7D%5D%27) is a *p*-vector of marker effects, and ![](https://latex.codecogs.com/gif.latex?%5Cboldsymbol%7B%5Cvarepsilon%7D_k%3D%5B%5Cvarepsilon_%7Bik%7D%2C...%2C%5Cvarepsilon_%7Bnk%7D%5D%27) is the *n*-vector of residuals.

This model can be represented as a G-BLUP model by setting ![](https://latex.codecogs.com/gif.latex?%5Ctextbf%7Bg%7D_k%3D%5Ctextbf%7BX%7D_k%5Cboldsymbol%7B%5Cbeta%7D_k), this is: 

<p align="center">
<img src="https://latex.codecogs.com/gif.latex?%5Ctextbf%7By%7D_k%3D%5Cmu_k%5Ctextbf%7B1%7D&plus;%5Ctextbf%7Bg%7D_k&plus;%5Cboldsymbol%7B%5Cvarepsilon%7D_k">
</p>

where ![](https://latex.codecogs.com/gif.latex?%5Ctextbf%7Bg%7D_k%3D%5Bg_%7B1k%7D%2C...%2Cg_%7Bnk%7D%5D%27) is a random effect assumed ![](https://latex.codecogs.com/gif.latex?%5Ctextbf%7Bg%7D_k%5Csim%20N%28%5Ctextbf%7B0%7D%2C%5Csigma%5E2_%7Bg_k%7D%5Ctextbf%7BG%7D_k%29) with the genomic relationship matrix estimated as ![](https://latex.codecogs.com/gif.latex?%5Ctextbf%7BG%7D_k%3D%5Ctextbf%7BX%7D_k%5Ctextbf%7BX%7D_k%27/p).

* **Across-environments model**

This model assumes that effects of markers are the same across environments, this is ![](https://latex.codecogs.com/gif.latex?%5Cboldsymbol%7B%5Cbeta%7D_1%3D%5Cboldsymbol%7B%5Cbeta%7D_2%3D...%3D%5Cboldsymbol%7B%5Cbeta%7D_s%3D%5Cboldsymbol%7B%5Cbeta%7D). The model above can be simultaneously fitted for all environments as (assume *s=3* environments)

<p align="center">
<img src="https://latex.codecogs.com/gif.latex?%5Cbegin%7Bbmatrix%7D%20%5Ctextbf%7By%7D_1%5C%5C%20%5Ctextbf%7By%7D_2%5C%5C%20%5Ctextbf%7By%7D_3%20%5Cend%7Bbmatrix%7D%3D%5Cbegin%7Bbmatrix%7D%20%5Cmu_1%5Ctextbf%7B1%7D%5C%5C%20%5Cmu_2%5Ctextbf%7B1%7D%5C%5C%20%5Cmu_3%5Ctextbf%7B1%7D%20%5Cend%7Bbmatrix%7D&plus;%5Cbegin%7Bbmatrix%7D%20%5Ctextbf%7BX%7D_1%5C%5C%20%5Ctextbf%7BX%7D_2%5C%5C%20%5Ctextbf%7BX%7D_3%20%5Cend%7Bbmatrix%7D%5Cboldsymbol%7B%5Cbeta%7D&plus;%20%5Cbegin%7Bbmatrix%7D%20%5Cboldsymbol%7B%5Cvarepsilon%7D_1%5C%5C%20%5Cboldsymbol%7B%5Cvarepsilon%7D_2%5C%5C%20%5Cboldsymbol%7B%5Cvarepsilon%7D_3%20%5Cend%7Bbmatrix%7D">
</p>

Similarly, it can be represented also as a G-BLUP model, by making ![](https://latex.codecogs.com/gif.latex?%5Ctextbf%7Bg%7D_k%3D%5Ctextbf%7BX%7D_k%5Cboldsymbol%7B%5Cbeta%7D), as: 

<p align="center">
<img src="https://latex.codecogs.com/gif.latex?%5Cbegin%7Bbmatrix%7D%20%5Ctextbf%7By%7D_1%5C%5C%20%5Ctextbf%7By%7D_2%5C%5C%20%5Ctextbf%7By%7D_3%20%5Cend%7Bbmatrix%7D%3D%5Cbegin%7Bbmatrix%7D%20%5Cmu_1%5Ctextbf%7B1%7D%5C%5C%20%5Cmu_2%5Ctextbf%7B1%7D%5C%5C%20%5Cmu_3%5Ctextbf%7B1%7D%20%5Cend%7Bbmatrix%7D&plus;%5Cbegin%7Bbmatrix%7D%20%5Ctextbf%7Bg%7D_1%5C%5C%20%5Ctextbf%7Bg%7D_2%5C%5C%20%5Ctextbf%7Bg%7D_3%20%5Cend%7Bbmatrix%7D&plus;%20%5Cbegin%7Bbmatrix%7D%20%5Cboldsymbol%7B%5Cvarepsilon%7D_1%5C%5C%20%5Cboldsymbol%7B%5Cvarepsilon%7D_2%5C%5C%20%5Cboldsymbol%7B%5Cvarepsilon%7D_3%20%5Cend%7Bbmatrix%7D">
</p>

where the random effect ![](https://latex.codecogs.com/gif.latex?%5Ctextbf%7Bg%7D_0%3D%5B%5Ctextbf%7Bg%7D%27_1%2C%5Ctextbf%7Bg%7D%27_2%2C%5Ctextbf%7Bg%7D%27_3%5D%27) is assumed ![](https://latex.codecogs.com/gif.latex?%5Ctextbf%7Bg%7D_0%5Csim%20N%28%5Ctextbf%7B0%7D%2C%5Csigma%5E2_%7Bg_0%7D%5Ctextbf%7BG%7D_0%29) with ![](https://latex.codecogs.com/gif.latex?%5Ctextbf%7BG%7D_0) being the marker-derived genomic relationship calculated as

<p align="center">
<img src="https://latex.codecogs.com/gif.latex?%5Ctextbf%7BG%7D_0%3D%5Cfrac%7B1%7D%7Bp%7D%20%5Cbegin%7Bbmatrix%7D%20%5Ctextbf%7BX%7D_1%5Ctextbf%7BX%7D_1%27%20%26%20%5Ctextbf%7BX%7D_1%5Ctextbf%7BX%7D_2%27%20%26%20%5Ctextbf%7BX%7D_1%5Ctextbf%7BX%7D_3%27%5C%5C%20%5Ctextbf%7BX%7D_2%5Ctextbf%7BX%7D_1%27%20%26%20%5Ctextbf%7BX%7D_2%5Ctextbf%7BX%7D_2%27%20%26%20%5Ctextbf%7BX%7D_2%5Ctextbf%7BX%7D_3%27%5C%5C%20%5Ctextbf%7BX%7D_3%5Ctextbf%7BX%7D_1%27%20%26%20%5Ctextbf%7BX%7D_3%5Ctextbf%7BX%7D_2%27%20%26%20%5Ctextbf%7BX%7D_3%5Ctextbf%7BX%7D_3%27%20%5Cend%7Bbmatrix%7D">
</p>

* **MxE model**

It models GxE interaction using a marker-by-environment (MxE) approach in which the effect of the *j*<sup>th</sup> marker on the *k*<sup>th</sup> environment, ![](https://latex.codecogs.com/gif.latex?%5Cbeta_%7Bjk%7D), is descomposed into an effect that is common to all environments (*b*<sub>j0</sub>) and an effect that is specific to each environment (*b*<sub>jk</sub>), this is ![](https://latex.codecogs.com/gif.latex?%5Cbeta_%7Bjk%7D%3Db_%7Bj0%7D&plus;b_%7Bjk%7D). Thus, the multi-environmental model is

<p align="center">
<img src="https://latex.codecogs.com/gif.latex?%5Cbegin%7Bbmatrix%7D%20%5Ctextbf%7By%7D_1%5C%5C%20%5Ctextbf%7By%7D_2%5C%5C%20%5Ctextbf%7By%7D_3%20%5Cend%7Bbmatrix%7D%3D%5Cbegin%7Bbmatrix%7D%20%5Cmu_1%5Ctextbf%7B1%7D%5C%5C%20%5Cmu_2%5Ctextbf%7B1%7D%5C%5C%20%5Cmu_3%5Ctextbf%7B1%7D%20%5Cend%7Bbmatrix%7D&plus;%5Cbegin%7Bbmatrix%7D%20%5Ctextbf%7BX%7D_1%5C%5C%20%5Ctextbf%7BX%7D_2%5C%5C%20%5Ctextbf%7BX%7D_3%20%5Cend%7Bbmatrix%7D%5Ctextbf%7Bb%7D_0&plus;%20%5Cbegin%7Bbmatrix%7D%20%5Ctextbf%7BX%7D_1%20%26%20%5Ctextbf%7B0%7D%20%26%20%5Ctextbf%7B0%7D%5C%5C%20%5Ctextbf%7B0%7D%20%26%20%5Ctextbf%7BX%7D_2%20%26%20%5Ctextbf%7B0%7D%5C%5C%20%5Ctextbf%7B0%7D%20%26%20%5Ctextbf%7B0%7D%20%26%20%5Ctextbf%7BX%7D_3%20%5Cend%7Bbmatrix%7D%20%5Cbegin%7Bbmatrix%7D%20%5Ctextbf%7Bb%7D_1%5C%5C%20%5Ctextbf%7Bb%7D_2%5C%5C%20%5Ctextbf%7Bb%7D_3%20%5Cend%7Bbmatrix%7D%20&plus;%5Cbegin%7Bbmatrix%7D%20%5Cboldsymbol%7B%5Cvarepsilon%7D_1%5C%5C%20%5Cboldsymbol%7B%5Cvarepsilon%7D_2%5C%5C%20%5Cboldsymbol%7B%5Cvarepsilon%7D_3%20%5Cend%7Bbmatrix%7D">
</p>

Likewise, by making ![](https://latex.codecogs.com/gif.latex?%5Ctextbf%7Bg%7D_%7B0,k%7D%3D%5Ctextbf%7BX%7D_k%5Ctextbf%7Bb%7D_0) and ![](https://latex.codecogs.com/gif.latex?%5Ctextbf%7Bg%7D_%7B1,k%7D%3D%5Ctextbf%7BX%7D_k%5Ctextbf%7Bb%7D_k), the G-BLUP representation is

<p align="center">
<img src="https://latex.codecogs.com/gif.latex?%5Cbegin%7Bbmatrix%7D%20%5Ctextbf%7By%7D_1%5C%5C%20%5Ctextbf%7By%7D_2%5C%5C%20%5Ctextbf%7By%7D_3%20%5Cend%7Bbmatrix%7D%3D%5Cbegin%7Bbmatrix%7D%20%5Cmu_1%5Ctextbf%7B1%7D%5C%5C%20%5Cmu_2%5Ctextbf%7B1%7D%5C%5C%20%5Cmu_3%5Ctextbf%7B1%7D%20%5Cend%7Bbmatrix%7D&plus;%20%5Cbegin%7Bbmatrix%7D%20%5Ctextbf%7Bg%7D_%7B0%2C1%7D%5C%5C%20%5Ctextbf%7Bg%7D_%7B0%2C2%7D%5C%5C%20%5Ctextbf%7Bg%7D_%7B0%2C3%7D%20%5Cend%7Bbmatrix%7D&plus;%20%5Cbegin%7Bbmatrix%7D%20%5Ctextbf%7Bg%7D_%7B1%2C1%7D%5C%5C%20%5Ctextbf%7Bg%7D_%7B1%2C2%7D%5C%5C%20%5Ctextbf%7Bg%7D_%7B1%2C3%7D%20%5Cend%7Bbmatrix%7D&plus;%20%5Cbegin%7Bbmatrix%7D%20%5Cboldsymbol%7B%5Cvarepsilon%7D_1%5C%5C%20%5Cboldsymbol%7B%5Cvarepsilon%7D_2%5C%5C%20%5Cboldsymbol%7B%5Cvarepsilon%7D_3%20%5Cend%7Bbmatrix%7D">
</p>

where the random effects 
<p align="center">
<img src="https://latex.codecogs.com/gif.latex?%5Ctextbf%7Bg%7D_0%3D%20%5Cbegin%7Bbmatrix%7D%20%5Ctextbf%7Bg%7D_%7B0%2C1%7D%5C%5C%20%5Ctextbf%7Bg%7D_%7B0%2C2%7D%5C%5C%20%5Ctextbf%7Bg%7D_%7B0%2C3%7D%20%5Cend%7Bbmatrix%7D%20%5Cqquad%5Ctext%7Band%7D%5Cqquad%20%5Ctextbf%7Bg%7D_1%3D%20%5Cbegin%7Bbmatrix%7D%20%5Ctextbf%7Bg%7D_%7B1%2C1%7D%5C%5C%20%5Ctextbf%7Bg%7D_%7B1%2C2%7D%5C%5C%20%5Ctextbf%7Bg%7D_%7B1%2C3%7D%20%5Cend%7Bbmatrix%7D">
</b>

are assumed ![](https://latex.codecogs.com/gif.latex?%5Ctextbf%7Bg%7D_0%5Csim%20N%28%5Ctextbf%7B0%7D%2C%5Csigma%5E2_%7Bg_0%7D%5Ctextbf%7BG%7D_0%29) and ![](https://latex.codecogs.com/gif.latex?%5Ctextbf%7Bg%7D_1%5Csim%20N%28%5Ctextbf%7B0%7D%2C%5Ctextbf%7BG%7D_1%29), respectively. Here, ![](https://latex.codecogs.com/gif.latex?%5Ctextbf%7BG%7D_0) is described previously and

<p align="center">
<img src="https://latex.codecogs.com/gif.latex?%5Ctextbf%7BG%7D_1%3D%5Cfrac%7B1%7D%7Bp%7D%20%5Cbegin%7Bbmatrix%7D%20%5Csigma%5E2_%7Bg_1%7D%5Ctextbf%7BX%7D_1%5Ctextbf%7BX%7D_1%27%20%26%20%5Ctextbf%7B0%7D%20%26%20%5Ctextbf%7B0%7D%5C%5C%20%5Ctextbf%7B0%7D%20%26%20%5Csigma%5E2_%7Bg_2%7D%5Ctextbf%7BX%7D_2%5Ctextbf%7BX%7D_2%27%20%26%20%5Ctextbf%7B0%7D%5C%5C%20%5Ctextbf%7B0%7D%20%26%20%5Ctextbf%7B0%7D%20%26%20%5Csigma%5E2_%7Bg_3%7D%5Ctextbf%7BX%7D_3%5Ctextbf%7BX%7D_3%27%20%5Cend%7Bbmatrix%7D">
</b>


* **Reaction Norm model**

The G-BLUP model presented in the across-environment approach can be extended to incorporate GxE by introducing covariance structures as a funcion of the marker information. This model, assumes the environment (*E*<sub>k</sub>) as a random effect as ![](https://latex.codecogs.com/gif.latex?%5Ctextbf%7BE%7D%3D%5BE_1%2C...%2CE_s%5D%5Csim%20N%28%5Ctextbf%7B0%7D%2C%5Csigma%5E2_e%5Ctextbf%7BI%7D%29). These random effects are conected with individuals through the design matrix **Z**<sub>e</sub>. Thus,  the main effect of environments, **e** = **Z**<sub>e</sub>**E**, is

<img src="https://latex.codecogs.com/gif.latex?%5Ctextbf%7Be%7D%3D%20%5Cbegin%7Bbmatrix%7D%20%5Ctextbf%7Be%7D_1%5C%5C%20%5Ctextbf%7Be%7D_2%5C%5C%20%5Ctextbf%7Be%7D_3%20%5Cend%7Bbmatrix%7D%3D%20%5Ctextbf%7BZ%7D_e%5Ctextbf%7BE%7D%5Csim%20N%28%5Ctextbf%7B0%7D%2C%5Csigma%5E2_e%5Ctextbf%7BZ%7D_e%5Ctextbf%7BZ%7D_e%27%29">

where ![](https://latex.codecogs.com/gif.latex?%5Ctextbf%7Be%7D_k%3D%5Be_%7B1k%7D%2C...%2Ce_%7Bnk%7D%5D%27) are the effects of the environment *k* for the *i*<sup>th</sup> individual (*i=1,...,n*). The reaction norm model, incorporates GxE by introducing the interaction terms ![](https://latex.codecogs.com/gif.latex?g_%7Bik%7De_%7Bik%7D), as 

<p align="center">
<img src="https://latex.codecogs.com/gif.latex?%5Cbegin%7Bbmatrix%7D%20%5Ctextbf%7By%7D_1%5C%5C%20%5Ctextbf%7By%7D_2%5C%5C%20%5Ctextbf%7By%7D_3%20%5Cend%7Bbmatrix%7D%3D%20%5Cbegin%7Bbmatrix%7D%20%5Cmu_1%5Ctextbf%7B1%7D%5C%5C%20%5Cmu_2%5Ctextbf%7B1%7D%5C%5C%20%5Cmu_3%5Ctextbf%7B1%7D%20%5Cend%7Bbmatrix%7D&plus;%20%5Cbegin%7Bbmatrix%7D%20%5Ctextbf%7Bg%7D_1%5C%5C%20%5Ctextbf%7Bg%7D_2%5C%5C%20%5Ctextbf%7Bg%7D_3%20%5Cend%7Bbmatrix%7D&plus;%20%5Cbegin%7Bbmatrix%7D%20%5Ctextbf%7Bg%7D_1%5Ctextbf%7Be%7D_1%5C%5C%20%5Ctextbf%7Bg%7D_2%5Ctextbf%7Be%7D_2%5C%5C%20%5Ctextbf%7Bg%7D_3%5Ctextbf%7Be%7D_3%20%5Cend%7Bbmatrix%7D&plus;%20%5Cbegin%7Bbmatrix%7D%20%5Cboldsymbol%7B%5Cvarepsilon%7D_1%5C%5C%20%5Cboldsymbol%7B%5Cvarepsilon%7D_2%5C%5C%20%5Cboldsymbol%7B%5Cvarepsilon%7D_3%20%5Cend%7Bbmatrix%7D">
</p>

where ![](https://latex.codecogs.com/gif.latex?%5Ctextbf%7Bg%7D_k%5Ctextbf%7Be%7D_k) is the interaction term between the random effect of markers and environments. It can be proved that, aproximately, 

<p align="center">
<img src="https://latex.codecogs.com/gif.latex?%5Ctextbf%7Bg%7D%5Ctextbf%7Be%7D%3D%20%5Cbegin%7Bbmatrix%7D%20%5Ctextbf%7Bg%7D_1%5Ctextbf%7Be%7D_1%5C%5C%20%5Ctextbf%7Bg%7D_2%5Ctextbf%7Be%7D_2%5C%5C%20%5Ctextbf%7Bg%7D_3%5Ctextbf%7Be%7D_3%20%5Cend%7Bbmatrix%7D%20%5Csim%20N%28%5Ctextbf%7B0%7D%2C%5Csigma%5E2_%7Bge%7D%5Ctextbf%7BG%7D_0%5Ccirc%20%5B%5Ctextbf%7BZ%7D_e%5Ctextbf%7BZ%7D_e%27%5D%20%29">
</b>

where ![](https://latex.codecogs.com/gif.latex?%5Ctextbf%7BG%7D_0) was described above and ![](https://latex.codecogs.com/gif.latex?%5Ccirc) denotes the element by element product (known as Hadamard product).

**Note**. The Hadamard product above will yield a matrix similar to ![](https://latex.codecogs.com/gif.latex?%5Ctextbf%7BG%7D_1) above described but with ![](https://latex.codecogs.com/gif.latex?%5Csigma%5E2_%7Bg_1%7D%3D%5Csigma%5E2_%7Bg_2%7D%3D%5Csigma%5E2_%7Bg_3%7D%3D1). This is, reaction norm model estimates environment-specific effects with constant variance ![](https://latex.codecogs.com/gif.latex?%5Csigma%5E2_%7Bge%7D) in contrast wth MxE model that accomodates environment-specific variance, ![](https://latex.codecogs.com/gif.latex?%5Csigma%5E2_%7Bg_k%7D), *k=1,...,s*.

## Model assessment
Performance of the *Reaction Norm* and *MxE model* will be compared with that of the *across-environments* model that ignores GxE modeling and the *single-environment* model which is fitted within each environment.

### Training-Testing random partitions.
The prediction power of the model will be assessed using the training-testing (TRN-TST) random partitions approach. 
Data is randomly splitted into training and testing sets. Model parameters are estimated in training set and model is tested in TST set.  Two main estimations problems are addressed using the multi-environments models. 

*   **Cross Validation 1 (CV1)**. Represent a scheme of prediction of lines that have not been evaluated in any field
trials.

*   **Cross Validation 2 (CV2)**. Represent a scheme of prediction of lines that have been evaluated in some but all target environments. Thus, prediction of non-evaluated lines benefits from borrowing of information from lines that were evaluated in other environments.

<p align="center">
<img src="https://github.com/MarcooLopez/Genomic-Selection-Demo/blob/master/CV1_2_scheme.png" width="450">
</b>

In our case, we will use 70% of the data for training set and the remaining 30% for the testing set.
For CV1, we will create a scheme in which 30% of the lines are missing in all environments. 
CV2 scheme is created by having 30% of the entries missing in one environment but present in all the rest of environments.

This procedure of TRN-TST can be repeated many times to allow for estimation of standard errors (SE).

## Data preparation
### Load data, generate G-matrix
The following R code, [prepareData_multi.R](https://github.com/MarcooLopez/Genomic-Selection/blob/master/prepareData_multi.R), can be used to prepare the data for analizes of the multi-environmental models.

```
setwd("/mnt/home/lopezcru/GS")
rm(list=ls())
library(BGLR)

# Load data
data(wheat)
X <- wheat.X
Y <- wheat.Y

# Select environments. For instance, environments 2,4, and 5
Y <- Y[,c(2,3,4)]

# Genomic relationship matrix
M <- scale(X)
G <- tcrossprod(M)/ncol(X)

# Design matrix for individuals. It connects individuals with environments
GID <- factor(rep(rownames(Y),ncol(Y)),levels=rownames(Y))
Zg <- model.matrix(~GID-1)

# Design matrix for environments. Used in the multi-environment R-Norm model
envID <- factor(rep(colnames(Y),each=nrow(Y)),levels=colnames(Y))
ZE <- model.matrix(~envID-1)   

#  Covariance structure for effects
ZgGZgt <- Zg%*%G%*%t(Zg)    # Genetic effect  
ZEZEt <- tcrossprod(ZE)     # Environmental effect
GE <- ZgGZgt*ZEZEt          # GxE interaction term (R-Norm model)

# Eigen decomposition (to speed computational time)
eigen_G <- eigen(G)
eigen_G0 <- eigen(ZgGZgt)
eigen_GE <- eigen(GE)

# Interaction terms (MxE model)
MxE_eigen <- vector("list",ncol(Y))
for(env in 1:ncol(Y)){ 
    tmp <- rep(0,ncol(Y)) ; tmp[env] <- 1; G1 <- kronecker(diag(tmp),G)
    MxE_eigen[[env]] <- eigen(G1)
}

# Save prepared data
dir.create("multiEnvironment")
save(Y,envID,eigen_G,eigen_G0,eigen_GE,MxE_eigen,file="multiEnvironment/prepData_multi.RData")
```

## Running models

### 1. Variance components estimation 
Code below, [get_VarComps_multi.R](https://github.com/MarcooLopez/Genomic-Selection/blob/master/get_VarComps_multi.R) script, can be used after 'data preparation' part to fit all the models and to extract variance components.

```
setwd("/mnt/home/lopezcru/GS")
rm(list=ls())
library(BGLR)

load("multiEnvironment/prepData_multi.RData")
n <- nrow(Y)
nEnv <- ncol(Y)
y <- as.vector(Y)

set.seed(123)

# Matrix to store results. It will save variance components for each model
outVAR <- matrix(NA,ncol=4,nrow=1+2*nEnv)
dimnames(outVAR) <- list(c("Main",rep(paste0("Env ",colnames(Y)),2)),c("Single","Across","MxE","R-Norm"))

# Number of iterations and burn-in for Bayesian models
nIter <- 30000; burnIn <- 2000

#--------------------------------------------------------
# 1. Single environment (within-environment) model
#--------------------------------------------------------
ETA <- list(G=list(V=eigen_G$vectors,d=eigen_G$values,model='RKHS'))
for(env in 1:nEnv){
    fm <-BGLR(y=Y[,env],ETA=ETA,nIter=nIter,burnIn=burnIn)
    outVAR[env+1,1] <- fm$ETA[[1]]$varU
    outVAR[env+4,1] <- fm$varE
}

#--------------------------------------------------------
# 2. Across-environments model
#--------------------------------------------------------
ETA <- list(list(~envID-1,model="FIXED"))
ETA[[2]] <- list(V=eigen_G0$vectors,d=eigen_G0$values,model='RKHS')

# Model Fitting
fm <- BGLR(y=y,ETA=ETA,nIter=nIter,burnIn=burnIn)
outVAR[1,2] <- fm$ETA[[2]]$varU
outVAR[(1:nEnv)+4,2] <- fm$varE

#--------------------------------------------------------
# 3. MxE interaction model
#--------------------------------------------------------
ETA <- list(list(~envID-1,model="FIXED"))
ETA[[2]] <- list(V=eigen_G0$vectors,d=eigen_G0$values,model='RKHS')

# Adding interaction terms
for(env in 1:nEnv){
    eigen_G1 <- MxE_eigen[[env]]
    ETA[[(env+2)]] <- list(V=eigen_G1$vectors,d=eigen_G1$values,model='RKHS')
}

# Model Fitting
fm <- BGLR(y=y,ETA=ETA,nIter=nIter,burnIn=burnIn)
outVAR[1,3] <- fm$ETA[[2]]$varU
for(env in 1:nEnv) outVAR[env+1,3] <- fm$ETA[[env+2]]$varU
outVAR[(1:nEnv)+4,3] <- fm$varE

#--------------------------------------------------------
# 4. Reaction-Norm model
#--------------------------------------------------------
ETA <- list(list(~envID-1,model="FIXED"))
ETA[[2]] <- list(V=eigen_G0$vectors,d=eigen_G0$values,model='RKHS')
ETA[[3]] <- list(V=eigen_GE$vectors,d=eigen_GE$values,model="RKHS")

# Model Fitting
fm <- BGLR(y=y,ETA=ETA,nIter=nIter,burnIn=burnIn)
outVAR[1,4] <- fm$ETA[[2]]$varU
outVAR[(1:nEnv)+1,4] <- fm$ETA[[3]]$varU
outVAR[(1:nEnv)+4,4] <- fm$varE
outVAR

# Save results
write.table(outVAR,file="multiEnvironment/varComps.csv",sep=",",row.names=F)
```

#### Results
The following table is the output of the code above for `nIter=30000` and `burnIn=2000`.
<img src="https://github.com/MarcooLopez/Genomic-Selection/blob/master/varComp.png" width="360">

##
### 2. Replicates of partitions to obtain standard deviations of predictions
Using a GBLUP approach, the prediction power of the multi-environment models (MxE and Reaction Norm) will be compared with that that ignores GxE (across-environment) and with the GBLUP model fitted within environment. 

#### 2.1. Training-Testing partitions
After running the 'data preparation' part, it can be chosen either to perform CV1 or CV2 aproaches

* **Cross Validation 1 (CV1)**

Code below will generate a matrix YNA containing "NA" values for the entries corresponding to the TST set mimicing the CV1 prediction problem. It generates a 'list' with 'm' matrices containing the TRN-TST partitions

```
setwd("/mnt/home/lopezcru/GS")
rm(list=ls())

#=========================================================
# User specifications
#=========================================================
# Number of replicates
m <- 100

# Percentage of the data assigned to Testing set
percTST <- 0.3
#=========================================================

# Load data
load("multiEnvironment/prepData_multi.RData")
n <- nrow(Y)

# Creation of seed for repeated randomizations
set.seed(123)
seeds <- round(seq(1E3,1E6,length=m))

nTST <- round(percTST*n)
YNA <- vector("list",m)

for(k in 1:m)
{
    set.seed(seeds[k])
    indexTST <- sample(1:n,size=nTST,replace=FALSE)
    YNA0 <- Y
    YNA0[indexTST,] <- NA
    YNA[[k]] <- YNA0
}

# Save YNA matrix
save(YNA,file="multiEnvironment/YNA_CV1_multiEnv.RData")
```

* **Cross Validation 2 (CV2)**

Code below will generate a matrix YNA containing "NA" values for the entries corresponding to the TST set mimicing the CV2 prediction problem. It generates a 'list' with 'm' matrices containing the TRN-TST partitions

```
setwd("/mnt/home/lopezcru/GS")
rm(list=ls())

#=========================================================
# User specifications
#=========================================================
# Number of replicates
m <- 100

# Percentage of the data assigned to Testing set
percTST <- 0.3
#=========================================================

# Load data
load("multiEnvironment/prepData_multi.RData")
n <- nrow(Y)
nEnv <- ncol(Y)

# Creation of seed for repeated randomizations
set.seed(123)
seeds <- round(seq(1E3,1E6,length=m))

nTST <- round(percTST*n)
YNA <- vector("list",m)
nNA <- nEnv*nTST

for(k in 1:m)
{
    set.seed(seeds[k])
    YNA0 <- Y
    
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
    for(j in 1:nEnv) YNA0[indexNA[indexEnv==j],j] <- NA
    YNA[[k]] <- YNA0
}

# Save YNA matrix
save(YNA,file="multiEnvironment/YNA_CV2_multiEnv.RData")
```

After running the code to generate partitions for either CV1 or CV2 scenarios, the following script ([fitModels_multi.R](https://github.com/MarcooLopez/Genomic-Selection/blob/master/fitModels_multi.R)) can be run to fit the models repeatealy for all partitions. In all multi-environment models, main effect of 'environment' will be regarded as fixed effect.

The code runs a single partition for each model either for CV1 or CV2. These specifications need to be passed in variables `mod`, `CV`, and `part`. 

```
setwd("/mnt/home/lopezcru/GS")
rm(list=ls())
library(BGLR)

#=========================================================
# User specifications
#=========================================================
# Choose one model. 1: single; 2:across; 3:MxE; 4:R-Norm
mod <- 4

# Type of CV. 1:CV1; 2:CV2
CV <- 1

# Partition number
part <- 1
#=========================================================

# Read arguments passed from command line
args=(commandArgs(TRUE))
if(length(args)==0){
   cat('No args provided',"\n")
 }else{
    for(i in 1:length(args)){
         eval(parse(text=args[[i]]))
    }
}

# Load data
load("multiEnvironment/prepData_multi.RData")
load(paste0("multiEnvironment/YNA_CV",CV,"_multiEnv.RData"))
n <- nrow(Y);  nEnv <- ncol(Y)

# Models
models <- c("Single","Across","MxE","R-Norm")
model <- models[mod]

# Number of iterations and burn-in for Bayesian models
nIter <- 30000;  burnIn <- 2000

YNA0 <- YNA[[part]]
yNA <- as.vector(YNA0)
    
#--------------------------------------------------------
# 1. Single environment (within-environment) model
#--------------------------------------------------------
if(model=="Single")
{
    YHat <- matrix(NA,nrow=nrow(Y),ncol=ncol(Y))
    ETA <- list(G=list(V=eigen_G$vectors,d=eigen_G$values,model='RKHS'))
    for(env in 1:nEnv){
        fm <-BGLR(y=YNA0[,env],ETA=ETA,nIter=nIter,burnIn=burnIn)
        YHat[,env] <- fm$yHat
    }
}

#--------------------------------------------------------
# 2. Across-environments model
#--------------------------------------------------------
if(model=="Across")
{
    ETA <- list(list(~envID-1,model="FIXED"))
    ETA[[2]] <- list(V=eigen_G0$vectors,d=eigen_G0$values,model='RKHS')

    # Model Fitting
    fm <- BGLR(y=yNA,ETA=ETA,nIter=nIter,burnIn=burnIn)
    YHat <- matrix(fm$yHat,ncol=nEnv)
}
    
#--------------------------------------------------------
# 3. MxE interaction model
#--------------------------------------------------------
if(model=="MxE")
{
    ETA <- list(list(~envID-1,model="FIXED"))
    ETA[[2]] <- list(V=eigen_G0$vectors,d=eigen_G0$values,model='RKHS')

    # Adding interaction terms
    for(env in 1:nEnv){
        eigen_G1 <- MxE_eigen[[env]]
        ETA[[(env+2)]] <- list(V=eigen_G1$vectors,d=eigen_G1$values,model='RKHS')
    }

    # Model Fitting
    fm <- BGLR(y=yNA,ETA=ETA,nIter=nIter,burnIn=burnIn)
    YHat <- matrix(fm$yHat,ncol=nEnv)
}

#--------------------------------------------------------
# 4. Reaction-Norm model
#--------------------------------------------------------
if(model=="R-Norm")
{
    ETA <- list(list(~envID-1,model="FIXED"))
    ETA[[2]] <- list(V=eigen_G0$vectors,d=eigen_G0$values,model='RKHS')
    ETA[[3]] <- list(V=eigen_GE$vectors,d=eigen_GE$values,model="RKHS")

    # Model Fitting
    fm <- BGLR(y=yNA,ETA=ETA,nIter=nIter,burnIn=burnIn)
    YHat <- matrix(fm$yHat,ncol=nEnv)
}

# Save results
outfolder <- paste0("multiEnvironment/CV",CV,"/",model)
if(!file.exists(outfolder)) dir.create(outfolder,recursive=T)
save(YHat,file=paste0(outfolder,"/outPRED_multiEnv_partition_",part,".RData"))
```

#### 2.2 Running in parallel many jobs
Code above will run a single combination of partition-model-CV, thus when running, for instance, several models for both CV1 and CV2, some parallelzation of jobs is needed for speeding of computation. The bash code called *[run_jobs_multi.sh](https://github.com/MarcooLopez/Genomic-Selection/blob/master/run_jobs_multi.sh)* will submit many jobs depending of the core capacity of the computer. Jobs will be sent by chunks whose size is specified in variable `nb` (for instance, `nb=10` will run 10 jobs at the time). After jobs in the chunk are done, another chunk will be submited to be run. Variables `seq1=2`, `seq2=4`, and `seq3=100` specify the 2 CV types, 4 models, and 100 partitions, respectiely.

The bash shell script needs to be saved in the same directory as the R script 'fitModels_multi.R' and it can be run from command line as
```
sh run_jobs_multi.sh &
```

#### 2.3 Retrieving results

The code below will retrieve results for all models fitted previously showing the within-environment correlation for all fitted models

```
setwd("/mnt/home/lopezcru/GS")
rm(list=ls())
library(ggplot2)
library(reshape)

#=========================================================
# User specifications
#=========================================================
# Type of CV. 1:CV1; 2:CV2
CV <- 2
#=========================================================

models <- c("Single","Across","MxE","R-Norm")

# Load data
load("multiEnvironment/prepData_multi.RData")
load(paste0("multiEnvironment/YNA_CV",CV,"_multiEnv.RData"))

# Calculate within-environment correlation
outCOR <- vector("list",length(models))
names(outCOR) <- models
for(mod in seq_along(models))
{
    outcor <- c()
    for(part in 1:length(YNA)){
        filename <- paste0("multiEnvironment/CV",CV,"/",models[mod],"/outPRED_multiEnv_partition_",part,".RData")
        if(file.exists(filename))
        {
            load(filename)
            YNA0 <- YNA[[part]]
            tmp <- rep(NA,ncol(YNA0))
            for(env in 1:ncol(YNA0)){
                indexTST <- which(is.na(YNA0[,env]))
                tmp[env] <- cor(Y[indexTST,env],YHat[indexTST,env])
            }
            outcor <- rbind(outcor,tmp)
        }
    }    
    colnames(outcor) <- paste0("Env ",colnames(YNA0))
    rownames(outcor) <- NULL
    outcor <- data.frame(model=models[mod],outcor)
    outCOR[[mod]] <- outcor
}
outCOR <- outCOR[!sapply(outCOR,is.null)]

# Calculate means and SD's
(means <- t(do.call("rbind",lapply(outCOR,function(x)apply(x[,-1],2,mean)))))
(sds <- t(do.call("rbind",lapply(outCOR,function(x)apply(x[,-1],2,sd)))))

write.csv(rbind(means,colnames(sds),sds),file=paste0("multiEnvironment/Accuracy_avg_CV",CV,"_multiEnv.csv"))

toplot <- do.call("rbind",lapply(outCOR,function(x)melt(x,id="model")))
png(paste0("multiEnvironment/Accuracy_distn_CV",CV,"_multiEnv.png"),height=350)
ggplot(toplot,aes(x=model,y=value,fill=variable)) + geom_boxplot()+
labs(fill="Env",y="Accuracy",title=paste0("Correlation between observed and predicted values. CV",CV))
dev.off()
```

Tables below are the results of running 100 partitions with `nIter=30000` and `burnIn=2000`. The mean and standard deviation (in parenthesis) across partitions are presented.

**Cross Validation 1. CV1**
 
|       |Single-Env |Across-Env | MxE  | RNorm |
|-------|-------|--------|------|------|
|Env 2  | 0.485(0.049)  | 0.441(0.052)  | 0.460(0.050) | 0.461(0.049) |
|Env 4  | 0.377(0.055)  | 0.395(0.053)  | 0.382(0.055) | 0.382(0.055) |
|Env 5  | 0.441(0.056)  | 0.382(0.057)  | 0.412(0.054) | 0.409(0.055) |

&nbsp;

<img src="https://github.com/MarcooLopez/Genomic-Selection/blob/master/Accuracy_distn_CV1_multiEnv.png" width="400">

###
**Cross Validation 2. CV2**

|       |Single-Env |Across-Env | MxE  | RNorm |
|-------|-------|--------|------|-----|
|Env 2  | 0.485(0.049)  | 0.629(0.032)  | 0.647(0.033) | 0.642(0.033) |
|Env 4  | 0.375(0.058)  | 0.602(0.042)  | 0.591(0.043) | 0.585(0.044) |
|Env 5  | 0.442(0.048)  | 0.493(0.045)  | 0.529(0.042) | 0.528(0.043) |

&nbsp;

<img src="https://github.com/MarcooLopez/Genomic-Selection/blob/master/Accuracy_distn_CV2_multiEnv.png" width="400">

#
* **[back](https://github.com/MarcooLopez/Genomic-Selection-Demo/blob/master/README.md)**

#
# References
* Jarquín, D., Crossa, J., Lacaze, X., Du Cheyron, P., Daucourt, J., Lorgeou, J., … de los Campos, G. (2014). **A reaction norm model for genomic selection using high-dimensional genomic and environmental data**. Theoretical and Applied Genetics, 127(3), 595–607. 
* Lopez-Cruz, M., Crossa, J., Bonnett, D., Dreisigacker, S., Poland, J., Jannink, J.-L., … de los Campos, G. (2015). **Increased Prediction Accuracy in Wheat Breeding Trials Using a Marker × Environment Interaction Genomic Selection Model**. G3: Genes, Genomes, Genetics, 5(4), 569–582. 



