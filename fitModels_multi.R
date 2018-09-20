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
