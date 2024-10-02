########################################################################################
# Simulation-based power analysis
# LIVE BETTER Trial Protein biomarker study
# Aim 1 - linear regression model w/ interaction terms
# Date: 2024.09.2x
# Technical notes: 
# 1. Main-interaction hierarchy structure will be followed. 
#   for penalization method: 
#         protein and its interaction will form a group
#         group selection will be used
#   for marginal method:
#         main effect of the interaction term < 0.05
#         interaction term will be selected if < 0.1 after FDR correction
########################################################################################


# clean objects
rm(list=ls())

# input parameters
ntp = 20          # number of true positives
nrep = 100        # number of replicates

# dimension from LIVE BETTER Trial
n = 640           # sample size
n1 = n2 = n/2     # equal allocation to treatment
p = 3072          # number of proteins

########################################################################################
# The following codes about specifications for generating protein data are copied from
# Bruce S. sbruce23/REHAB-HFpEF-simulated-power-analysis: Initial release (1.0). 
# Zenodo. 2023. doi: 10.5281/zenodo.8388571. 
# Bruce obtained those specifications by examining the pilot data.
########################################################################################
# generate multivariate normal to be used to simulate 3,072 protein expression levels 
mus=rgamma(p,shape=6.3,rate=1.2)
vars=rlnorm(p,meanlog=-0.87,sdlog=0.76)
covs=rnorm((p^2-p)/2,mean=0.095,sd=0.10)

# find nearest positive definite covariance matrix
sigma=diag(vars)
sigma[lower.tri(sigma)]=covs
sigma[upper.tri(sigma)]=covs
sigma = (sigma + t(sigma))/2
sigma.evd=eigen(sigma)
########################################################################################
eS          <- eigen(sigma)
ev          <- eS$values

########################################################################################
# simulate data and fit gLASSO for each replicate
########################################################################################
# - trt variable will not be penalized
# - group structure will be incorperated 
# - tuning is selected by $lambda.min 
#   other option is $lambda.1se
# - all other settings are as default.
########################################################################################

# set up different coef

coefList     <- list()

coefList[[1]]<- list(
  main=c(runif(n = (ntp/2),0.4,0.6),
         runif(n = (ntp/2),min=-0.6,max=-0.4)),
  init=c(runif(n =(ntp/2),0.375, 0.6))
)

coefList[[2]]<- list(
  main=c(runif(n = (ntp/2),0.4,0.6),
         runif(n = (ntp/2),min=-0.6,max=-0.4)),
  init=c(runif(n =(ntp/2),0.575, 0.8))
)

coefList[[3]]<- list(
  main=c(runif(n = (ntp/2),0.4,0.6),
         runif(n = (ntp/2),min=-0.6,max=-0.4)),
  init=c(runif(n =(ntp/2),0.675, 0.95))
)


library(foreach)
library(doParallel)

cores <- detectCores(logical = T)
cl    <- makeCluster(cores)
registerDoParallel(cl, cores = cores)

# loop through replicates 

rList <- foreach(r = 1:nrep) %dopar%{
  
  library(grpreg)
  
  #  simulate protein expression levels from multivariate normal
  X1          <- matrix(rnorm(p * n), nrow=n)
  u1.matrix   <- matrix(rep(mus,times=n),ncol=p,byrow=TRUE)
  dt.prot     <- u1.matrix + X1 %*% diag(sqrt(pmax(ev, 0)), p) %*% t(eS$vectors)
  
  # set trt information
  trt = c(rep(1,n1),rep(0,n2) )
  dt.prot.trt = apply(dt.prot,2,function(x) x*trt)
  
  # construct the entire dataset
  dt = cbind(trt,dt.prot,dt.prot.trt)
  colnames(dt) = c("trt",paste("protein",1:p,sep="_"),
                   paste("proteinXtrt",1:p,sep="_"))
  
  # simulate responses for logistic model and calculate the estimates
  # define the function
  
  cal.fun <- function(coefinput){
    
    ## construct data
    p <- dim(dt.prot)[2]
    
    coefinput.main  <- coefinput[["main"]]
    coefinput.inter <- coefinput[["init"]]
    
    mean.xb   <- trt+dt.prot[,1:ntp] %*% matrix(coefinput.main,ncol=1) +
      dt.prot.trt[,1:(ntp/2)] %*% matrix(coefinput.inter,ncol=1)  
    
    y <- mean.xb + rnorm(n=n,mean = 0,sd = 1)
    y <- y/sd(y)*0.15
    
    ## fit lasso
    # use group LASSO to estimate model parameters
    # check   plot(fit.high) for the first few cases to make sure the 
    # turning is reasonably selected
    
    fit <- cv.grpreg(dt, y,
                     group = c("0",as.character(1:p),as.character(1:p)),
                     penalty = "grLasso",
                     family="gaussian")
    
    fit.sel <- names(coef(fit, lambda=fit$lambda.min))[which(coef(fit, lambda=fit$lambda.min)!=0)]
    
    
    fit.sel.main <- grep(pattern = "protein_",x = fit.sel,value = T)
    fit.selIndex.main <- as.numeric( gsub(pattern = "protein_",replacement = "",x = fit.sel.main))
    
    
    fit.sel.inter <- grep(pattern = "proteinXtrt_",x = fit.sel,value = T)
    fit.selIndex.inter <- as.numeric( gsub(pattern = "proteinXtrt_",replacement = "",x = fit.sel.inter))
    
    fit.selIndex <- c(fit.selIndex.inter,fit.selIndex.main)
    
    fit.selIndex.both <- fit.selIndex[duplicated(fit.selIndex)]
    
    ntp.pen <- sum(fit.selIndex.both<=(ntp/2))
    nfp.pen <- length(fit.selIndex.both) - ntp.pen
    
    ## fit marginal model
    margin.p.main <- sapply(X = 1:p,
                            FUN = function(i){summary(
                              lm(y~ trt+dt.prot[,i]+dt.prot.trt[,i]))$coef[3,"Pr(>|t|)"]})
    margin.p.inter <- sapply(X = 1:p,
                             FUN = function(i){summary(
                               lm(y~ trt+dt.prot[,i]+dt.prot.trt[,i]))$coef[4,"Pr(>|t|)"]})
    
    fit.selIndex.main  <- which(margin.p.main<0.05)
    fit.selIndex.inter <- which(p.adjust(margin.p.inter, method = "fdr")<0.1)
    fit.selIndex      <- c(fit.selIndex.inter,fit.selIndex.main)
    fit.selIndex.both <- fit.selIndex[duplicated(fit.selIndex)]
    
    ntp.mar = sum(fit.selIndex.both<=(ntp/2))
    nfp.mar = length(fit.selIndex.both)- ntp.mar
    
    return(data.frame(
      sdy = sd(y),
      ntp.pen = ntp.pen,
      nfp.pen = nfp.pen,
      ntp.mar = ntp.mar,
      nfp.mar = nfp.mar,
      rtp.pen = ntp.pen/(ntp/2),
      rfp.pen = nfp.pen/(p-ntp/2),
      rtp.mar = ntp.mar/(ntp/2),
      rfp.mar = nfp.mar/(p-ntp/2) 
    )
    )
    
  }
  
  lapply(coefList,FUN = cal.fun)
  
}

stopImplicitCluster()
stopCluster(cl)


# process simulation results

rltAvg <- sapply(1:length(coefList), function(j) {
  # Get the variable names from the first element of the list for a given j
  var_names <- names(rList[[1]][[j]])
  
  # For each variable, sum its values across i
  tmp <- lapply(var_names, function(var) {
    mean(sapply(1:nrep, function(i) rList[[i]][[j]][[var]]))
  })
  
  names(tmp) <- var_names
  
  tmp
})


rltSD <- sapply(1:length(coefList), function(j) {
  var_names <- names(rList[[1]][[j]])
  tmp <- lapply(var_names, function(var) {
    sd(sapply(1:nrep, function(i) rList[[i]][[j]][[var]]))
  })
  
  names(tmp) <- paste0(var_names, "SD")
  tmp
})


saveRDS(object = list(coefList=coefList,rList=rList,
                      rltAvg=rltAvg, rltSD=rltSD),
        file = paste("PowerAnalysis-aim1-EQ5D5L-ntp",ntp,".RDS",sep=""))