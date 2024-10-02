########################################################################################
# Simulation based power analysis
# LIVE BETTER Trial Protein biomarker study
# Aim 3 - logistic regression model
# Date: 2024.09.2x
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
# ########################################################################################
eS          <- eigen(sigma)
ev          <- eS$values


########################################################################################
# simulate data and fit LASSO for each replicate
########################################################################################
# - trt variable will not be penalized
# - tuning is selected by $lambda.min 
#   other option is $lambda.1se
# - 10-fold CV (default) is used
# - all other settings are as default.
########################################################################################

# set up different coef
coefList     <- list()

data_setting_ntp20_mat <- matrix(c(
  c(0.45, 0.58, -0.55, -0.45),
  c(0.82, 0.92, -0.85, -0.75),
  c(1.3,  1.45, -1.3,  -1.2)), nrow = 3, byrow = TRUE)

data_setting_ntp20 <- data.frame(data_setting_ntp20_mat)
colnames(data_setting_ntp20) <- paste0("range", 1:4)

for (i in 1:nrow(data_setting_ntp20)) {
  mtp <- i
  range1 <- as.numeric(data_setting_ntp20[i, c(1, 2)])
  range2 <- as.numeric(data_setting_ntp20[i, c(3, 4)])
  print(c(range1, range2))
  coefList[[i]]<- c(runif(n=(ntp/2), min=range1[1], max=range1[2]),
                    runif(n=(ntp/2), min=range2[1], max=range2[2]))
}

# set up for parallel computing
library(foreach)
library(doParallel)

cores <- detectCores(logical = T)
cl    <- makeCluster(cores)
registerDoParallel(cl, cores = cores)

# loop through replicates 

rList <- foreach(r = 1:nrep) %dopar%{
  
  library(glmnet)

  #  simulate protein expression levels from multivariate normal
  X1          <- matrix(rnorm(p * n), nrow=n)
  u1.matrix   <- matrix(rep(mus,times=n),ncol=p,byrow=TRUE)
  dt.prot     <- u1.matrix + X1 %*% diag(sqrt(pmax(ev, 0)), p) %*% t(eS$vectors)
  
  # set trt information
  trt = c(rep(1,n1),rep(0,n2) )
  
  # construct the entire dataset
  dt = cbind(trt,dt.prot)
  colnames(dt) = c("trt",paste("protein",1:3072,sep="_"))
  
  # simulate responses for logistic model and calculate the estimates
  cal.fun <- function(coefinput){
    ## construct data
    p <- dim(dt.prot)[2]
    mean.xb   <- trt+dt.prot[,1:ntp] %*% matrix(coefinput,ncol=1)
    p.logit   <- exp(mean.xb)/(1+exp(mean.xb))
    y.logit   <- rbinom(length(p.logit), size = 1, prob = p.logit)

    ## fit lasso 
    fit <- cv.glmnet(dt, y.logit,
                     penalty.factor = c(0,rep(1,p)),
                     penfamily=binomial())
    
    fit.sel <- rownames(coef(fit, lambda=fit$lambda.min))[which(coef(fit, lambda=fit$lambda.min)!=0)]
    fit.sel <- grep(pattern = "protein_",x = fit.sel,value = T)
    fit.selIndex <- as.numeric( gsub(pattern = "protein_",replacement = "",x = fit.sel))
    
    ntp.pen <- sum(fit.selIndex<=ntp)
    nfp.pen <- length(fit.selIndex) - ntp.pen
    
    ## fit marginal model
    margin.p <- sapply(X = 1:p,
                       FUN = function(i){summary(
                         glm(y.logit~ trt+dt.prot[,i],
                             family = "binomial"))$coef[3,"Pr(>|z|)"]})
    
    
    ntp.mar = sum(which(p.adjust(margin.p, method = "fdr")<0.1)<=ntp)
    nfp.mar = length(which(p.adjust(margin.p, method = "fdr")<0.1)) - ntp.mar
    
    return(data.frame(
      ntp.pen = ntp.pen,
      nfp.pen = nfp.pen,
      ntp.mar = ntp.mar,
      nfp.mar = nfp.mar,
      rtp.pen = ntp.pen/ntp,
      rfp.pen = nfp.pen/(p-ntp),
      rtp.mar = ntp.mar/ntp,
      rfp.mar = nfp.mar/(p-ntp))
    )
  }
  lapply(coefList,FUN = cal.fun)
}

stopImplicitCluster()
stopCluster(cl)


# process simulation results
rltAvg <- sapply(1:length(coefList), function(j) {
  var_names <- names(rList[[1]][[j]])
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
        file = paste("PowerAnalysis-aim3-ntp",ntp,".RDS",sep=""))