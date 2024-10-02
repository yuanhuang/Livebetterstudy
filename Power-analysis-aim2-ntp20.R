########################################################################################
# Simulation-based power analysis
# LIVE BETTER Trial Protein biomarker study
# Aim 2 - linear regression model (ordinal outcome)
# Date: 2024.09.2x
# Technical notes: 
# 1. As an ad-hoc approach, we fit the linear regression model for this simulation.
#   For marginal model, we also fit the ordered logistic model, 
#        whose results are very close to the linear regression model.
#      The reason is that we generated the outcome from a continuous one,
#        and then assign levels by the specified quartiles.
#      Results are not reported due to limited space, but the codes are provided.
#  For penalization method, we also adopt the linear regression. 
#     This is due to the findings from the marginal model for its similar
#       performance in terms of variable selection.
#     The package "ordinalNet" can be applied to the ordinal outcome, but 
#       it runs slowly. Hence it was not adopted for this simulation.
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
# simulate data and fit LASSO for each replicate
########################################################################################
# - trt variable will not be penalized
# - tuning is selected by $lambda.min 
#   other option is $lambda.1se
# - all other settings are as default.
########################################################################################


# set up different coef

coefList     <- list()

coefList[[1]]<- 
  c(runif(n = (ntp/2),0.2,0.3),
    runif(n = (ntp/2),min=-0.3,max=-0.2))

coefList[[2]]<- 
  c(runif(n = (ntp/2),0.3,0.4),
    runif(n = (ntp/2),min=-0.4,max=-0.3))

coefList[[3]]<- 
  c(runif(n = (ntp/2),0.4,0.5),
    runif(n = (ntp/2),min=-0.5,max=-0.4))

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
  colnames(dt) = c("trt",paste("protein",1:p,sep="_"))
  
  # simulate responses for logistic model and calculate the estimates
  # define the function
  
  cal.fun <- function(coefinput){
    
    ## construct data
    p <- dim(dt.prot)[2]
    
    # assume half of the magnitude for interaction terms
    mean.xb   <- trt+dt.prot[,1:ntp] %*% matrix(coefinput,ncol=1)
    y         <- mean.xb + rnorm(n=n,mean = 0,sd = 1.6)
    
    quant_vec <- quantile(y, probs =c(0,0.28,0.4+0.28,0.28+0.4+0.28,1))
    ind1 <- which(y <= quant_vec[2] & y >= quant_vec[1])
    ind2 <- which(y <= quant_vec[3] & y >= quant_vec[2])
    ind3 <- which(y <= quant_vec[4] & y >= quant_vec[3])
    ind4 <- which(y <= quant_vec[5] & y >= quant_vec[4])
    y[ind1] <- 1; y[ind2] <- 2; y[ind3] <- 3; y[ind4] <- 4
    
    ## fit lasso
    # ordinalNet package run slowly 
    # use LASSO to estimate model parameters
    # check   plot(fit.high) for the first few cases to make sure the
    # turning is reasonably selected
    
    fit <- cv.glmnet(dt, y,
                     penalty.factor = c(0, rep(1, p)),
                     penfamily = gussian())
    
    fit.sel <- rownames(coef(fit, lambda=fit$lambda.min))[which(coef(fit, lambda=fit$lambda.min)!=0)]
    
    fit.sel <- grep(pattern = "protein_",x = fit.sel,value = T)
    fit.selIndex <- as.numeric( gsub(pattern = "protein_",replacement = "",x = fit.sel))
    
    ntp.pen <- sum(fit.selIndex <= ntp)
    nfp.pen <- length(fit.selIndex) - ntp.pen
    
    ## fit marginal model
    ###### Note, we report the linear regression model output
    ###### However, we have tested under the ordered logistic model
    ######    whose results are very similar to those from linear model.
    
    margin.p <- sapply(X = 1:p,
                       FUN = function(i){summary(
                         lm(y ~ trt + dt.prot[, i]))$coef[3, "Pr(>|t|)"]})
    
    ########## Run the following codes for the ordered logistic model
    ## LR test would be preferred in real data
    # sapply(X = 1:p,
    #        FUN = function(i){
    #          fit1=MASS::polr(factor(y) ~ trt + dt.prot[, i],Hess= T)
    #          fit2=MASS::polr(factor(y) ~ trt,Hess= T)
    #          anova(fit1,fit2)[["Pr(Chi)"]][2]
    #        })
    # Wald test
    # margin.p  <- sapply(X = 1:p,
    #              FUN = function(i){
    #                fit1=MASS::polr(factor(y) ~ trt + dt.prot[, i],Hess= T)
    #                pnorm(abs(coef(summary(fit1))[2,"t value"]), 
    #                     lower.tail = FALSE) * 2
    #              })
    
    ntp.mar = sum(which(p.adjust(margin.p, method = "fdr") < 0.1) <= ntp)
    nfp.mar = length(which(p.adjust(margin.p, method = "fdr") < 0.1)) - ntp.mar
    
    return(data.frame(
      ntp.pen = ntp.pen,
      nfp.pen = nfp.pen,
      ntp.mar = ntp.mar,
      nfp.mar = nfp.mar,
      rtp.pen = ntp.pen/ntp,
      rfp.pen = nfp.pen/(p-ntp),
      rtp.mar = ntp.mar/ntp,
      rfp.mar = nfp.mar/(p-ntp) 
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
        file = paste("PowerAnalysis-aim2-ntp",ntp,".RDS",sep=""))
