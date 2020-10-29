## ----------------------------------------------------------------- ##
## post_experiment.R ----------------------------------------------- ##
## Author: Peter Norwood, NC State University ---------------------- ##
## Purpose: take the trial information and assign estimated -------- ##
## optimal interventions to a group of randomly selected new ------- ##
## subjects to evaluate the strength of information within --------- ##
## the experiment -------------------------------------------------- ##
## ----------------------------------------------------------------- ##

## load functions
setwd("~/Research/Written Prelim/WrittenPrelim")
source("funcs.R")

library(foreach)
library(doParallel)

## post_experiment
## Purpose: evaluate post experiment performance
## param exp_data: data from experiment
## param post_data: post experiment sample
## param p: dimension of contexts (not including intercept)
## param K: number of arms
## param theta: true mean parameter vector
## return output: list with (1) dataframe with correct and regret
##                (2) final norm between theta and theta_hat
post_experiment <- function(exp_data,post_data,p,K,theta){
  
  ## fit the model on the experiment data
  X_temp <- exp_data[1:N_post,1:p]
  A_temp <- exp_data[1:N_post,p+1]
  Y <- exp_data[1:N_post,p+3]
  
  temp <- data.frame(X_temp,A=A_temp,Y)
  
  fit <- lm(Y~-1+as.factor(A)+as.factor(A):.-
            as.factor(A):A,
            data=temp)
  
  ## assign intercention to out of experiment subjects
  perf <- foreach(i=1:nrow(post_data),.combine=rbind) %dopar% {
    context <- post_data[i,1:p]
    dat <- matrix(NA,ncol=p+3,nrow=K)
    #dat <- as.data.frame(dat)
    colnames(dat) <- c(colnames(context),"A","mu_hat","mu")
    for(k in 1:K){
      dat[k,1:p] <- as.matrix(context)
      dat[k,p+1] <- k
      ## predicted outcome
      dat[k,p+2] <- predict(fit,newdat=data.frame(t(dat[k,1:p]),A=k,Y=0))
      ## true outcome               
      dat[k,p+3] <- mean_outcome(X=dat[k,1:p],A=c(1:K),a=k,theta=theta)
    }
    dat <- as.data.frame(dat)
    ## save data
    perf_temp <- c()
    ## subject
    perf_temp[1] <- i
    ## did we get it right
    perf_temp[2] <- ifelse(which.max(dat$mu_hat)==which.max(dat$mu),1,0)
    ## regret
    perf_temp[3] <- max(dat$mu) - dat$mu[which.max(dat$mu_hat)]
    print(perf_temp)
    
    ## rbind
    perf_temp
  }
  
  perf <- as.data.frame(perf)
  colnames(perf) <- c("sub","correct","regret")
  
  return(perf)
  
}