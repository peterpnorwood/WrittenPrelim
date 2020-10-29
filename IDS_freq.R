## ----------------------------------------------------------------- ##
## IDS_freq.R ------------------------------------------------------ ##
## Author: Peter Norwood, NC State University ---------------------- ##
## Purpose: run an experiment with a frequentist analog of IDS ----- ##
## ----------------------------------------------------------------- ##

## load functions
setwd("~/Research/Written Prelim/WrittenPrelim")
source("funcs.R")

## IDS_freq
## Purpose: run an experiment with IDS (frequentist analog)
## param train_set: dataset with context for N individuals
## param burn_in: sample size of simple randomization
## param A: vector of possible treatments
## param theta: true mean outcome parameter vector
## param sd_Y: standard deviation for response
## return dat: dataframe with X,A,mu,Y,regret,norm
IDS_freq <- function(train_set,burn_in,A,theta,sd_Y){
  
  ## number of subjects
  N <- nrow(train_set)
  ## dimension of context
  p <- ncol(train_set)-3
  ## number of arms
  K <- length(A)
  
  ## trial dataset
  dat <- matrix(NA,nrow=N,ncol=p+5)
  ## context
  dat[1:N,1:p] <- as.matrix(train_set)[1:N,1:p]
  ## first burn_in interventions
  dat[1:burn_in,p+1] <- train_set$A[1:burn_in]
  ## first burn_in means
  dat[1:burn_in,p+2] <- train_set$mu[1:burn_in]
  ## first burn_in outcomes
  dat[1:burn_in,p+3] <- train_set$Y[1:burn_in]
  ## name the same colnames
  colnames(dat) <- c(colnames(train_set),"regret","norm")
  
  ## loop through the new patients
  for(i in (burn_in+1):N){
    
    ## fit the outcome model
    X_temp <- dat[1:(i-1),1:p]
    A_temp <- dat[1:(i-1),p+1]
    Y <- dat[1:(i-1),p+3]
    
    temp <- data.frame(X_temp,A=A_temp,Y)
    
    fit <- lm(Y~-1+as.factor(A)+as.factor(A):.-
                as.factor(A):A,
              data=temp)
    
    ## gather parameter convergence information
    coef_fit <- coef(fit)
    theta_hat <- c()
    ## put them in the same format as the theta vector
    tik <- 1
    for(ii in 1:K){
      for(jj in 0:p){
        theta_hat[tik] <- coef_fit[ii+(K)*jj]
        tik=tik+1
      }
    }
    
    ## measure the euclidean norm between theta and theta_hat
    dat[i,ncol(dat)] <- norm(matrix(theta-theta_hat),type="F")
  
    ## loop through interventions to find greedy intevention
    info <- matrix(NA,nrow=length(A),ncol=4)
    
    ## determinant before addition
    prev_X <- model.matrix(fit)
    det_prev_XtX <- det(t(prev_X) %*% prev_X)
    
    tick=1
    for(a in A){
      ## gather det if a is assigned
      temp2 <- data.frame(t(dat[i,1:p]),A=a,Y=0)
      temp_X <- model.matrix(fit,data=rbind(temp,temp2))
      det_XtX <- det(t(temp_X) %*% temp_X)
      info_gain <- det_XtX/det_prev_XtX
      ## reward + (alpha/t)*d_k
      mu_hat <- predict(fit,temp2)
      ## true mean outcome given a
      mu <- mean_outcome(X=dat[i,1:p],A=A,a=a,theta=theta)
      ## save info
      info[tick,] <- c(a,mu_hat,info_gain,mu)
      tick=tick+1
    }
    ## save info as dataframe
    info <- data.frame(info)
    colnames(info) <- c("A","mu_hat","info_gain","mu")
    info$cost <- (max(info$mu_hat)+1)-info$mu_hat
    info$ir <- info$cost/info$info_gain
    
    ## assign intervention
    dat[i,p+1] <- info$A[which.min(info$ir)]
    ## find mean outcome
    dat[i,p+2] <- info$mu[which.min(info$ir)]
    ## find outcome
    dat[i,p+3] <- rnorm(1,dat[i,p+2],sd_Y)
    ## find regret
    dat[i,p+4] <- max(info$mu) - dat[i,p+2]
  }
  
  dat <- data.frame(dat)
  dat$sub <- 1:nrow(dat)
  return(dat)
  
}


# p <- 5
# K=5
# theta <- rnorm((p+1)*K,0,1)
# train_set <- gen_data(N=500,p=p,sd_X=0.5,A=1:K,sd_Y=1,theta=theta)
# test_IDS <- IDS_freq(train_set=train_set,burn_in=(p+1)*K*3,A=1:K,theta=theta,sd_Y=1)
# #
# # hist(test_greedy$regret)
# # #
# ggplot(data=test_IDS[((p+1)*K*3+1):nrow(test_IDS),])  +
#     geom_line(aes(x=sub,y=norm))
