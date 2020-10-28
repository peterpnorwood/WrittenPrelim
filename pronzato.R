## ----------------------------------------------------------------- ##
## pronzato.R ------------------------------------------------------ ##
## Author: Peter Norwood, NC State University ---------------------- ##
## Purpose: run an experiment where -------------------------------- ##
## A = argmax mu_hat + alpha*g, from Pronzato (2000) --------------- ##
## ----------------------------------------------------------------- ##

## load functions
setwd("~/Research/Written Prelim/WrittenPrelim")
source("funcs.R")

## alpha_t
## Purpose: generate alpha_t for the pronzato algorithm
##          that satisfies theorem conditions
## param alpha: some value > 0
## param t: time step
alpha_t <- function(alpha,t){
  alpha_t <- alpha*sqrt(t)*log(t)
  return(alpha_t)
}

## pronzato
## Purpose: run a experiment following pronzato (2000)'s method
## param train_set: dataset with context for N individuals
## param burn_in: sample size of simple randomization
## param A: vector of possible treatments
## param theta: true mean outcome parameter vector
## param sd_Y: standard deviation for response
## param al: parameter that represents compromise
## return dat: dataframe with X,A,mu,Y,regret,norm
pronzato <- function(train_set,burn_in,A,theta,sd_Y,al){
  
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
    
    ## get "I-inv" matrix
    sigma_hat <- sigma(fit)
    Iinv <- (i-1)*vcov(fit)/I(sigma_hat**2)
    
    ## alpha
    alpha <- alpha_t(al,i)
    
    ## loop through interventions to find greedy intevention
    info <- matrix(NA,nrow=length(A),ncol=3)
    tick=1
    for(a in A){
      
      ## gather det if a is assigned
      temp <- data.frame(t(dat[i,1:p]),A=a,mu=0,Y=0,regret=0,norm=0)
      fx <- model.matrix(fit,data=temp)
      ## reward + (alpha/t)*d_k
      mu_hat <- predict(fit,temp)
      z <- mu_hat + (alpha/(i-1))*( fx %*% (Iinv %*% t(fx)) )
      ## true mean outcome given a
      mu <- mean_outcome(X=dat[i,1:p],A=A,a=a,theta=theta)
      ## save info
      info[tick,] <- c(a,z,mu)
      tick=tick+1
    }
    ## save info as dataframe
    info <- data.frame(info)
    colnames(info) <- c("A","z","mu")
    
    ## assign intervention
    dat[i,p+1] <- info$A[which.max(info$z)]
    ## find mean outcome
    dat[i,p+2] <- info$mu[which.max(info$z)]
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
# test_pronzato <- pronzato(train_set=train_set,burn_in=(p+1)*K*3,
#                         A=1:K,theta=theta,sd_Y=1,al=0.5)
# #
# # hist(test_greedy$regret)
# # #
# ggplot(data=test_pronzato[((p+1)*K*3+1):nrow(test_pronzato),])  +
#     geom_line(aes(x=sub,y=norm))
