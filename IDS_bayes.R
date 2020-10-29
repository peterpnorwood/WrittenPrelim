## ----------------------------------------------------------------- ##
## IDS_bayes.R ----------------------------------------------------- ##
## Author: Peter Norwood, NC State University ---------------------- ##
## Purpose: run an experiment with IDS ----------------------------- ##
## ----------------------------------------------------------------- ##

## load functions
setwd("~/Research/Written Prelim/WrittenPrelim")
source("funcs.R")
library(rstanarm)
library(MASS)


## calc_ir
## Purpose: calculate information ratio
## param M: samples to take from posterior distribution
## param mu: posterior mean
## param Sigma: posterior covariance matrix
## param A: posterior of interventions
## param context: context for subject
## return ir: column of information ratios for each individual
calc_ir <- function(M,mu,Sigma,A,context){
  
  ## sample from posterior
  posterior_sample <- mvrnorm(n=M,mu=mu,Sigma=Sigma)
  
  ## create mu_hat
  mu_hat <- (1/M)*colSums(posterior_sample)
  
  ## see how often each piece is maximized
  maxs <- c()
  for(m in 1:M){
    ## get estimated value for each a
    ests <- c()
    for(a in A){
      ## rearrange context
      ctx <- make_design(K=K,p=p,a=a,X=context)
      ## estimated mean response
      ests[a] <- ctx %*% posterior_sample[m,]
    }
    ## save max
    maxs[m] <- which.max(ests)
  }
  max_df <- data.frame(max=maxs,m=1:M)
  
  ## create p_star, mu_starm, L, rho_star
  p_star <- c()
  mu_a <- list()
  L <- matrix(0,nrow=length(mu),ncol=length(mu))
  rho_star <- 0
  for(a in A){
    p_star[a] <- length(maxs[maxs==a])/M
    ## create mu_a and L
    if(p_star[a]==0){
      mu_a[[a]] <- 0
    }else if(p_star[a]==1/M){
      ## subset the ms
      ms <- max_df$m[max_df$max==a]
      mu_a[[a]] <- posterior_sample[ms,]
      ## create L
      L <- L + p_star[a]*( (mu_a[[a]]-mu_hat) %*% t((mu_a[[a]]-mu_hat)))
      ## create rho_star
      ctx <- make_design(K=K,p=p,a=a,X=context)
      rho_star <- rho_star + p_star[a]*(ctx %*% mu_a[[a]])
    }else{
      ## subset the ms
      ms <- max_df$m[max_df$max==a]
      mu_a[[a]] <- colSums(posterior_sample[ms,])/length(ms)
      ## create L
      L <- L + p_star[a]*( (mu_a[[a]]-mu_hat) %*% t((mu_a[[a]]-mu_hat)))
      ## create rho_star
      ctx <- make_design(K=K,p=p,a=a,X=context)
      rho_star <- rho_star + p_star[a]*(ctx %*% mu_a[[a]])
    }

  }
  
  ## ceate v_a, delta_a, and information ratio
  v_a <- c()
  delta_a <- c()
  ir <- c()
  for(a in A){
    ctx <- make_design(K=K,p=p,a=a,X=context)
    v_a[a] <- t(ctx) %*% (L %*% ctx) 
    delta_a[a] <- rho_star - ctx %*% mu_hat
    ir[a] <- (delta_a[a]**2)/v_a[a]
  }
  
  return(ir)
  
}

## IDS_bayes
## Purpose: run an experiment with IDS (frequentist analog)
## param train_set: dataset with context for N individuals
## param burn_in: sample size of simple randomization
## param A: vector of possible treatments
## param theta: true mean outcome parameter vector
## param sd_Y: standard deviation for response
## param M: samples to take from posterior distribution
## return dat: dataframe with X,A,mu,Y,regret,norm
IDS_bayes <- function(train_set,burn_in,A,theta,sd_Y,M){
  
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
    
    ## fit linear model with default priors
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
    
    ## calcualte information ratio
    ir <- calc_ir(M=M,mu=theta_hat,Sigma=vcov(fit),A=1:K,context=dat[i,1:p])
    ## calculate true mean outcomes
    mean_outcomes <- c()
    for(a in A){
      mean_outcomes[a] <- mean_outcome(X=dat[i,1:p],A=A,a=a,theta=theta)
    }
    ## save info as dataframe
    info <- data.frame(A=1:K,ir=ir,mu=mean_outcomes)
    
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
# test_IDS <- IDS_bayes(train_set=train_set,burn_in=(p+1)*K*3,A=1:K,theta=theta,sd_Y=1,
#                       M=100)
# #
# # hist(test_greedy$regret)
# # #
# ggplot(data=test_IDS[((p+1)*K*3+1):nrow(test_IDS),])  +
#     geom_line(aes(x=sub,y=norm))