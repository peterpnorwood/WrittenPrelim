## ----------------------------------------------------------------- ##
## greedy_first.R -------------------------------------------------- ##
## Author: Peter Norwood, NC State University ---------------------- ##
## Purpose: run an experiment where A = argmax mu_hat -------------- ##
## ----------------------------------------------------------------- ##

## load functions
setwd("~/Research/Written Prelim/WrittenPrelim")
source("funcs.R")

## greedy_first
## Purpose: run an experiment with the GF algoirthm
## param train_set: dataset with context for N individuals
## param burn_in: sample size of simple randomization
## param A: vector of possible treatments
## param theta: true mean outcome parameter vector
## param sd_Y: standard deviation for response
## param t0: number of obsrvations (after burn in)
##           before we start checking to move from greedy
## param eps: P(explore) if we switch to eps greedy
## return dat: dataframe with X,A,mu,Y,regret,norm,R
greedy_first <- function(train_set,burn_in,A,theta,sd_Y,t0,eps){
  
  
  ## number of subjects
  N <- nrow(train_set)
  ## dimension of context
  p <- ncol(train_set)-3
  ## number of arms
  K <- length(A)
  
  ## trial dataset
  dat <- matrix(NA,nrow=N,ncol=p+6)
  ## context
  dat[1:N,1:p] <- as.matrix(train_set)[1:N,1:p]
  ## first burn_in interventions
  dat[1:burn_in,p+1] <- train_set$A[1:burn_in]
  ## first burn_in means
  dat[1:burn_in,p+2] <- train_set$mu[1:burn_in]
  ## first burn_in outcomes
  dat[1:burn_in,p+3] <- train_set$Y[1:burn_in]
  ## name the same colnames
  colnames(dat) <- c(colnames(train_set),"regret","norm","R")
  
  ## parameter to see if we stick with greedy or move to e-greedy
  R=0
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
    dat[i,ncol(dat)-1] <- norm(matrix(theta-theta_hat),type="F")
    
    
    ## loop through interventions to find greedy intevention
    info <- matrix(NA,nrow=length(A),ncol=3)
    tick=1
    for(a in A){
      ## gather det if a is assigned
      temp <- data.frame(t(dat[i,1:p]),A=a,mu=0,Y=0,regret=0,norm=0)
      ## estiamted mean outcome given a
      mu_hat <- predict(fit,temp)
      ## true mean outcome given a
      mu <- mean_outcome(X=dat[i,1:p],A=A,a=a,theta=theta)
      ## save info
      info[tick,] <- c(a,mu_hat,mu)
      tick=tick+1
    }
    ## save info as dataframe
    info <- data.frame(info)
    colnames(info) <- c("A","mu_hat","mu")
    
    ## conditional logic to select action
    if(i < t0){
      
      ## assign intervention (greedy)
      dat[i,p+1] <- info$A[which.max(info$mu_hat)]
      ## find mean outcome
      dat[i,p+2] <- info$mu[which.max(info$mu_hat)]
      ## find outcome
      dat[i,p+3] <- rnorm(1,dat[i,p+2],sd_Y)
      ## find regret
      dat[i,p+4] <- max(info$mu) - dat[i,p+2]
      
    }else if(i == t0){
      
      ## assign intervention (greedy)
      dat[i,p+1] <- info$A[which.max(info$mu_hat)]
      ## find mean outcome
      dat[i,p+2] <- info$mu[which.max(info$mu_hat)]
      ## find outcome
      dat[i,p+3] <- rnorm(1,dat[i,p+2],sd_Y)
      ## find regret
      dat[i,p+4] <- max(info$mu) - dat[i,p+2]
      
      ## get lambda0
      X_t0 <- model.matrix(fit)
      M_t0 <- t(X_t0) %*% X_t0
      eig_M <- eigen(M_t0)
      lambda0 <- min(eig_M$values)/(2*t0)
      
    } else if(i > t0){
      
      ## if there was covariate diversity at the last time step
      ## then check again
      ## if there wasn't then no need to check
      if(R==0){
        ## assign intervention (greedy)
        dat[i,p+1] <- info$A[which.max(info$mu_hat)]
        ## find mean outcome
        dat[i,p+2] <- info$mu[which.max(info$mu_hat)]
        ## find outcome
        dat[i,p+3] <- rnorm(1,dat[i,p+2],sd_Y)
        ## find regret
        dat[i,p+4] <- max(info$mu) - dat[i,p+2]
        dat[i,p+6] <- R
      
        ## check to see if we have covariate diversity
        X_t1 <- model.matrix(fit)
        M_t1 <- t(X_t0) %*% X_t0
        eig_M_t1 <- eigen(M_t1)
        lambda1 <- min(eig_M_t1$values)
        R <- R+ifelse(lambda1 < lambda0*(i-1)*0.25,1,0)
      
      } else{ 
        ## randomize via e-greedy
        prob_vec <- rep(eps/(length(A)-1),length(A))
        prob_vec[which.max(info$mu_hat)] <- 1-eps
        
        ## assign intervention (e-greedy)
        dat[i,p+1] <- sample(A,1,prob=prob_vec)
        ## find mean outcome
        dat[i,p+2] <- info$mu[dat[i,p+1]]
        ## find outcome
        dat[i,p+3] <- rnorm(1,dat[i,p+2],sd_Y)
        ## find regret
        dat[i,p+4] <- max(info$mu) - dat[i,p+2]
        dat[i,p+6] <- R
      }
      
    }  

  }
  
  dat <- data.frame(dat)
  dat$sub <- 1:nrow(dat)
  return(dat)
  
}

# p <- 2
# K=8
# A=1:K
# theta <- rnorm((p+1)*K,0,1)
# N=2000
# sd_Y <- 0.5*sqrt(p)
# burn_in=(p+1)*K*3
# t0 = burn_in+100
# eps=0.05
# train_set <- gen_data(N=2000,p=p,sd_X=0.5,A=1:K,sd_Y=sd_Y,theta=theta)
# test_gf <- greedy_first(train_set=train_set,burn_in=(p+1)*K*3,A=1:K,
#                         theta=theta,t0=100,sd_Y=1,eps=0.1)
# # #
# hist(test_greedy$regret)
# # # #
# ggplot(data=test_gf[((p+1)*K*2+1):nrow(test_gf),])  +
#      geom_line(aes(x=sub,y=norm))

     