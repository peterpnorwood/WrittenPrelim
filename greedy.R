## ----------------------------------------------------------------- ##
## greedy.R -------------------------------------------------------- ##
## Author: Peter Norwood, NC State University ---------------------- ##
## Purpose: run an experiment where A = argmax mu_hat -------------- ##
## ----------------------------------------------------------------- ##

## load functions
setwd("~/Research/Written Prelim/WrittenPrelim")
source("funcs.R")

## greedy
## Purpose: run a greedy d-optimal experiment
## param train_set: dataset with context for N individuals
## param burn_in: sample size of simple randomization
## param A: vector of possible treatments
## param theta: true mean outcome parameter vector
## param sigma: standard deviation for response
## return dat: dataframe with X,A,mu,Y,regret,norm
greedy <- function(train_set,burn_in,A,theta,sigma){
  
  ## create trial dataset
  N <- nrow(train_set)
  dat <- matrix(NA,nrow=N,ncol=6)
  ## context
  dat[,1] <- train_set$X
  ## first burn_in interventions
  dat[1:burn_in,2] <- train_set$A[1:burn_in]
  ## first burn_in means
  dat[1:burn_in,3] <- train_set$mu[1:burn_in]
  ## first burn_in outcomes
  dat[1:burn_in,4] <- train_set$Y[1:burn_in]
  ## name the same colnames
  colnames(dat) <- c(colnames(train_set),"regret","norm")
  
  ## loop through the new patients
  for(i in (burn_in+1):N){
    
    ## fit the outcome model
    prev <- data.frame(dat[1:(i-1),])
    fit <- lm(Y~-1+as.factor(A) + as.factor(A):(X+I(X**2)),
              data=prev)
    
    ## gather parameter convergence information
    coef_fit <- coef(fit)
    theta_hat <- c(coef_fit[1],coef_fit[6],coef_fit[11],
                   coef_fit[2],coef_fit[7],coef_fit[12],
                   coef_fit[3],coef_fit[8],coef_fit[13],
                   coef_fit[4],coef_fit[9],coef_fit[14],
                   coef_fit[5],coef_fit[10],coef_fit[15])
    
    ## measure the euclidean norm between theta and theta_hat
    dat[i,6] <- norm(matrix(theta-theta_hat),type="F")
    
    
    ## loop through interventions to find greedy intevention
    info <- matrix(NA,nrow=length(A),ncol=4)
    tick=1
    for(a in A){
      ## gather det if a is assigned
      temp <- data.frame(X=dat[i,1],A=a,Y=0,mu=0,regret=0,norm=0)
      ## estiamted mean outcome given a
      mu_hat <- predict(fit,temp)
      ## true mean outcome given a
      mu <- mean_outcome(X=dat[i,1],A=a,theta=theta)
      ## save info
      info[tick,] <- c(a,NA,mu_hat,mu)
      tick=tick+1
    }
    ## save info as dataframe
    info <- data.frame(info)
    colnames(info) <- c("A","det","mu_hat","mu")
    
    ## assign intervention
    dat[i,2] <- info$A[which.max(info$mu_hat)]
    ## find mean outcome
    dat[i,3] <- info$mu[which.max(info$mu_hat)]
    ## find outcome
    dat[i,4] <- rnorm(1,dat[i,3],sigma)
    ## find regret
    dat[i,5] <- max(info$mu) - dat[i,4]
  }
  
  dat <- data.frame(dat)
  dat$sub <- 1:nrow(dat)
  return(dat)
  
}


theta1 <- c(3.25,0.01,0.0)
theta2 <- c(2.25,-0.5,0.01)
theta3 <- c(3.5,0.05,-0.50)
theta4 <- c(3.15,-0.40,-0.15)
theta5 <- c(2.5,0.5,-0.05)
theta <- c(theta1,theta2,theta3,theta4,theta5)
sigma=0.2

train_set <- gen_data(N=500,lower=-3.0,upper=3.0,A=1:5,theta=theta,sigma=sigma)
test_greedy <- greedy(train_set=train_set,burn_in=50,A=1:5,
               theta=theta,sigma=sigma)

ggplot(data=test_greedy[51:nrow(test_greedy),])  +
  geom_line(aes(x=sub,y=norm))
