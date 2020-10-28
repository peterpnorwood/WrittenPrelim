## ----------------------------------------------------------------- ##
## funcs.R --------------------------------------------------------- ##
## Author: Peter Norwood, NC State University ---------------------- ##
## Purpose: general functions used for the study ------------------- ##
## ----------------------------------------------------------------- ##

library(MASS)


## make_design
make_design <- function(K,p,a,X){
  
  ## create design
  design <- c()
  for(j in 1:K){
    if(a==j){
      x <- c(1,X)
    }else{
      x <- rep(0,p+1)
    }
    design <- c(design,x)
  }
  
  return(design)
  
}

## mean_outcome1
## Purpose: generate mean outcome
## param X: p-dimensional context
## param a: intervention
## param A: possible interventions
## param theta: mean parameter vector
## return mu: mean response
mean_outcome1 <- function(X,a,A,theta){
  
  ## dimensions
  K <- length(A)
  p <- length(X)
  
  design <- make_design(K=K,p=p,a=a,X=X)
  
  mu <- design %*% theta
  
  ## return the mean outcome
  return(mu)
  
}

## "vectorize" mean outcome for X and A
mean_outcome <- function(X,a,A,theta){
  N <- max(1,nrow(X))
  
  if(N==1){
    mu <- mean_outcome1(X=X,a=a,A=A,theta=theta)
  } else{
    mu <- c()
    ## loop through rows of x
    for(n in 1:N){
      mu[n] <- mean_outcome1(X=X[n,],a=a[n],A=A,theta=theta)
    }
  }
  return(mu)
}

## gen_data
## Purpose: generate datasets
## param N: sample size
## param p: dimension of the context
## param A: vector of possible interventions
## param theta: mean parameter vector
## param sigma: standard deviation of random error
## return dat: dataset with columns: X, A, mu, Y
gen_data <- function(N,p,sd_X,A,sd_Y,theta){
  
  ## generate context
  sd_X <- sd_X*sqrt(p+1)
  X <- mvrnorm(N,rep(0,p),diag(p))
  
  ## create randomly assigned A
  A_vec <- sample(A,N,replace=TRUE)
  
  ## mean outcome
  mu <- mean_outcome(X=X,a=A_vec,A=A,theta=theta)
  ## true outcome
  Y <- rnorm(N,mu,sd_Y)
  
  ## dataset to return
  dat <- data.frame(X,A=A_vec,mu,Y)
  
  return(dat)
  
}

#test <- gen_data(N=1000,p=5,sd_X=0.5,A=1:5,sd_Y=1,theta=mvrnorm(n=1,rep(0,25),diag(25)))

