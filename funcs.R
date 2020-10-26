## ----------------------------------------------------------------- ##
## funcs.R --------------------------------------------------------- ##
## Author: Peter Norwood, NC State University ---------------------- ##
## Purpose: general functions used for the study ------------------- ##
## ----------------------------------------------------------------- ##


## mean_outcome
## Purpose: generate mean outcome
## param X: context
## param A: intervention
## param theta: mean parameter vector
## return mu: mean response
mean_outcome <- function(X,A,theta){
  
  ## logic to decide mean response based on A
  if(A==1){
      mu <- theta[1] + theta[2]*X + theta[3]*I(X**2)
  } else if(A==2){
    mu <- theta[4] + theta[5]*X + theta[6]*I(X**2)
  } else if(A==3){
    mu <- theta[7] + theta[8]*X + theta[9]*I(X**2)
  } else if(A==4){
    mu <- theta[10] + theta[11]*X + theta[12]*I(X**2)
  } else{
    mu <- theta[13] + theta[14]*X + theta[15]*I(X**2)
  }
  
  ## return the mean outcome
  return(mu)
  
}

## vectorize mean outcome for X and A
mean_outcome <- Vectorize(mean_outcome,vectorize.args = c("X","A"))

## gen_data
## Purpose: generate datasets
## param N: sample size
## param lower,upper: context x ~ unif(lower,upper)
## param A: vector of possible interventions
## param theta: mean parameter vector
## param sigma: standard deviation of random error
## return dat: dataset with columns: X, A, mu, regret, Y
gen_data <- function(N,lower,upper,A,theta,sigma){
  
  ## generate context
  X <- runif(N,lower,upper)
  ## create randomly assigned A
  A_vec <- sample(A,N,replace=TRUE)
  ## mean outcome
  mu <- mean_outcome(X,A_vec,theta)
  ## true outcome
  Y <- rnorm(N,mu,sigma)
  
  ## dataset to return
  dat <- data.frame(X,A=A_vec,mu,Y)
  
  return(dat)
  
}


# theta1 <- c(3.25,0.01,0.0)
# theta2 <- c(2.25,-0.5,0.01)
# theta3 <- c(3.5,0.05,-0.50)
# theta4 <- c(3.15,-0.40,-0.15)
# theta5 <- c(2.5,0.5,-0.05)
# theta <- c(theta1,theta2,theta3,theta4,theta5)
# sigma=0.2
# 
# train_set <- gen_data(N=500,lower=-3.0,upper=3.0,A=1:5,theta=theta,sigma=sigma)
# ggplot(data=dat) +
#   geom_line(aes(x=X,y=Y,color=as.factor(A),group=as.factor(A)))
# 
# fit <- lm(Y~-1+as.factor(A) + as.factor(A):(X+I(X**2)),data=dat)
# X <- model.matrix(fit)
# XtX <- t(X) %*% X
# lst[i] <- det(XtX)
# 
# 
# coef_fit <- coef(fit)
# theta_hat <- c(coef_fit[1],coef_fit[6],coef_fit[11],
#                coef_fit[2],coef_fit[7],coef_fit[12],
#                coef_fit[3],coef_fit[8],coef_fit[13],
#                coef_fit[4],coef_fit[9],coef_fit[14],
#                coef_fit[5],coef_fit[10],coef_fit[15])
# norm(theta-theta_hat,"2")
