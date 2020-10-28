## ----------------------------------------------------------------- ##
## sims.R ---------------------------------------------------------- ##
## Author: Peter Norwood, NC State University ---------------------- ##
## Purpose: simulation experiments with different methods ---------- ##
## and store their results in a user-friendly format --------------- ##
## ----------------------------------------------------------------- ##

## load functions
setwd("~/Research/Written Prelim/WrittenPrelim")
source("funcs.R")
source("pronzato.R")
source("IDS_freq.R")
source("greedy_first.R")
source("e_greedy.R")
source("d_opt.R")
source("greedy.R")

## load packages
library(parallel)
library(tidyverse)
library(MASS)


## sim
## Purpose: simulate all methods on the same dataset
## param N: total sample size
## param lower, upper: X ~ unif(lower,upper)
## param burn_in: length of simple rand until 
##                 adaptive methods begin
## param A: vector of treatments
## param t0: time for greedy first to stay with 
##           greedy before checking
## param eps: P(explore) = eps in greedy first if it switches
## param al: value for weight put on penalty in pronzato method
## return output: list of trial information from each method
sim <- function(N,mean,sd,
                burn_in,A,sigma,t0,eps,al){
  
  #theta <- mvrnorm(n=1,mu=rep(0,15),Sigma=diag(15)
  
  ## generate training set
  train_set <- gen_data(N,mean,sd,A,theta,sigma)
  
  ## run pronzato experiment
  pronzato_sim <- try(pronzato(train_set=train_set,burn_in=burn_in,
                      A=A,theta=theta,sigma=sigma,al=al))
  
  ## run d-opt experiment
  d_opt_sim <- try(d_opt(train_set=train_set,burn_in=burn_in,
                   A=A,theta=theta,sigma=sigma))
  
  ## run e_greedy experiment
  e_greedy_sim <- try(e_greedy(train_set=train_set,burn_in=burn_in,
                      sigma=sigma,A=A,theta=theta,eps=eps))
  
  
  ## run greedy experiment
  greedy_sim <- try(greedy(train_set=train_set,burn_in=burn_in,
                           sigma=sigma,A=A,theta=theta))
  
  ## run GF experiment
  greedy_first_sim <- try(greedy_first(train_set=train_set,burn_in=burn_in,
                           sigma=sigma,A=A,theta=theta,t0=t0,eps=eps))
  
  ## run IDS experiment
  IDS_freq_sim <- try(IDS_freq(train_set=train_set,burn_in=burn_in,
                               sigma=sigma,A=A,theta=theta))
  
  ## save output in a list
  output <- list(pronzato=pronzato_sim,
                 d_opt=d_opt_sim,
                 e_greedy=e_greedy_sim,
                 greedy=greedy_sim,
                 greedy_first=greedy_first_sim,
                 IDS_freq=IDS_freq_sim)
  
  return(output)
  
}

## run the simulations

## parameters
# theta1 <- c(3.25,0.01,0.0)
# theta2 <- c(2.25,-0.5,0.01)
# theta3 <- c(3.5,0.05,-0.50)
# theta4 <- c(3.15,-0.40,-0.15)
# theta5 <- c(2.5,0.5,-0.05)

# theta1 <- c(3.25,0.01,0.0)
# theta2 <- c(2.25,-0.25,0.01)
# theta3 <- c(3.75,0.05,-0.20)
# theta4 <- c(3.15,-0.40,-0.15)
# theta5 <- c(2.5,0.05,-0.05)
# theta <- c(theta1,theta2,theta3,theta4,theta5)
theta <- c(1.310,-1.139,-0.483,0.614,1.868,-0.650,  
           0.505,-0.308,-1.263,0.970,0.345,0.184,
           0.148,0.931,-0.625)
A=1:5
sigma=0.25*sqrt(5)
lower=-1
upper=1
eps=0.1
N=500
burn_in=50
t0=100
al=0.15
sd=2.5

## how many simulations to run
r <- 8
## Simulate the experiments
start <- Sys.time()
sims <- mclapply(X=1:r, 
                 function(X){sim(N=N,mean=0,sd=sd,
                                 burn_in=burn_in,A=A,
                                 sigma=sigma,t0=t0,eps=eps,al=al)},
                 mc.cores = 1
)
end <- Sys.time()

## collect experiment information into one dataframe
exps <- data.frame()
for(i in 1:r){
  
  ## pronzato
  pronzato <- sims[[i]]$pronzato
  pronzato$R <- NA
  pronzato$method <- "pronzato"
  pronzato$rep <- i
  
  ## d_opt
  d_opt <- sims[[i]]$d_opt
  d_opt$R <- NA
  d_opt$method <- "d_opt"
  d_opt$rep <- i
  
  ## e_greedy
  e_greedy <- sims[[i]]$e_greedy
  e_greedy$R <- NA
  e_greedy$method <- "e_greedy"
  e_greedy$rep <- i
  
  ## greedy
  greedy <- sims[[i]]$greedy
  greedy$R <- NA
  greedy$method <- "greedy"
  greedy$rep <- i
  
  ## greedy first
  greedy_first <- sims[[i]]$greedy_first
  greedy_first$method <- "greedy_first"
  greedy_first$rep <- i
  
  ## IDS_freq
  IDS_freq <- sims[[i]]$IDS_freq
  IDS_freq$R <- NA
  IDS_freq$method <- "IDS_freq"
  IDS_freq$rep <- i
  
  ## attach these datasets to the bigger one
  exps <- rbind(exps,pronzato,d_opt,greedy,e_greedy,
                greedy_first,IDS_freq)
  
  
}


exps %>% 
  filter(sub>burn_in) %>%
  group_by(method) %>%
  summarise(mean(regret))

exps %>% 
  filter(sub>t0 & method=="greedy_first") %>%
  summarise(mean(R))

convergence <- exps %>% 
  filter(sub>burn_in) %>%
  group_by(method,sub) %>%
  summarise(norm=mean(norm))


regret <- exps %>%
  filter(sub>burn_in) %>%
  group_by(method,rep) %>%
  mutate(cum_regret=cumsum(regret)) %>%
  group_by(method,sub) %>%
  summarise(cum_regret=mean(cum_regret))

ggplot(data=convergence) +
  geom_line(aes(x=sub,y=norm,color=method))

ggplot(data=regret) +
  geom_line(aes(x=sub,y=cum_regret,color=method))


theta1 <- c(3.25,0.01,0.0)
theta2 <- c(2.25,-0.25,0.01)
theta3 <- c(3.75,0.05,-0.20)
theta4 <- c(3.15,-0.40,-0.15)
theta5 <- c(2.5,0.05,-0.05)
theta <- c(theta1,theta2,theta3,theta4,theta5)

theta <- mvrnorm(n=1,mu=rep(0,15),Sigma=diag(15))


no_sigma <- gen_data(N=1000,-1,1,A,theta,sigma=0)
ggplot(data=no_sigma) +
  geom_line(aes(x=X,y=Y,color=as.factor(A)))

