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
#source("d_opt.R")
source("greedy.R")

## load packages
library(parallel)
library(tidyverse)
library(MASS)


## sim
## Purpose: simulate all methods on the same dataset
## param N: total sample size
## param p: number of features in the context (not including intercept)
## param K: number of arms
## param sd_X: standard deviation of the context
## param sd_Y: scaling factor sd(Y)=(sd_Y)*sqrt(p+1)
## param t0: time for greedy first to stay with 
##           greedy before checking
## param eps: P(explore) = eps in greedy first if it switches
## param al: value for weight put on penalty in pronzato method
## return output: list of trial information from each method
sim <- function(N,p,K,sd_X,sd_Y,
                t0,eps,al){
  
  ## standard deviation of Y
  sd_Y <- sd_Y*sqrt(p)
  ## burn in period is (# of parameters in full model)*3
  burn_in <- (p+1)*K*3
  ## vector of treatments
  A <- 1:K
  
  ## randomly sample the mean parameters
  theta <- rnorm((p+1)*K,0,1)
  
  ## generate training set
  train_set <- gen_data(N=N,p=p,sd_X=sd_X,A=1:K,sd_Y=sd_Y,theta=theta)
  
  ## run pronzato experiment
  pronzato_sim <- try(pronzato(train_set=train_set,burn_in=burn_in,
                      A=A,theta=theta,sd_Y=sd_Y,al=al))
  
  ## run e_greedy experiment
  e_greedy_sim <- try(e_greedy(train_set=train_set,burn_in=burn_in,
                      sd_Y=sd_Y,A=A,theta=theta,eps=eps))
  
  
  ## run greedy experiment
  greedy_sim <- try(greedy(train_set=train_set,burn_in=burn_in,
                           sd_Y=sd_Y,A=A,theta=theta))
  
  ## run GF experiment
  greedy_first_sim <- try(greedy_first(train_set=train_set,burn_in=burn_in,
                           sd_Y=sd_Y,A=A,theta=theta,t0=t0,eps=eps))
  
  ## run IDS experiment
  IDS_freq_sim <- try(IDS_freq(train_set=train_set,burn_in=burn_in,
                               sd_Y=sd_Y,A=A,theta=theta))
  
  ## save output in a list
  output <- list(pronzato=pronzato_sim,
                 e_greedy=e_greedy_sim,
                 greedy=greedy_sim,
                 greedy_first=greedy_first_sim,
                 IDS_freq=IDS_freq_sim)
  
  return(output)
  
}

## run the simulations

## parameters
N=500
p=2
K=5
sd_X=1
sd_Y=0.5
t0=(p+1)*K*3+100
eps=0.1
al=0.05

## how many simulations to run
r <- 3
## Simulate the experiments
start <- Sys.time()
sims <- mclapply(X=1:r, 
                 function(X){sim(N=N,p=p,K=K,sd_X=sd_X,sd_Y=sd_Y,
                                 t0=t0,eps=eps,al=al)},
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
  exps <- rbind(exps,pronzato,greedy,e_greedy,
                greedy_first,IDS_freq)
  
  
}

burn_in=(p+1)*K*3


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

