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
source("IDS_bayes.R")
source("greedy_first.R")
source("e_greedy.R")
source("greedy.R")
source("post_experiment.R")

## load packages
library(parallel)
library(tidyverse)
library(MASS)
library(foreach)
library(doParallel)

#setup parallel backend to use many processors
cores=detectCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl)
clusterExport(cl,c('mean_outcome','mean_outcome1',
                   'make_design','gen_data'))

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
## param N_post: sample to take out of the experiment
## param M: samples to take out of posterior distribution for IDS
## return output: list of trial information from each method
##                and out of trial information
sim <- function(N,p,K,sd_X,sd_Y,
                t0,eps,al,N_post,M){
  
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
  
  ## run IDS experiment
  IDS_bayes_sim <- try(IDS_bayes(train_set=train_set,burn_in=burn_in,
                                 sd_Y=sd_Y,A=A,theta=theta,M=M))
  
  
  ## out of experiment information
  
  ## out of experiment data
  post_data <- gen_data(N=N_post,p=p,sd_X=sd_X,A=1:K,sd_Y=sd_Y,theta=theta)
  
  post_pronzato <- try(post_experiment(exp_data=pronzato_sim,
                                       post_data=post_data,
                                       p=p,K=K,theta=theta))
  
  post_e_greedy <- try(post_experiment(exp_data=e_greedy_sim,
                                       post_data=post_data,
                                       p=p,K=K,theta=theta))
  
  post_greedy <- try(post_experiment(exp_data=greedy_sim,
                                       post_data=post_data,
                                p=p,K=K,theta=theta))
  
  post_greedy_first <- try(post_experiment(exp_data=greedy_first_sim,
                                     post_data=post_data,
                                     p=p,K=K,theta=theta))
  
  post_IDS_freq <- try(post_experiment(exp_data=IDS_freq_sim,
                                     post_data=post_data,
                                     p=p,K=K,theta=theta))
  
  post_IDS_bayes <- try(post_experiment(exp_data=IDS_bayes_sim,
                                       post_data=post_data,
                                       p=p,K=K,theta=theta))
  
  post_simple <- try(post_experiment(exp_data=train_set,
                                        post_data=post_data,
                                        p=p,K=K,theta=theta))
  
  
  
  ## save output in a list
  exp_output <- list(train_set=train_set,
                     pronzato=pronzato_sim,
                     e_greedy=e_greedy_sim,
                     greedy=greedy_sim,
                     greedy_first=greedy_first_sim,
                     IDS_freq=IDS_freq_sim,
                     IDS_bayes=IDS_bayes_sim)
  
  post_output <- list(pronzato=post_pronzato,
                      e_greedy=post_greedy,
                      greedy=post_greedy,
                      greedy_first=post_greedy_first,
                      IDS_freq=post_IDS_freq,
                      IDS_bayes=post_IDS_bayes,
                      simple=post_simple)
  
  output <- list(exp=exp_output,post=post_output,theta=theta)
  
  return(output)
  
}
