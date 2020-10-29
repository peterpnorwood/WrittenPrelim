## ----------------------------------------------------------------- ##
## run_sims.R ------------------------------------------------------ ##
## Author: Peter Norwood, NC State University ---------------------- ##
## Purpose: run a simulation experiment with certain inputs -------- ##
## ----------------------------------------------------------------- ##

## load functions
setwd("~/Research/Written Prelim/WrittenPrelim")
source("sims.R")
library(parallel)

## run the simulations

## parameters
p=2
K=8

N=1000
sd_X=1
sd_Y=0.5
t0=(p+1)*K*3+100
eps=0.05
al=0.01
N_post=500
M=100

## how many simulations to run
r <- 100
## Simulate the experiments
start <- Sys.time()
sims <- mclapply(X=1:r, 
                 function(X){sim(N=N,p=p,K=K,sd_X=sd_X,sd_Y=sd_Y,
                                 t0=t0,eps=eps,al=al,N_post=N_post,M=M)},
                 mc.cores = 1
)
end <- Sys.time()


## check bad sims
good_sets <- c()
tick=1
for(i in 1:r){
  check <- c()
  for(j in 1:6){
    check[j] <- is.data.frame(sims[[i]]$exp[[j]])
  }
  if(sum(check)==6){
    good_sets[tick]=i
    tick=tick+1
  }else{
    ## no nothing
  }
}

## collect experiment information into one dataframe
exps <- data.frame()
post <- data.frame()
tick <- 1
for(i in good_sets){
  
  ## check convergence for simple randomization
  exp_simple <- sims[[i]]$exp$train_set
  fit <- lm(Y~-1+as.factor(A)+as.factor(A):.-as.factor(A):A,
            data=exp_simple)
  
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
  
  ## simple norm
  simple_norm <- norm(matrix(sims[[i]]$theta-theta_hat),type="F")
  
  ## collect experiment information
  ## pronzato
  exp_pronzato <- sims[[i]]$exp$pronzato
  exp_pronzato$R <- NA
  exp_pronzato$method <- "pronzato"
  exp_pronzato$rep <- tick
  exp_pronzato$simple_norm <- simple_norm
  
  post_pronzato <- sims[[i]]$post$pronzato
  post_pronzato$method <- "pronzato"
  post_pronzato$rep <- tick
  
  ## e_greedy
  exp_e_greedy <- sims[[i]]$exp$e_greedy
  exp_e_greedy$R <- NA
  exp_e_greedy$method <- "e_greedy"
  exp_e_greedy$rep <- tick
  exp_e_greedy$simple_norm <- simple_norm
  
  post_e_greedy <- sims[[i]]$post$e_greedy
  post_e_greedy$method <- "e_greedy"
  post_e_greedy$rep <- tick
  
  ## greedy
  exp_greedy <- sims[[i]]$exp$greedy
  exp_greedy$R <- NA
  exp_greedy$method <- "greedy"
  exp_greedy$rep <- tick
  exp_greedy$simple_norm <- simple_norm
  
  post_greedy <- sims[[i]]$post$greedy
  post_greedy$method <- "greedy"
  post_greedy$rep <- tick
  
  ## greedy first
  exp_greedy_first <- sims[[i]]$exp$greedy_first
  exp_greedy_first$method <- "greedy_first"
  exp_greedy_first$rep <- tick
  exp_greedy_first$simple_norm <- simple_norm
  
  post_greedy_first <- sims[[i]]$post$greedy_first
  post_greedy_first$method <- "greedy_first"
  post_greedy_first$rep <- tick
  
  ## IDS_freq
  exp_IDS_freq <- sims[[i]]$exp$IDS_freq
  exp_IDS_freq$R <- NA
  exp_IDS_freq$method <- "IDS_freq"
  exp_IDS_freq$rep <- tick
  exp_IDS_freq$simple_norm <- simple_norm
  
  post_IDS_freq <- sims[[i]]$post$IDS_freq
  post_IDS_freq$method <- "IDS_freq"
  post_IDS_freq$rep <- tick
  
  ## IDS_bayes
  exp_IDS_bayes <- sims[[i]]$exp$IDS_bayes
  exp_IDS_bayes$R <- NA
  exp_IDS_bayes$method <- "IDS_bayes"
  exp_IDS_bayes$rep <- tick
  exp_IDS_bayes$simple_norm <- simple_norm
  
  post_IDS_bayes <- sims[[i]]$post$IDS_bayes
  post_IDS_bayes$method <- "IDS_bayes"
  post_IDS_bayes$rep <- tick
  
  ## simple
  post_simple <- sims[[i]]$post$simple
  post_simple$method <- "simple"
  post_simple$rep <- tick
  
  ## create out of trial information
  ## attach these datasets to the bigger one
  exps <- rbind(exps,exp_pronzato,exp_greedy,exp_e_greedy,
                exp_greedy_first,exp_IDS_freq,exp_IDS_bayes)
  
  post <- rbind(post,post_pronzato,post_greedy,post_e_greedy,
                post_greedy_first,post_IDS_freq,post_IDS_bayes,post_simple)
  
  
  ## update the tick
  tick=tick+1
  
}

## save trial information
exps$K<- K
exps$p <- p
post$K <- K
post$p <- p

save(exps,file=paste0("exps_K=",K,"_p=",p,".RData"))
save(post,file=paste0("post_K=",K,"_p=",p,".RData"))

# burn_in=(p+1)*K*3
# 
# exps %>%
#   filter(sub>burn_in) %>%
#   group_by(method) %>%
#   summarise(mean(regret))
# 
# exps %>%
#   filter(sub>t0 & method=="greedy_first") %>%
#   summarise(mean(R))
# 
# convergence <- exps %>%
#   filter(sub>burn_in) %>%
#   group_by(method,sub) %>%
#   summarise(norm=mean(norm))
# 
# ggplot(exps %>% filter(sub==N)) +
#   geom_boxplot(aes(x=method,y=norm/simple_norm))
# 
# regret <- exps %>%
#   filter(sub>burn_in) %>%
#   group_by(method,rep) %>%
#   mutate(cum_regret=cumsum(regret)) %>%
#   group_by(method,sub) %>%
#   summarise(cum_regret=mean(cum_regret))
# 
# opt <- exps %>%
#   filter(sub>burn_in) %>%
#   mutate(opt=ifelse(regret==0,1,0)) %>%
#   group_by(method,rep) %>%
#   mutate(cum_opt=cumsum(opt),
#          cum_prop=cum_opt/sub) %>%
#   group_by(method,sub) %>%
#   summarise(cum_opt=mean(cum_prop))
# 
# ggplot(data=opt) +
#   geom_line(aes(x=sub,y=cum_opt,color=method))
# 
# ggplot(data=convergence) +
#   geom_line(aes(x=sub,y=norm,color=method))
# 
# ggplot(data=regret) +
#   geom_line(aes(x=sub,y=cum_regret,color=method))
# 
# 
# post %>%
#   group_by(method) %>%
#   summarise(mean(correct))
# 
# post %>%
#   group_by(method) %>%
#   summarise(mean(regret))
