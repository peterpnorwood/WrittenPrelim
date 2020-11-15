## ----------------------------------------------------------------- ##
## analyze_sims.R -------------------------------------------------- ##
## Author: Peter Norwood, NC State University ---------------------- ##
## Purpose: create tables and figures for the sim study ------------ ##
## ----------------------------------------------------------------- ##

## load libraries
library(tidyverse)
library(xtable)

## load functions
setwd("~/Research/Written Prelim/Simiulation Data")

## load experiment data
load("exps_K=2_p=2.RData")
exps22 <- exps; exps <- NULL;
load("exps_K=2_p=8.RData")
exps28 <- exps; exps <- NULL;
load("exps_K=8_p=2.RData")
exps82 <- exps; exps <- NULL;
load("exps_K=8_p=8.RData")
exps88 <- exps; exps <- NULL;
## put them all together
exps <- rbind(exps22[,3:14],exps28[,9:20],exps82[,3:14],exps88[,9:20])

## create "scenario" variable
exps <- exps %>%
        mutate(scenario=paste0("K=",K," p=",p),
               burn_in=3*(p+1)*K) %>%
        filter(rep<=5000)

## analyze cumulative regret 
regret <- exps %>% filter(sub>burn_in) %>%
  group_by(scenario,method,rep) %>%
  mutate(cum_regret=cumsum(regret)) %>%
  group_by(scenario,method,sub) %>%
  summarise(mean_regret=mean(cum_regret),
            median_regret=median(cum_regret),
            se_regret=sd(cum_regret)/sqrt(5000))

save(regret,file="regret.RData")

## analyze convergence
theta <- exps %>% filter(sub>burn_in) %>%
  group_by(scenario,method,sub) %>%
  mutate(rel_efficiency=simple_norm/norm) %>%
  summarise(mean_norm=mean(norm),
            median_norm=median(norm),
            se_norm=sd(norm)/sqrt(5000),
            mean_simple_norm=mean(simple_norm)) %>%
  mutate(rel_eff=mean_simple_norm/mean_norm)

save(theta,file="theta.RData")

## ----------------------------------------------------------------- ##
## ----------------------------------------------------------------- ##

## post experiment data
load("post_K=2_p=2.RData")
post22 <- post; post <- NULL;
load("post_K=2_p=8.RData")
post28 <- post; post <- NULL;
load("post_K=8_p=2.RData")
post82 <- post; post <- NULL;
load("post_K=8_p=8.RData")
post88 <- post; post <- NULL;
## put them all together
post <- rbind(post22,post28,post82,post88)

## create scenario column
## create "scenario" variable
post <- post %>%
  mutate(scenario=paste0("K=",K," p=",p)) %>%
  filter(rep<5000)

## post trial correctimal number correct
correct_post <- post %>% 
      group_by(method,scenario,rep) %>%
      summarise(prop_correct=mean(correct)) %>%
      group_by(scenario,method) %>%
      summarise(mean_correct=mean(prop_correct),
                median_correct=median(prop_correct),
                se_correct=sd(prop_correct)/sqrt(5000))
View(correct_post)

save(correct_post,file="correct_post.RData")

regret_post <- post %>% 
  group_by(method,scenario,rep) %>%
  summarise(regret=sum(regret)) %>%
  group_by(scenario,method) %>%
  summarise(mean_cum_regret=mean(regret),
            median_cum_regret=median(regret),
            se_cum_regret=sd(regret)/sqrt(5000))

View(regret_post)

save(regret_post,file="correct_post.RData")


## ----------------------------------------------------------------- ##
## ----------------------------------------------------------------- ##

## calcualte ratios
regret <- regret %>% filter(sub==1000)
correct2 <- correct_post %>% filter(method!="simple")
norm2 <- theta %>% filter(sub==1000)
regret_ratio <- cbind(regret,correct2,norm2)
regret_ratio <- regret_ratio %>%
  mutate(ratio1=mean_regret/mean_correct,
         ratio2=mean_regret/rel_eff) %>%
  select(scenario,method,ratio1,ratio2)
View(regret_ratio)
