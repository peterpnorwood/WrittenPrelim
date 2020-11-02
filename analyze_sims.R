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
            se_regret=sd(cum_regret)/5000)

View(regret %>% filter(sub==1000)%>%
select(scenario,method,mean_regret,median_regret,se_regret))

ggplot(data=regret) +
  geom_line(aes(x=sub,y=mean_regret,color=method)) +
  facet_wrap(vars(scenario),scales="free") 


## analyze proportion of correctimal intervention
correct <- exps %>% filter(sub>burn_in) %>%
  group_by(scenario,method,rep) %>%
  mutate(correct=ifelse(regret==0,1,0),
         cum_correct=cumsum(correct),
         prop_correct=cum_correct/sub) %>%
  group_by(scenario,method,sub) %>%
  summarise(mean_correct=mean(prop_correct),
            median_correct=median(prop_correct),
            se_correct=sd(prop_correct)/5000)

View(correct %>% filter(sub==1000)%>%
  select(scenario,method,mean_correct,median_correct,se_correct))

ggplot(data=correct) +
  geom_line(aes(x=sub,y=mean_correct,color=method)) +
  facet_wrap(vars(scenario),scales="free") 

## analyze convergence
theta <- exps %>% filter(sub>burn_in) %>%
  group_by(scenario,method,sub) %>%
  mutate(rel_efficiency=simple_norm/norm) %>%
  summarise(mean_norm=mean(norm),
            median_norm=median(norm),
            se_norm=sd(norm)/5000,
            mean_simple_norm=mean(simple_norm)) %>%
  mutate(rel_eff=mean_simple_norm/mean_norm)

View(theta %>% filter(sub==1000)%>%
  select(scenario,method,mean_norm,median_norm,se_norm,rel_eff))

ggplot(data=theta) +
  geom_line(aes(x=sub,y=mean_norm,color=method)) +
  facet_wrap(vars(scenario),scales="free") 

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
                se_correct=sd(prop_correct)/5000)

regret_post <- post %>% 
  group_by(method,scenario,rep) %>%
  summarise(regret=sum(regret)) %>%
  group_by(scenario,method) %>%
  summarise(mean_cum_regret=mean(regret),
            median_cum_regret=median(regret),
            sd_cum_regret=sd(regret)/5000)
View(regret_post)
