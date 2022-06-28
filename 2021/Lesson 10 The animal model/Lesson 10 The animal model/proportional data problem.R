
#################################################
#### DATASET #####
#################################################
setwd("G:/21.02.2020/Final article 2/Simulation/Model debarking")

# database with the corrected data
debark <- read.csv("debarking data - Kyle course.csv", header = TRUE, sep = ";", row.names = NULL)


# loading phylogenetic matrix "phylo_cor"
load("G:/21.02.2020/Final article 2/Simulation/Phylogeny/phylo_cor.Rdata") #phylo_cor from phylogeny


#################################################
#### DATA EXPLORATION #####
#################################################
head(debark$Ring.P, n=31)
head(debark$Ring.S, n=31)
head(debark$Prop.P, n=31)
head(debark$Prop.S, n=31)

boxplot(debark$Ring.P~debark$Strategy)
boxplot(debark$Ring.S~debark$Strategy)
boxplot(debark$Prop.P~debark$Strategy)
boxplot(debark$Prop.S~debark$Strategy)

plot(density(debark$Ring.S))
plot(density(debark$Ring.P))
plot(density(debark$Prop.S))
plot(density(debark$Prop.P))

#################################################
#### MODELS #####
#################################################
library(brms)

#Full model: proportion of debarking for a squirrel
debark$Prop.S.Beta=(debark$Prop.S)/100

model_Flog <-brm(Prop.S.Beta~Strategy+(1|Phylo), 
                 cov_ranef = list(Phylo = phylo_cor),
                 data = debark, family = Beta(), sample_prior = TRUE,
                 iter = 20000,warmup = 10000, chains = 4, cores = 4,
                 thin = 10, save_all_pars = TRUE, seed = T,
                 control=list(adapt_delta=0.95, max_treedepth=15))

saveRDS(model_Flog, "Prop.S.beta_Phyl.rds")

#Reduced model: proportion of debarking for a squirrel
model_Flog <-brm(Prop.S.Beta~Strategy,
                 data = debark, family = Beta(), sample_prior = TRUE,
                 iter = 20000,warmup = 10000, chains = 4, cores = 4,
                 thin = 10, save_all_pars = TRUE, seed = T,
                 control=list(adapt_delta=0.95, max_treedepth=15))

saveRDS(model_Flog, "Prop.S.beta.rds")

#################################################
#### MODELS ASSUMPTIONS AND RESULTS #####
#################################################

#Models 
model_Flog1 <- readRDS("Prop.S.beta_Phyl.rds")
model_Flog2 <- readRDS("Prop.S.beta.rds")

#Website to explore models characteristics
#launch_shinystan(model_Flog1)

#Summary: good convergence of the model 
#(effective sample size measures, good chain convergence,...)
summary(model_Flog1)
summary(model_Flog2)

#Posterior predictive distribution: good data fitting
pp_check(model_Flog1, nsamples =100)
pp_check(model_Flog2, nsamples =100)

#Posterior predictive intervals: good data fitting
pp_check(model_Flog1, type = "intervals_grouped", group = "Strategy")
pp_check(model_Flog2, type = "intervals_grouped", group = "Strategy")

#Posterior predictive fit: good data fitting
library(ggplot2)
pp_check(model_Flog1, type = "scatter_avg_grouped", group = "Strategy") + 
  geom_abline(intercept = 0, slope = 1 , color = "red", lty = 2)
pp_check(model_Flog2, type = "scatter_avg_grouped", group = "Strategy") + 
  geom_abline(intercept = 0, slope = 1 , color = "red", lty = 2)

####More check if needed to improve models (more can be found)
#Effective number of independent simulation draws: good model running
mcmc_plot(model_Flog1, type = "neff")
mcmc_plot(model_Flog2, type = "neff")

#No-U-turn: efficient sampling from the posterior
mcmc_plot(model_Flog1, type = "nuts_acceptance")
mcmc_plot(model_Flog2, type = "nuts_acceptance")

#Autocorrelation
library(rstan)
stan_ac(model_Flog1$fit)

####Model selection
#Leave-one-out cross-validation: compare models performance
loo(model_Flog1,model_Flog2)

#R-squared: variance explained by model
bayes_R2(model_Flog1)
bayes_R2(model_Flog2)

#Interclass Correlation Coefficient: variance partition coefficient
library(performance)
library(dplyr)
icc(model_Flog1, ppd=TRUE) %>% print(prob = .95, digits = 3)
icc(model_Flog2, ppd=TRUE) %>% print(prob = .95, digits = 3)

#Boxplot fitted values and original data: good to compare model and data
yrep<-as.data.frame(fitted(model_Flog1))
yrep$Strategy<-debark$Strategy

library(ggplot2)
ggplot()+
  geom_boxplot(aes(x = debark$Strategy, y=debark$Prop.S.Beta))+
  geom_point(aes(x = debark$Strategy, y=debark$Prop.S.Beta))+
  geom_boxplot(aes(x = yrep$Strategy, y=yrep$Estimate),col=2)+
  geom_point(aes(x =yrep$Strategy, y=yrep$Estimate), col=2) 

####Results
#Posterior estimates and confidence intervals: test of significance
mcmc_plot(model_Flog1, pars = c("^b_", "^sd_"))
mcmc_plot(model_Flog2, pars = c("^b_", "^sd_"))

#Confidence intervals
library(broom)
tidyMCMC(model_Flog1, conf.int = TRUE, conf.method = "HPDinterval")
tidyMCMC(model_Flog2, conf.int = TRUE, conf.method = "HPDinterval")

#Posterior estimates and uncertainty interval of the response
plot(conditional_effects(model_Flog1),points=T)
plot(conditional_effects(model_Flog2),points=T)

#Phylogenetic signal lambda: indicate by estimate

#Case for gaussian distribution when phylogeny and repeated measures are 
#random effects
hyp <- paste(
  "sd_Phylo__Intercept^2 /", 
  "(sd_Phylo__Intercept^2 + sd_Species__Intercept^2+sigma^2) = 0"
)
(hyp1 <- hypothesis(model_Flog1, hyp, class = NULL))
plot(hyp1)

#for other distribution: I found one paper using bernoulli and they change
#sigma by ??2/3. But no reference or explanation are provided.