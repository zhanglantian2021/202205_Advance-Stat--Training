
library(brms)
library(ggplot2)
library(ape)
library(caper)
library(performance)


#################################################
#### load trait data #####
#################################################
getwd()

setwd("/Users/kyletomlinson/Dropbox/Teaching/XTBG advanced stats 2021/Lectures/Lesson 10 The animal model")

# database with the corrected data
bark <- read.csv("bark/bark_nutrition.csv", header = TRUE, sep = ",", row.names = NULL)
summary(bark)

bark$Species <- as.factor(bark$Species)
bark$Family <- as.factor(bark$Family)
bark$Month <- as.factor(bark$Month)
bark$Function <- as.factor(bark$Function)

summary(bark)
dim(bark)

length(unique(bark$Species)) #53 species

#################################################
#### load phylogeny data #####
#################################################


phylo <- read.tree("bark/scenario.2_run.1.tre")
plot(phylo)
summary(phylo)  #31 species

# there are no polytomies
is.binary.tree(phylo) #TRUE

# check tree is ultrametric
is.ultrametric(phylo) # TRUE


#################################################
### reduce data types to set of common species ##
#################################################

# B. reduce both data types to the minimum set of common species

#try comparative data
combine<-comparative.data(phylo,bark,Species,vcv=TRUE,na.omit=F)
### doesnt work because the function cannot handle repeats

#different route
#first reduce data to list of common taxa with tree
n1<-unique(phylo$tip.label)
str(n1)
bark2 <- bark[bark$Species %in% n1,]

dim(bark2) 
length(unique(bark2$Species)) #only 27, less than the tree

#so, second reduce tree to match common taxa with data
n2<-unique(bark2$Species)
phylo2<-keep.tip(phylo,phylo$tip.label[match(n2, phylo$tip.label)])
summary(phylo2)

#################################################
### make covariance matrix ##
#################################################

#make the covariance matrix from the phylogeny
phylo_cor <- ape::vcv.phylo(phylo2)
summary(phylo_cor)

#make an additional factor to reference the species in the phylogeny
bark2$phylo2 <- bark2$Species

summary(bark2)

#################################################
#### response data distribution #####
#################################################
#we will concentrate on D

#Distribution
par(mfrow=c(1,3))
plot(density((bark2$D)))
plot(density(sqrt(bark2$D)))
plot(density(log(bark2$D)))

bark2 $Dlog=log(bark2 $D)

#bark2 $int=interaction(bark2 $Function,bark2$Month)

summary(bark2)
#################################################
#### MODELS #####
#################################################

#Full model: Crude protein available in bark
model_Flog <-brm(Dlog~Function+(1|phylo2)+(1|Species), 
                  cov_ranef = list(phylo2 = phylo_cor),
                  data = bark2, family = gaussian(), sample_prior = TRUE,
                  iter = 10000,warmup = 5000, chains = 4, cores = 4,
                  thin = 10, save_all_pars = TRUE,seed=T,
                  control=list(adapt_delta=0.99))

# model_Flog <-brm(Dlog~Function+(1|phylo2)+(1|Species), 
#                  cov_ranef = list(phylo2 = phylo_cor),
#                  data = bark2, family = gaussian(), sample_prior = TRUE,
#                  iter = 10000,warmup = 5000, chains = 4, cores = 2,
#                  thin = 10, save_all_pars = TRUE,seed=T,
#                  control=list(adapt_delta=0.8, max_treedepth=10))

saveRDS(model_Flog, "bark/Dlog_Int_Phyl_Spec.rds")

model_Flog <- readRDS("bark/Dlog_Int_Phyl_Spec.rds")


#use summary() to check convergence of the model and the parameter estimates
summary(model_Flog)

#now look at the posterior distributions generated for each estimate and the chain value (this should form nice horizontal bands called "hairy caterpillar")
plot(model_Flog, N = 6, ask = FALSE)

#plot the posterior fixed effect estimates with credible intervals over the data
plot(conditional_effects(model_Flog), points = TRUE)

#Posterior estimates and confidence intervals: tests of significance
mcmc_plot(model_Flog, pars = c("^b_", "^sd_"))

#estimate proportional variation explained by phylogenetic signal
hyp <- "(sd_phylo2__Intercept^2)/ (sd_phylo2__Intercept^2 + sd_Species__Intercept^2 + sigma^2) = 0"
(hyp <- hypothesis(model_Flog, hyp, class = NULL))
plot(hyp)

#R-squared: variance explained by model
bayes_R2(model_Flog)


#output table including conditional and marginal R2
devtools::install_github("strengejacke/strengejacke")
library(sjPlot)
tab_model(model_Flog, show.intercept = T, show.r2 = TRUE, 
          transform = NULL,show.re.var=T)


#################################################

#what follows are some tests evaluating the performance of the analysis

#Posterior predictive distribution: good data fitting
pp_check(model_Flog, nsamples =1000)
#pp_check shows the distributions of many replicated data sets drawn from the posterior predictive distribution (thin light curves) compared with the empirical distribution of the observed outcome (the thick dark curve).


#Posterior predictive intervals: good data fitting
pp_check(model_Flog, type = "intervals_grouped", group = "Function")


#Posterior predictive fit: good data fitting
pp_check(model_Flog, type = "scatter_avg_grouped", group = "Function") + 
  geom_abline(intercept = 0, slope = 1 , color = "red", lty = 2)
####  Scatterplots (and similar visualizations) of the data y vs. individual simulated datasets (rows) in yrep, or vs. the average value of the distributions of each data point (columns) in yrep.


####More check if needed to improve models (more can be found)
#Effective number of independent simulation draws: good model running
mcmc_plot(model_Flog, type = "neff")
#Ratios of effective sample size to total sample size as either points or a histogram. Values are colored using different shades (lighter is better). The chosen thresholds are somewhat arbitrary, but can be useful guidelines in practice.



#Autocorrelation
library(rstan)
stan_ac(model_Flog$fit)

#Interclass Correlation Coefficient: variance partition coefficient
library(performance)
library(dplyr)
icc(model_Flog) %>% print(prob = .95, digits = 3)
##proportion of non-fixed variation explained by the random effects (phylogeny and species)
## ICC (intra-class correlation coefficients) (????ЧӦ??ռ?ķ????ڳ??̶?ЧӦ֮?????ܷ???????ռ????)

#now just estimate variance explained by phylogeny
# icc(model_Flog,ppd=T,re.form = ~(1| phylo2),robust=T)
icc(model_Flog)

#################################################

#Reduced model without phylogeny

#Reduced model without random phylogeny
model_Flog2 <-brm(Dlog~Function+(1|Species),
                  data = bark2, family = gaussian(), sample_prior = TRUE,
                  iter = 10000,warmup = 5000, chains = 4, cores = 2,
                  thin = 10, save_all_pars = TRUE,seed=T, 
                  control=list(adapt_delta=0.8, max_treedepth=10))

saveRDS(model_Flog2, "bark/Dlog_Int_Spec.rds")
model_Flog2 <- readRDS("bark/Dlog_Int_Spec.rds")

summary(model_Flog2)


#plot the posterior fixed effect estimates with credible intervals over the data
plot(conditional_effects(model_Flog), points = TRUE)

#Posterior estimates and confidence intervals: tests of significance
mcmc_plot(model_Flog, pars = c("^b_", "^sd_"))




####Model selection
#Leave-one-out cross-validation: compare models performance
loo(model_Flog,model_Flog2,reloo=F)
##reloo=T takes a lot of time!!
##Models with lower LOO are considered to be better models, particularly if the difference in LOO between two models (ΔLOO) is more than twice its standard error.

