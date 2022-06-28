setwd("E:/academic_resources/Advance Stat Training/2022/Lesson 10 The animal model")


install.packages("brms")
library(brms)
library(ggplot2)
library(ape)
library(caper)
install.packages("sjPlot")
library(sjPlot)
install.packages("MLmetrics")
library(MLmetrics)
library(MuMIn)
library(rstan)


#devtools::install_version('rstan', version = '2.19.3', repos = "http://cran.us.r-project.org")

#################################################
#### load trait data #####
#################################################

# database with the corrected data
bark <- read.csv("bark/bark_nutrition.csv", header = TRUE, sep = ",", row.names = NULL)
summary(bark)

bark$Species <- as.factor(bark$Species)
bark$Family <- as.factor(bark$Family)
bark$Month <- as.factor(bark$Month)
bark$Function <- as.factor(bark$Function)

summary(bark) #note that there are multiple individuals per species!!! (not possible to have this in many comparative packages)
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
#if your tree has polytomies (i.e. "FALSE") then you need to trick it into having fake binary splits using multi2di()

# check tree is ultrametric
is.ultrametric(phylo) # TRUE


#################################################
### reduce data types to set of common species ##
#################################################

# B. reduce both data types to the minimum set of common species

#try comparative data
combine<-comparative.data(phylo,bark,Species,vcv=TRUE,na.omit=F)
#doesnt work because the function cannot handle species repeats

#different route
#first reduce data to list of common taxa with tree
n1<-unique(phylo$tip.label)
str(n1)
bark2 <- bark[bark$Species %in% n1,]
dim(bark)
dim(bark2) 
length(unique(bark2$Species)) #only 27, less than the tree

#so, second reduce tree to match common taxa with data
n2<-unique(bark2$Species)
phylo2<-keep.tip(phylo,phylo$tip.label[match(n2, phylo$tip.label)])
summary(phylo2) #now 27 tips, which match the ones you made above


#################################################
### make covariance matrix ##
#################################################

#make the covariance matrix from the phylogeny
phylo_cor <- ape::vcv.phylo(phylo2)
summary(phylo_cor)
dim(phylo_cor)

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

#take a look at what it says about the function
?brm

#lets run the model using 
set.seed(1)
model_Flog <-brm(Dlog~Function+(1|phylo2)+(1|Species), 
                 cov_ranef = list(phylo2 = phylo_cor),
                 data = bark2, family = gaussian(link="identity"), sample_prior = TRUE,
                 iter = 10000,warmup = 5000, chains = 4, cores = 4,
                 thin = 10, save_all_pars = TRUE,seed=T,
                 control=list(adapt_delta=0.99))

#learn more about brms
#vignette("brms_overview")
#vignette("brms_multilevel")

#https://stackoverflow.com/questions/63025052/brm-model-compiling-but-returning-model-object
model_Flog

#NOTE: you can save the output (a good idea when your model has a lot of data and takes a long time to run)
#saveRDS(model_Flog, "bark/Dlog_Int_Phyl_Spec.rds")

#NOTE: the next command allows you to read in the saved output
#model_Flog <- readRDS("bark/Dlog_Int_Phyl_Spec.rds")


#use summary() to check convergence of the model (Rhat and ESS) and the parameter estimates (coefficients and betas)
summary(model_Flog)

#now look at the posterior distributions generated for each estimate and the chain value (this should form nice horizontal bands called "hairy caterpillars"; if not flat then sampling possibly insufficient or model underspecified)
quartz(); plot(model_Flog, N = 6, ask = FALSE)

#plot the posterior fixed effect estimates with credible intervals over the data
plot(conditional_effects(model_Flog), points = TRUE)

#Posterior estimates and confidence intervals: tests of significance
mcmc_plot(model_Flog, pars = c("^b_", "^sd_"))
#output shows 99% (wider) and 95% (narrower) credible intervals

#estimate proportional variation explained by phylogenetic signal
hyp <- "(sd_phylo2__Intercept^2)/ (sd_phylo2__Intercept^2 + sd_Species__Intercept^2 + sigma^2) = 0"
(hyp <- hypothesis(model_Flog, hyp, class = NULL))


#R-squared: variance explained by model
bayes_R2(model_Flog)


#output table including conditional and marginal R2 
library(sjPlot)
tab_model(model_Flog, show.intercept = T, show.r2 = TRUE, 
          transform = NULL,show.re.var=T)


#################################################

#what follows are some tests evaluating the performance of the analysis
#library(performance)
#Posterior predictive distribution: good data fitting
quartz()
brms::pp_check(model_Flog, nsamples =100)
#pp_check shows the distributions of many replicated data sets drawn from the posterior predictive distribution (thin light curves) compared with the empirical distribution of the observed outcome (the thick dark curve).


#Posterior predictive intervals: good data fitting
pp_check(model_Flog, type = "intervals_grouped", group = "Function")
#in this case, the model predicts the raw response data reasonably well

#Posterior predictive fit: good data fitting
pp_check(model_Flog, type = "scatter_avg_grouped", group = "Function") + 
  geom_abline(intercept = 0, slope = 1 , color = "red", lty = 2)
####  Scatterplots (and similar visualizations) of the data y vs. individual simulated datasets (rows) in yrep, or vs. the average value of the distributions of each data point (columns) in yrep.
#in this dataset, its clear that the "spear" and "spiny" groups are well-fitted but the "harrow" group is more poorly fitted



#################################################
# "MODEL SELECTION"
#Reduced model without phylogeny

#Reduced model without random phylogeny
model_Flog2 <-brm(Dlog~Function+(1|Species),
                  data = bark2, family = gaussian(), sample_prior = TRUE,
                  iter = 20000,warmup = 10000, chains = 4, cores = 4,
                  thin = 10, save_all_pars = TRUE,seed=T, 
                  control=list(adapt_delta=0.99))


#saveRDS(model_Flog2, "bark/Dlog_Int_Spec.rds")
#model_Flog2 <- readRDS("bark/Dlog_Int_Spec.rds")

summary(model_Flog2)


#plot the posterior fixed effect estimates with credible intervals over the data
quartz()
plot(conditional_effects(model_Flog), points = TRUE)

#Posterior estimates and confidence intervals: tests of significance
mcmc_plot(model_Flog)




####Model selection
#Leave-one-out cross-validation: compare models performance
loo(model_Flog,model_Flog2)
##Models with higher elpd (expected log posterior density) are considered to be better models
#in our case the models are generating very similar results, so simpler model better

#also look up the function 'loo_compare()' (requires adding additional specifications into the model fitting)
#loo_compare(model_Flog,model_Flog2,criterion="waic")

#more reading:
#https://www.rdocumentation.org/packages/brms/versions/2.15.0


#############
# a different evaluation method: 
# RMSPE (root mean square percentage error)

#get fitted values
Yrep.Flog<- as.data.frame(fitted (model_Flog))
summary(Yrep.Flog)
Yrep.Flog2<- as.data.frame(fitted (model_Flog2))

obs <- bark2$Dlog

y.fit.Flog <- Yrep.Flog$Estimate
y.fit.Flog2 <- Yrep.Flog2$Estimate
comb <- data.frame(obs,y.fit.Flog,y.fit.Flog2)
summary(comb)


#RMSPE tests
RMSPE(y_pred=y.fit.Flog,y_true=obs)
RMSPE(y_pred=y.fit.Flog2,y_true=obs)

#a smaller RMSPE means smaller differences between fitted and observed, so smaller is better!
#first model fits the data slightly better, but the advantage is minute (~1%)



#################################################
### A proportional data problem
#################################################

#this is a dataset where squirrels are trying to remove bark from the surface of a tree trunk and how this is impeded by spines on the trunk


debark <- read.csv("bark/debarking.csv", header = TRUE, sep = ",", row.names = NULL)
summary(debark)

debark$Species <- as.factor(debark $Species)
debark $Family <- as.factor(debark $Family)
debark $Strategy <- as.factor(debark $Strategy)

summary(debark)

#please note the form of the response data. it is percentage data for the percentage of bark removed off the trunk of a particular species by a squirrel
#this data is therefore bound below and above, so could be well described by a binomial distribution or a beta distribution
#the data is really based on a continuous variable, so the beta distribution might be better than the binomial

summary(debark)


#now lets get the phylogeny and make sure the data and phylo match up

#different route
#first reduce data to list of common taxa with tree
n1<-unique(phylo$tip.label)
str(n1)
debark2 <- debark[debark$Species %in% n1,]
dim(debark)
dim(debark2) 
length(unique(debark2$Species)) 

#so, second reduce tree to match common taxa with data
n2<-unique(debark2$Species)
phylo3<-keep.tip(phylo,phylo$tip.label[match(n2, phylo$tip.label)])
summary(phylo3)

#make the covariance matrix

#make the covariance matrix from the phylogeny
phylo_cor3 <- ape::vcv.phylo(phylo3)
summary(phylo_cor3)

#make an additional factor to reference the species in the phylogeny
debark2$phylo3 <- debark2$Species

dim(debark2)


#OK now we can run the model
#i want to compare the beta model against the binomial model
debark2$Prop.S.int <- as.integer(debark2$Prop.S*10)

summary(debark2)

#so lets force the percentage data into proportional data
debark2$Prop.S.01=as.numeric(debark2$Prop.S.int)/1000
summary(debark2$Prop.S.01)

#now for the model
m_debark_beta <-brm(Prop.S.01~Strategy+(1|phylo3), 
                 cov_ranef = list(phylo3 = phylo_cor3),
                 data = debark2, family = Beta(), sample_prior = TRUE,
                 iter = 2000,warmup = 1000, chains = 4, cores = 4,
                 thin = 1, save_all_pars = TRUE, seed = T,
                 control=list(adapt_delta=0.95, max_treedepth=15))
#please note we are using very few samples here in order to run the model to completion quickly. in reality you should sample much more

summary(m_debark_beta)



#saveRDS(m_debark_beta, "bark/debark_squirrel_beta.rds")

#m_debark_beta <- readRDS("bark/debark_squirrel_beta.rds")


########################################
#now do model evaluation

summary(debark2)

#R-squared: variance explained by model
bayes_R2(m_debark_beta)


#now look at the posterior distributions generated for each estimate and the chain value (this should form nice horizontal bands called "hairy caterpillars"; if not flat then sampling possibly insufficient or model underspecified)
plot(m_debark_beta, N = 6, ask = FALSE)


#Posterior estimates and confidence intervals: tests of significance
mcmc_plot(m_debark_beta)
#output shows 99% (wider) and 95% (narrower) credible intervals


#proceed with other tests and model comparisons as before..





