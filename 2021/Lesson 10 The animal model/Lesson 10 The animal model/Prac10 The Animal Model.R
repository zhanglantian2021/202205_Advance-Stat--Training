
library(brms)
library(ggplot2)
library(ape)
library(caper)
library(sjPlot)

#devtools::install_version('rstan', version = '2.19.3', repos = "http://cran.us.r-project.org")

#################################################
#### load trait data #####
#################################################
setwd("/Users/kyletomlinson/Dropbox/Teaching/XTBG advanced stats 2021/Lectures/Lesson 10 The animal model")

# database with the corrected data
bark <- read.csv("bark/bark_nutrition.csv", header = TRUE, sep = ",", row.names = NULL)
summary(bark)

bark$Species <- as.factor(bark$Species)
bark$Family <- as.factor(bark$Family)
bark$Month <- as.factor(bark$Month)
bark$Function <- as.factor(bark$Function)

summary(bark) #note that there are multiple individuals per species
dim(bark)

length(unique(bark$Species)) #53 species

#################################################
#### load phylogeny data #####
#################################################


phylo <- read.tree("/Users/kyletomlinson/Dropbox/Teaching/XTBG advanced stats 2021/Lectures/Lesson 10 The animal model/bark/scenario.2_run.1.tre")
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
#doesnt work because the function cannot handle repeats

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
#model_Flog <-brm(Dlog~Function+(1|phylo2)+(1|Species), cov_ranef = list(phylo2 = phylo_cor),data = bark2, family = gaussian(link="identity"),sample_prior = TRUE,iter = 10000,warmup = 5000, chains = 4, cores = 4,thin = 10, save_all_pars = TRUE,seed=T,control=list(adapt_delta=0.99))

model_Flog <-brm(Dlog~Function+(1|phylo2)+(1|Species), 
                 cov_ranef = list(phylo2 = phylo_cor),
                 data = bark2, family = gaussian(link="identity"), sample_prior = TRUE,
                 iter = 10000,warmup = 5000, chains = 4, cores = 4,
                 thin = 10, save_all_pars = TRUE,seed=T,
                 control=list(adapt_delta=0.99))

#https://stackoverflow.com/questions/63025052/brm-model-compiling-but-returning-model-object
model_Flog

#you can save the output (a good idea when your model has a lot of data )
#saveRDS(model_Flog, "bark/Dlog_Int_Phyl_Spec.rds")

#the next comman allows you to read in the saved output
#model_Flog <- readRDS("bark/Dlog_Int_Phyl_Spec.rds")


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

#R-squared: variance explained by model
bayes_R2(model_Flog)


#output table including conditional and marginal R2 
library(sjPlot)
tab_model(model_Flog, show.intercept = T, show.r2 = TRUE, 
          transform = NULL,show.re.var=T)


#################################################

#what follows are some tests evaluating the performance of the analysis

#Posterior predictive distribution: good data fitting
pp_check(model_Flog, nsamples =100)
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
#sum of random effect variances/sum of random plus residual variances

#now just estimate variance explained by phylogeny
#icc(model_Flog,ppd=T,re.form = ~(1| phylo2),robust=T)



#################################################

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
plot(conditional_effects(model_Flog), points = TRUE)

#Posterior estimates and confidence intervals: tests of significance
mcmc_plot(model_Flog, pars = c("^b_", "^sd_"))




####Model selection
#Leave-one-out cross-validation: compare models performance
loo(model_Flog,model_Flog2,reloo=F)
##reloo=T takes a lot of time!!
##Models with lower LOO are considered to be better models, particularly if the difference in LOO between two models (ΔLOO) is more than twice its standard error.



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

summary(debark2)


#OK now we can run the model
#i want to compare the beta model against the binomial model
debark2$Prop.S.int <- as.integer(debark2$Prop.S*10)

summary(debark2)

#so lets force the percentage data into proportional data
debark2$Prop.S.01=as.numeric(debark2$Prop.S.int)/1000


m_debark_beta <-brm(Prop.S.01~Strategy+(1|phylo3), 
                 cov_ranef = list(phylo3 = phylo_cor3),
                 data = debark2, family = Beta(), sample_prior = TRUE,
                 iter = 2000,warmup = 1000, chains = 4, cores = 4,
                 thin = 10, save_all_pars = TRUE, seed = T,
                 control=list(adapt_delta=0.95, max_treedepth=15))

summary(m_debark_beta)

saveRDS(m_debark_beta, "bark/debark_squirrel_beta.rds")

m_debark_beta <- readRDS("bark/debark_squirrel_beta.rds")


debark2$Prop.S.int



m_debark_binom <-brm(Prop.S.int|trials(1000)~Strategy+(1|phylo3), 
                 cov_ranef = list(phylo3 = phylo_cor3),
                 data = debark2, family = binomial(), sample_prior = TRUE,
                 iter = 2000,warmup = 1000, chains = 4, cores = 4,
                 thin = 10, save_all_pars = TRUE, seed = T,
                 control=list(adapt_delta=0.95, max_treedepth=15))

summary(m_debark_binom)

saveRDS(m_debark_binom, "bark/debark_squirrel_binom.rds")

m_debark_binom <- readRDS("bark/debark_squirrel_binom.rds")


########################################
#now lets compare the two models


#now do model evaluation

library(MLmetrics)
library(MuMIn)
library(ggplot2)

summary(debark2)

#R-squared: variance explained by model
bayes_R2(m_debark_beta)
#R-squared: variance explained by model
bayes_R2(m_debark_binom)



#Posterior predictive intervals: good data fitting
pp_check(m_debark_beta, type = "intervals_grouped", group = "Strategy")
quartz()
pp_check(m_debark_binom, type = "intervals_grouped", group = "Strategy")


#Posterior predictive fit: good data fitting
pp_check(m_debark_beta, type = "scatter_avg_grouped", group = "Strategy") + 
  geom_abline(intercept = 0, slope = 1 , color = "red", lty = 2)

pp_check(m_debark_binom, type = "scatter_avg_grouped", group = "Strategy") + 
  geom_abline(intercept = 0, slope = 1 , color = "red", lty = 2)



#get fitted values
Yrep.beta<- as.data.frame(fitted (m_debark_beta))
Yrep.binom<- as.data.frame(fitted (m_debark_binom))


obs <- debark2$Prop.S.01
y.beta <- Yrep.beta$Estimate
y.binom <- Yrep.binom$Estimate/1000

comb <- data.frame(obs,y.beta,y.binom)

#plot the models

ggplot(data=comb, aes(x=obs, y=y.beta)) + geom_point(col="red")+
  geom_point(data=comb,aes(x=obs, y=y.binom),
              stat='identity')


par(mfrow=c(1,2))  
plot(y.beta~obs,data=comb,col="red")
abline(a=0,b=1,col="black")
plot(y.binom~obs,data=comb,col="blue")
abline(a=0,b=1,col="black")
#beta model does not make sense


#RMSPE tests
RMSPE(y_pred=y.beta,y_true=obs)
RMSPE(y_pred=y.binom,y_true=obs)




#########################################################
#what about dropping the random effect?
m_debark_binom2 <-brm(Prop.S.int|trials(1000)~Strategy, 
                 data = debark2, family = binomial(), sample_prior = TRUE,
                 iter = 2000,warmup = 1000, chains = 4, cores = 4,
                 thin = 10, save_all_pars = TRUE, seed = T,
                 control=list(adapt_delta=0.95, max_treedepth=15))

summary(m_debark_binom2)

#Posterior predictive intervals: good data fitting
quartz()
pp_check(m_debark_binom2, type = "intervals_grouped", group = "Strategy")
#generates junk!!






