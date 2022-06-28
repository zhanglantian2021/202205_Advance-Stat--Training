#############################
## Inference with LMMs
#############################

## CODE adapted from Robert Bagchi (University of Connecticut)

rm(list=ls())
options(digits=3, width=60)


## load the relevant libraries
library(lme4)
library(ggplot2)
library(lattice)
library(arm)
library(utils)
library(car)
library(MuMIn)


?pvalues


##********* BACK TO POWERPOINT *********##

######################################
#### Code 5.1: inferences in lmer and p values
#######################################
getwd()

biodepth <- read.csv("exampledata/Biodepth.csv")
summary(biodepth)
dim(biodepth)

biodepth$location<-as.factor(biodepth$location)
biodepth$block<-as.factor(biodepth$block)

mod1 <- lmer(sqrt(biomass)~log.diversity + (1|location), data=biodepth)
summary(mod1)

anova(mod1) ## no denominator dfs or p values


######################################################
#likelihood ratio method
#LRT?
mod0 <- update(mod1, ~.-log.diversity)
summary(mod0)


#LRT
anova(mod0, mod1)

#######################################################
# now apply denominator df estimation.


## so lets make a stab at the denominator df

getME(mod1, 'devcomp')$dims 
## number of random effect parameters = nth = 1
?getME

#N (# observations/ sampling units)
#p (# fixed parameters)= 2, 
#nth (# random effect variances) = 1,
#q (# random effect levels) = 8
# REML = 2

#first, we will be non-conservative
#if we assume fixed effects parameters only plus only one df for the random effects' variance (nth) then: 
den.df <- 476 - 2 - 1 - 2   # n - p - 1 random variance - 2 REML
den.df  #471

## now get F-value for fixed effect
anova(mod1)   #F=47.5

## perform an f-test for log.diversity effect
pf(47.5, df1=1, df2=den.df, lower.tail=FALSE) ## highly significant

## more conservatively, let's assume each random level takes a df. Then: 

den.df <- 476  -2 - 8 - 2  # n - p - q random levels - 2 REML
den.df  #464
pf(47.5, 1, den.df, lower.tail=FALSE) 
## in this case there is not much difference here from the previous calculation. Why? 
## (think about the total df..)


#there are packages that run each of the major estimation methods

#first, sattherthwaite estimation

#just load the lmerTest package

library(lmerTest)

#rerun the basic codes:
mod1 <- lmer(sqrt(biomass)~log.diversity + (1|location), data=biodepth)
summary(mod1)

anova(mod1)


#second, Kenward-Rogers estimation

## can estimate them from the number of replicates
summary(mod1) ## 476 observations


## Using the Kenward-Rodgers approximation of den DFs
library(pbkrtest)

KRmodcomp(mod0, mod1)
#check the ddf number
464+(471-464)/2  
#468. This random effects model is quite simple. Hence our ability to calculate the answer. The calculations are more convoluted for complex random effects (nesting, crossed)

## the 'car' library also allows one to test mixed models with KR df
Anova(mod1, test='F')





##********* BACK TO POWERPOINT *********##


#############################
## Code 5.2
##############################

## finding confidence intervals for parameter estimates
## using bootstrapping is one way to evaluate the significance
## of parameter estimates (both means and variances)

##To understand bootstrapping you might find this useful

##http://influentialpoints.com/Training/nonparametric-or-parametric_bootstrap.htm


summary(mod1)
## can use the confint function to conduct bootstraps
confint.result <- confint(mod1, method='boot', oldNames=F, nsim = 2999) ## has a number of options
confint.result
## if parameter's CIs cross 0, then they are not significantly different from 0
## don't worry about the warnings; these are evaluation methods

## We can also do this 'by hand'  ;-)

## produce a number of models
mods <- replicate(499, {
  newresp <- simulate(mod1)
  newmod <- refit(mod1, newresp)}, 
  simplify=FALSE)

fixef.sims <- sapply(mods, fixef) 
fixef.sims
## will give you the fixed effect for all models
class(fixef.sims)

## to get confidence intervals, used ordered quantiles
apply(fixef.sims, 1, quantile, c(0.025, 0.975))


#compare this to the confint() result
confint.result

#Plant damage codes







