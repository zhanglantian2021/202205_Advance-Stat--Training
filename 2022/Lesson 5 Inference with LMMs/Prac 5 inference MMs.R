#############################
## Inference with LMMs
#############################

## CODE adapted from Robert Bagchi (University of Connecticut)
setwd("E:/academic_resources/Advance Stat Training/2022/Lesson 5 Inference with LMMs")


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


biodepth <- read.csv("/Users/kyletomlinson/Dropbox/Teaching/XTBG advanced stats 2022/Lectures/Lesson 5 Inference with LMMs/exampledata/Biodepth.csv")
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

#plantdamage dataset
damage <- read.csv("exercises/plantdamage.csv")
summary(damage)
dim(damage)

damage$shadehouse <- as.factor(damage$shadehouse)
damage$light <- as.factor(damage$light)

quartz(); hist((damage$growth))

damage$growth.cubrt <- (damage$growth)^(1/3)

head(damage,n=15)

unique(damage$shadehouse)

lmm1 <- lmer(growth~light*damage+(1|shadehouse),damage)
summary(lmm1)

anova(lmm1)

Anova(lmm1, test='F')

getME(lmm1, 'devcomp')$dims 
n=103
p=4 #(intercept,light,damage,lightXdamage)
q=10 #shadehouses
nth=1 #random effect variance
REML=4

#conservative
#shadehouse level
10-2=8
#inside shadehouse
103-4-10-4 = 85

#这两个的区别：
Anova(lmm1)
Anova(lmm1, test='F')

#light (at level of greenhouse)
pf(38.2, 1, 8, lower.tail=FALSE) 

#damage (inside greenhouse)
pf(47.6, 1, 85, lower.tail=FALSE) 
#damage:light (inside greenhouse)
pf(47.6, 1, 85, lower.tail=FALSE) 

#everything is significant
#now look at summary to understand the direction of effects

summary(lmm1)

#############################
## Exercide 5.2
##############################

summary(damage)

lmm2 <-lmer(growth~light*damage+(1+damage|shadehouse),damage)
summary(lmm2)

#first evaluate random effects

confint.result <- confint(lmm2, method='boot', oldNames=F, nsim = 2999) ## has a number of options
confint.result
#model does not need to be simplified because:
#1. teh random slopes variance is different from zero
#2. the interaction in the fixed effects is significant
#so either report the bootstrapped results or rerun lmer() and get the expected coefficients and their t-probabilities.


summary(lmm2)
Anova(lmm2,test="F")

#also report variance explained
r.squaredGLMM(lmm2)




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








