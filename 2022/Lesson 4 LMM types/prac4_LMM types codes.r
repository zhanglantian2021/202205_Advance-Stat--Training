## Fitting different random effects models with lmer

## CODE from Robert Bagchi (University of Connecticut)

## preliminaries
rm(list=ls())  

## Set the working directory
setwd("E:/academic_resources/Advance Stat Training/2022/Lesson 4 LMM types")

getwd()
## Load the libraries
library(ggplot2)
library(lme4) ## library for linear mixed effects models in R
library(car)
library(arm)
library(MuMIn)


## load the data
radon <- read.csv('exampledata/Radon_Data_RB.csv', h=TRUE)

## Look at your data
dim(radon)
summary(radon)
head(radon,n=20)
str(radon)

radon$floor <- as.factor(radon$floor)
radon$county <- as.factor(radon$county)
radon$cgroup <- as.factor(radon$cgroup)


par(mfrow=c(1,1))
boxplot(radon~floor, data=radon, notch=T)


## now add a regression line
#group=1?
ggplot(data=radon, aes(x=floor, y=radon,group=1)) + geom_point() +
  geom_smooth(method='lm')

## separate into different panels for each county
ggplot(data=radon, aes(x=floor, y=radon, group=county)) +
  geom_point() +
  geom_smooth(method='lm') +
  facet_wrap(~county)  
##note how the patterns differ by county
##note how the confidence intervals differ by county, partly due to differences in sample sizes
  

##################################
## now fit a simple linear model
##################################
mod.radon.lm1 <- lm(radon~floor, data=radon)

par(mfrow=c(2,2))
plot(mod.radon.lm1)  

summary(mod.radon.lm1) ## quick look at model


##********* BACK TO POWERPOINT *********##


################################################

## CODE 4.1
## Fitting a mixed-effects model with lmer
###############################################
mod.radon.lmer1 <- lmer(radon~floor + (1|county), data=radon)

par(mfrow=c(1,2))
hist(resid(mod.radon.lmer1))  

str(mod.radon.lmer1)
#residuals should be normally distributed
hist(ranef(mod.radon.lmer1)[[1]][,1])
#random variables (intercepts) should be normally distributed
ranef(mod.radon.lmer1)

summary(mod.radon.lmer1) ## summary of the model

coef(summary(mod.radon.lm1))
coef(summary(mod.radon.lmer1))

## explanation in powerpoint for 1 slide.

fixef(mod.radon.lmer1) ## fixed effect coefficients
summary(mod.radon.lmer1)$coef

VarCorr(mod.radon.lmer1) ## The variance components


## the differences between the overall intercept b0 and the intercept for each county are provided by the random effects or BLUPs (best linear unbiased estimators)
ranef(mod.radon.lmer1) ## rather long
str(ranef(mod.radon.lmer1)) ## a list


res1 <- ranef(mod.radon.lmer1)$county[,1]
summary(res1) ## mean ~ 0
hist(res1) ## histogram
abline(v=mean(res1), col='red')


## Random effects should be normally distributed
par(mfrow=c(1,1)); qqPlot(res1)


## to get the expected radiation level for the ground floor in a given county we can simply add the intercept fixed effect to the random effect values
gfloor <- fixef(mod.radon.lmer1)["(Intercept)"] +
            ranef(mod.radon.lmer1)$county$"(Intercept)"
gfloor

# this means real fixed effect
mean(gfloor) ## note the mean is not 0 anymore

hist(gfloor) 
abline(v=mean(gfloor), col=2)  #here you can see the mean is simply the estimate

## How would we calculate the expected radiation level on the first floor?

firstfloor <- fixef(mod.radon.lmer1)["(Intercept)"] +fixef(mod.radon.lmer1)["floor"]+ ranef(mod.radon.lmer1)$county$"(Intercept)"

mean(firstfloor) #equals beta0+beta1

##********* BACK TO POWERPOINT *********##

#######################################
##  CODE 4.2
## Fit a random slope model
######################################

## recall:
## separate into different panels for each county
ggplot(data=radon, aes(x=floor, y=radon, group=county)) +
  geom_point() +
  geom_smooth(method='lm') +
  facet_wrap(~county)  
  ##note how the slopes differ by county


#so add random slopes into the model
mod.radon.lmer2 <- lmer(radon ~ floor + (1+floor|county), data=radon)

summary(mod.radon.lmer2) ## summary of the model



## compare the models
rbind(fixef(mod.radon.lmer1),fixef(mod.radon.lmer2)) ## fixed effects change slightly.


VarCorr(mod.radon.lmer2) # Note the extra row for floor
# Also extra column for correlation between the random effects

summary(ranef(mod.radon.lmer2)$county) ## now two columns here
dim(ranef(mod.radon.lmer2)$county)
##First column is difference between overall intercept and intercept in that county
## Second column is the difference between overall slope and the slope in that county.

head(ranef(mod.radon.lmer2)$county)

## how would we calculate the difference in radon concentration between the ground and first floor in county 2?

fixef(mod.radon.lmer2)["floor1"]+ranef(mod.radon.lmer2)$county$"floor1"

## what is the radon concentration in county 1 on the ground floor?

## What is the radon concentration in county 1 on the first floor?

#either
floor1 <- fixef(mod.radon.lmer2)["(Intercept)"] +ranef(mod.radon.lmer2)$county$"(Intercept)" +fixef(mod.radon.lmer2)["floor1"]+ranef(mod.radon.lmer2)$county$"floor1"

#or
            fixef(mod.radon.lmer2)["(Intercept)"] +
			fixef(mod.radon.lmer2)["floor1"]+
            ranef(mod.radon.lmer2)[[1]][1,1]+
            ranef(mod.radon.lmer2)[[1]][1,2]
#or
            fixef(mod.radon.lmer2)["(Intercept)"] +
			fixef(mod.radon.lmer2)["floor1"]+
            ranef(mod.radon.lmer2)$county[1,1]+
            ranef(mod.radon.lmer2)$county[1,2]



##********* BACK TO POWERPOINT *********##


############################################################
## Code 4.3: assessing mixed effects model fits
##########################################################

## Do they show any trend?
## plot residuals against fitted values
resids <- resid(mod.radon.lmer2, type='pearson')
length(resids)

#homogeneity
plot(resids~fitted(mod.radon.lmer2))
lines(lowess(resids~fitted(mod.radon.lmer2)), col='red')

## are they comparable across the range of predictors?
boxplot(resids~radon$floor)

## Are the residuals homoscedastic?
## plot the sqrt of the absolute residuals against fitted values
plot(sqrt(abs(resids))~ fitted(mod.radon.lmer2))
lines(lowess(sqrt(abs(resids))~
               fitted(mod.radon.lmer2)), col='red')

## and normally distributed?
qqPlot(resids)

## But what about the random effects - are they normally distributed?
## and if we have random slopes they should be normally distributed too!
par(mfrow=c(1,2))
qqPlot(ranef(mod.radon.lmer2)$county$'(Intercept)')
qqPlot(ranef(mod.radon.lmer2)$county$floor)

## to get an R^2 type statistic (not equivalent to OLS R-square)
r.squaredGLMM(mod.radon.lmer2)

#What is it doing?
?r.squaredGLMM
## returns two values: 
#marginal: explanatory power of just fixed effects component
#conditional: measure of explanatory power of random AND fixed effects
## paper by Johnson 2014 



##********* BACK TO POWERPOINT *********##

############################################################
## 4.4 Nested models
############################################################

## are my data nested?
summary(radon)
table(radon$cgroup, radon$county) ## each county is only occurring once

?table

## just to be sure
colSums(table(radon$cgroup, radon$county)!=0) ## all values are == 1
## or (assuming lme4 is loaded)
isNested(radon$county, radon$cgroup)
## county is nested in cgroup
isNested(radon$cgroup, radon$county) # cgroup is not nested in county!

### fitting a  nested model
mod.radon.lmer3 <- lmer(radon~floor + (1|cgroup/county), data=radon)
summary(mod.radon.lmer3)

## alternatively the model can be formulated like this
mod.radon.lmer3b <- lmer(radon~floor + (1|cgroup) +  (1|cgroup:county),data=radon)
rbind(fixef(mod.radon.lmer3),fixef(mod.radon.lmer3b)) ## identical
cbind(VarCorr(mod.radon.lmer3), 
      VarCorr(mod.radon.lmer3b)) ## also identical




##********* BACK TO POWERPOINT *********##











##********* BACK TO POWERPOINT *********##

plantsoil <- read.csv('exercises/plantsoil.csv', h=TRUE)
summary(plantsoil)

#convert to factors
plantsoil$Plot <- as.factor(plantsoil$Plot)
plantsoil$soil <- as.factor(plantsoil$soil)
plantsoil$species <- as.factor(plantsoil$species)
summary(plantsoil)
length(unique(plantsoil$species))  #14 species

#now check distribution of response variable
par(mfrow=c(1,1)); hist(plantsoil$growth)
#looks pretty normal

#now make mixed model
mod1 <- lmer(growth~soil+(1|Plot)+(1|species),plantsoil)
summary(mod1)
#fixed的结果显示就是第一节课的多个x，鱼的例子，就是基准线不一样，系数方程不变化，但是两两线之间比较发生变化。

plantsoil$soil <- relevel(plantsoil$soil,ref="M")


plot(growth~soil,plantsoil)

#check diagnostics
resids <- resid(mod1, type='pearson')
length(resids)


## Are the residuals homoscedastic?
## plot the sqrt of the absolute residuals against fitted values
plot(sqrt(abs(resids))~ fitted(mod1))
lines(lowess(sqrt(abs(resids))~
               fitted(mod1)), col='red')

## and normally distributed?
qqPlot(resids)

## But what about the random effects - are they normally distributed?
## and if we have random slopes they should be normally distributed too!
par(mfrow=c(1,2))
qqPlot(ranef(mod1)$species$'(Intercept)')
qqPlot(ranef(mod1)$Plot$'(Intercept)')

## to get an R^2 type statistic (not equivalent to OLS R-square),因为是线性模型，所以可以用R方来表明每个解释部分占多少。
r.squaredGLMM(mod1)

#last question is a random slopes model
mod2 <- lmer(growth~soil+(1|Plot)+(1+soil|species),plantsoil)
summary(mod2)
#这个模型上面的结果显示species这个随机效应影响y很多（看R2m和R2c），所以看看它是否会影响斜率，而这个模型结果显示species但不影响随机效应的斜率（也就是该随机效应不影响y和固定效应的关系），因为只用0.003和0.019


