##########################################
## Part of code from Robert Bagchi
## 13-19/7/2014
##########################################


## Prac 9: Generalised least squares

setwd("E:/academic_resources/Advance Stat Training/2022/Lesson 8 GLS")


rm(list=ls())
library(nlme)
library(ggplot2)


####################################################
Code 8.1 - weighting residuals to improve behaviour
####################################################

#here we will use some simulation to show some simple error structures that can be dealt with using generalised least squares
#we will simulate a plant growth experiment with irrigation and soil type treatments

## Let's make up some data :-)
irr2 <- data.frame(water=runif(100, 0, 10), soil=sample(c('loam', 'sand'), 100, replace=T))

dim(irr2)         
summary(irr2)  

irr2$soil <- as.factor(irr2$soil)

## make up a heterscedastic error where the variance increases
## with the fitted values
beta <- c(5, 0.4)   #these are the coefficients of our regression relationship ;-)
X <- model.matrix(~water, data=irr2)
X   #now it is set up as the components of a regression model ~ 1 + B*Water

#predict y values using X and the betas: y = X*beta
irr2$y <- as.vector(X %*% beta)
summary(irr2)

plot(irr2$water,irr2$y)

#even variance
#irr2$y_even <- irr2$y + rnorm(100, mean=0, sd=sqrt(0.1))
#plot(irr2$water,irr2$y_even)


## make the variance increase with the mean
irr2$y <- irr2$y + rnorm(100, mean=0, sd=sqrt(0.01*irr2$y^3))
plot(irr2$water,irr2$y)

#y2 <- as.vector(X %*% beta)
#irr2
#plot(irr2$y~y2)

#now lets analyse our fake data using gls

#first use a model assuming no covariance error structure
modlm <- lm(y ~ water, data=irr2)
summary(modlm)

mod1 <- gls(y ~ water, data=irr2)
summary(mod1)

coef(summary(modlm))
coef(summary(mod1))

plot(mod1, (abs(resid(.)))^(1/2)~fitted(.,type='pearson'), type=c('smooth', 'p'),abline=0) 

## obvious trend

## now model the variance as a power function of the fitted values   
mod1b <- gls(y ~ water, data=irr2, weights=varPower(form=~fitted(.)))
# mod1c <- gls(y ~ water, data=irr2, weights=varFixed())
# AIC(mod1b,mod1c)
            
#look at output, note additional section for variance function and parameter estimate            
summary(mod1b)

plot(mod1b, abs(resid(., type='pearson'))~fitted(.), type=c('smooth', 'p'), abline=0) ## better


## ......................................

## can also have a relationship between the variance and a covariate

## make the variance increase with the irrigation amount
irr2$y <- as.vector(X %*% beta)
irr2$y <- irr2$y + rnorm(100, mean=0, sd=irr2$water^0.5)

plot(y~water, data=irr2) 

## initially you might diagnose this as a exponential growth problem
## but a histogram of y does not yield an exponential pattern
hist(irr2$y)



mod2 <- gls(y~water, data=irr2)
plot(mod2, abs(resid(., type='p'))~fitted(.), type=c('smooth', 'p'))
## not good, so lets test whether there is sime relationship with the predictor
mod2a <- update(mod2, weights=varPower(form=~fitted(.)))

## improves some diagnostics
plot(mod2a, abs(resid(., type='p'))~fitted(.), type=c('p', 'r'))

summary(mod2a) 
## check the power

## ......................................



beta1 <- c(5, 0.4,0.1)   #these are the coefficients of our regression relationship ;-)
X1 <- model.matrix(~water+soil, data=irr2)
X1   #now it is set up as the components of a regression model ~ 1 + B*Water

irr2$y <- as.vector(X1 %*% beta1)
irr2$y <- irr2$y + rnorm(100, mean=0, sd=irr2$water^0.5)

mod2b <- gls(y~water+soil, data=irr2)
plot(mod2b, abs(resid(., type='p'))~fitted(.), type=c('smooth', 'p'))
mod2b2 <- gls(y~water+soil, data=irr2, weights=varPower(form=~fitted(.)))

## improves some diagnostics
plot(mod2b2, abs(resid(., type='p'))~fitted(.), type=c('smooth', 'p'))




### so, lets make it that the sd is different between soil types - 2x in sand
irr2$y <- rnorm(100, mean=X %*% beta, sd=ifelse(irr2$soil=='loam', 1, 2))
plot(y~soil, data=irr2)

mod3 <- gls(y~water, data=irr2)
plot(mod3, (abs(resid(., type='p')))^(1/2)~fitted(.),type=c('smooth', 'p'))
plot(mod3, (abs(resid(., type='p')))^(1/2)~fitted(.)|soil, type=c('smooth', 'p')) ## variance larger in sand
## now allow different variances for soil types
mod3a <- update(mod3, weights=varIdent(form=~1|soil))
plot(mod3a, (abs(resid(., type='p')))^(1/2)~fitted(.)|soil,type=c('smooth', 'p')) ## fixed!


#plot(mod3a, (abs(resid(., type='p')))^(1/2)~fitted(.)|soil,type=c('p', 'r')) ## fixed!

quartz()

summary(mod3a) ## variance greater in sand 
## please note!~data is simulated so actual values will change 



########################################
## Example with random effects

#make a grouping predictor
irr2$grp <- factor(sample(1:10, nrow(irr2),replace=T) )
unique(irr2$grp)
irr2
summary(irr2)

#remember: y = X*beta+ Z*b + e

#introduce variation into data produced by random effect levels
sigma.b <- 0.2
b1 <- rnorm(10, sd=sigma.b)
b1
Z <- model.matrix(~ 0+grp, data=irr2) 

irr2$y <- as.vector(X%*%beta + Z %*% b1)

par(mfrow=c(1,1))
plot(y~water,irr2)
dim(irr2)
irr2$y <- irr2$y + 0.05*rnorm(100, mean=0, sd=irr2$y^2)
summary(irr2)
plot(y~water,irr2)

#now run the model using lme
#first without assuming heteroscedastic errors
mod4 <- lme(y~water, random=~1|grp, data=irr2)

plot(mod4, (abs(resid(., type='p')))^(1/2)~fitted(.),type=c('smooth', 'p')) ## increasing variance
     
#now assuming heteroscedastic errors     
mod4a <- update(mod4, weights=varPower(form=~fitted(.)))
plot(mod4a, abs(resid(., type='p'))~fitted(.), type=c('p', 'r')) ## better
coef(summary(mod4a))
coef(summary(mod4))


AIC(mod4a,mod4)
## this exercise shows you how to construct and test the models
## and how to understand systematic ways errors can accumulate in your models that are not easily diagnosed

#evaluate whether the model is better with the covariance structure?
anova(mod4a,mod4)



####################################################
			     RETURN TO PPT
####################################################


##################################
## Code 8.2.  Spatial correlation structures
######################################

#the data is for total basal area in about 48 forest plots distributed in forest fragments around Menglun
#there are different recognised forest types in the area
#the question is: do forest types differ in basal area?


treebase <- read.csv('exampledata/Foresttype_BA.csv')
summary(treebase)
head(treebase,n=10)

#convert characters to factors
treebase $Plot <- as.factor(treebase $Plot)
treebase $Fragment <- as.factor(treebase $Fragment)
treebase $forest.type <- as.factor(treebase $forest.type)
summary(treebase)

ggplot(data= treebase, aes(x=longitude, y=latitude, size=BA.ha, colour= forest.type)) + geom_point()

hist(treebase $BA.ha)
hist(log(treebase $BA.ha))

m.rs1 <- gls(log(BA.ha)~ forest.type, data= treebase)

summary(m.rs1) #some group differences significant
plot(m.rs1) ## residual homoscedasticity looks okay


## We can assess spatial autocorrelation using a variogram
## which shows how data variance is related to distance 
## between observations
## if there is no spatial autocorrelation, 
## then we expect to see a flat line
plot(Variogram(m.rs1, form=~longitude+latitude, resType='n')) 
#our line doesnt look too bad, just a bit cranky


anova(m.rs1)

## increase in variance with distance...
## so now include information on the spatial distrbution between points
m.rs2 <- update(m.rs1, correlation=corExp(form=~ longitude + latitude))
anova(m.rs1, m.rs2) ## improves model fit
plot(Variogram(m.rs2, form=~longitude+latitude, resType='n')) 
#line has flattened a bit
plot(m.rs2) ## residual homoscedasticity looks okay
summary(m.rs2)
anova(m.rs2) #groups no longer different
#不需要relevel的方法：
contrast(emmeans(m.rs2, specs="forest.type"),method="pairwise")


AIC(m.rs1,m.rs2)

##include a nugget variance
## a nugget variance is an additional variance term 
## included to account for fine-scale variation at zero lag distance
## or observation error. 

m.rs3 <- update(m.rs2, correlation=corExp(form=~longitude + latitude, nugget=T))
anova(m.rs1, m.rs2, m.rs3) #nugget makes no change in this case
plot(Variogram(m.rs3, form=~longitude + latitude, resType='n')) 
## variogram barely moves

AIC(m.rs1,m.rs2,m.rs3)
#so the model with spatial correction and no nugget seems to fit the data best of the 3 models tested
#and we conclude that forest type has no effect on plot basal area

library(emmeans)


##OTHER USEFUL DATASETS
http://plantecology.syr.edu/fridley/bio793/mixed2.html

#interpretation of variograms
http://geog.uoregon.edu/GeogR/topics/variograms.html









