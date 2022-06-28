#GLMM EXERCISES

rm(list=ls())
library(lme4)
library(ggplot2)

setwd("E:/academic_resources/Advance Stat Training/2022/Lesson 7 GLMMs")

polls <- read.csv('exercises/pollsdata_day2.csv')
summary(polls)

polls$state <- as.factor(polls$state)
polls$edu <- as.factor(polls$edu)
polls$black <- as.factor(polls$black)
polls$female <- as.factor(polls$female)

glm1 <- glm(bush~black+female,polls,family=binomial())
summary(glm1)


############################

#plot the data by states

ggplot(data=polls, aes(x=black, y=bush,group=female)) +
  geom_point() +
  geom_smooth(method='lm') +
  facet_wrap(~state)  

glmm1 <- glmer(bush~black+female+(1|state),polls,family=binomial())
summary(glmm1)


##################################

glmm2 <- glmer(bush~black+female+(1+female|state),polls,family=binomial())
summary(glmm2)

confint(glmm2,method='boot', oldNames=F, nsim = 199)
#dont do lots unless it's your final analysis!! because glm bootstraps take ages!!!!!   and ages !!!!!     and ages !!!!!!
#if you do enough bootstraps (i did 499), you find that the slope random effect IS significant
#if you do too few bootstraps then youn might find the slopes are non-significant... so always use large boostrap numbers for final analysis

#########################

#diagnostics

## plot the sqrt of the absolute residuals against fitted values
plot(sqrt(abs(resid(glmm2)))~ fitted(glmm2))
lines(lowess(sqrt(abs(resid(glmm2)))~
               fitted(glmm2)), col='red')

#overdispersion??
#binned plot
library(arm)
x <- predict(glmm2)
y <- resid(glmm2)
par(mfrow=c(1,1))
quartz();binnedplot(x,y)

#overdispersion test
sum(resid(glmm2, type='pearson')^2)/df.residual(glmm2)

#check normality of random effects
par(mfrow=c(1,2))
qqPlot(ranef(glmm2)$state$'(Intercept)')
qqPlot(ranef(glmm2)$state$'female1')


##########################

pests <- read.csv('exercises/Pesticide.csv')
summary(pests)
head(pests,n=10)
dim(pests)

pests$station <- as.factor(pests$station)
pests$trt <- as.factor(pests$trt)

#response variable is N (number of seedling recorded)


ggplot(data=pests, aes(x=trt, y=N,col=trt)) +
  geom_point() +
  geom_smooth(method='lm') +
  facet_wrap(~station)  

#probably we need random slopes.. but for the sake of learning we will ignore them for now and just have a random intercepts model

glmp1 <- glmer(N~trt+(1|station),pests,family=poisson())
summary(glmp1)

#formally test with confint..

#check diagnostics

#overdispersion test
sum(resid(glmp1, type='pearson')^2)/df.residual(glmp1)

#try a covariate: seed rain
hist(pests$Nseeds)
pests$Nseeds.log <- log(pests$Nseeds)


glmp2 <- glmer(N~trt+Nseeds.log+(1|station),pests,family=poisson())
summary(glmp2)

#overdispersion test
sum(resid(glmp2, type='pearson')^2)/df.residual(glmp2)


#try an observation  level random effect (OLRE)
pests$indx <- 1:nrow(pests)
#designed to draw off additional variance that is assumed to be normal

glmp3 <- glmer(N~trt+(1|station) + (1|indx),pests,family=poisson())
summary(glmp3)

#overdispersion test
sum(resid(glmp3, type='pearson')^2)/df.residual(glmp3)

#check normality of random effects
par(mfrow=c(1,2))
qqPlot(ranef(glmp3)$station$'(Intercept)')
qqPlot(ranef(glmp3)$indx$'(Intercept)')


#lets model this with random slopes instead
glmp4 <- glmer(N~trt+(1+trt|station),pests,family=poisson())


#overdispersion test
sum(resid(glmp4, type='pearson')^2)/df.residual(glmp4)

#also good


summary(glmp4)

#random slopes model is probably better than an OLRE model because it can be easily explained

#now bootstrap
confint(glmp4,method='boot', oldNames=F, nsim = 199)

#########################
#########################

#open the seedling survival data

getwd()
survive <- read.csv('exercises/plantdamage2.csv')
summary(survive)
dim(survive)
survive$shadehouse <- as.factor(survive$shadehouse)
survive$light <- as.factor(survive$light)
length(unique(survive$shadehouse))

head(survive)

glmm1 <- glmer(cbind(survs,4-survs)~light*damage+(1|shadehouse),survive,family=binomial())
summary(glmm1)

glmm2 <- glmer(cbind(survs,4-survs)~light+damage+(1|shadehouse),survive,family=binomial())
summary(glmm2)
anova(glmm1,glmm2) #LRT

glmm3 <- glmer(cbind(survs,4-survs)~light+(1|shadehouse),survive,family=binomial())
summary(glmm3)
anova(glmm2,glmm3) #LRT

confint(glmm3,method='boot', oldNames=F, nsim = 199)

#so light survives!!

#any lets just use the full model, glmm1


#assuming light was significant

preddat <- expand.grid(light=c('D', 'L'),
                       damage=c(0, 0.1, 0.25),
                       shadehouse =1:10)
summary(preddat)
dim(preddat)
preddat$shadehouse <-as.factor(preddat$shadehouse)

## prediction with NO random effects
preddat$pred.fix.T <- 
  predict(glmm1, newdata=preddat, re.form=~0) 
## prediction WITH random effects for shadehouse 
preddat$pred.sh.T <- predict(glmm1, newdata=preddat, re.form=NULL) 

head(preddat,n=10)

#now backtransform (plogis)

preddat$pred.fix <-plogis(preddat$pred.fix.T)
preddat$pred.sh <-plogis(preddat$pred.sh.T)

#fixed effect
ggplot(data=survive, aes(x=damage, y=survs/4)) + geom_point()+
  geom_smooth(data=preddat,aes(x=damage, y=pred.fix, col=light),stat='identity')

#shadehouses
ggplot(data=survive, aes(x=damage, y=survs/4)) + geom_point()+
  geom_smooth(data=preddat,aes(x=damage, y=pred.sh, col=light),stat='identity')+facet_wrap(~shadehouse)