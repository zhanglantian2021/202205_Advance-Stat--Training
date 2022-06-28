setwd("E:/academic_resources/Advance Stat Training/2022/Lesson 7 GLMMs")

library(lme4)
library(ggplot2)
?qqplot

polls <- read.csv("exercises/pollsdata_day2.csv")

polls$state <- as.factor(polls$state )
polls$edu <- as.factor(polls$edu)

#问题1
mod.polls1 <- glm(bush~black+female,
                data=polls,family=binomial)
summary(mod.polls1)  
anova(mod.polls1,test='Chisq')

#问题2
mod.polls2 <- glmer(bush~black+female+(1|state),
                data=polls,family=binomial)
summary(mod.polls2)  

library(ggplot2)
ggplot(polls, aes(x=black,y=bush, group=female))+
    geom_point()+
    geom_smooth(method = 'lm')+
    facet_wrap(~state)

#问题3,评估random effect的slope效应应该用bootstrap
mod.polls3 <- glmer(bush~black+female+(1+female|state),
                data=polls,family=binomial)
summary(mod.polls3) 
confint(mod.polls3, method='boot',oldNames=F,nsim = 99)

#问题4
anova(mod.polls3)
#check diagnostics (residual properties)
qqplot()
#residual normality
qqPlot(resid(mod.polls3))

#residual homogeneity

## plot the sqrt of the absolute residuals against fitted values
plot(sqrt(abs(resid(mod.polls3)))~ fitted(mod.polls3))
lines(lowess(sqrt(abs(resid(mod.polls3)))~
               fitted(mod.polls3)), col='red')


#random effect normality
qqPlot(ranef(mod.polls3)$state$'(Intercept)')

#练习pest

pests <- read.csv("exercises/Pesticide.csv")

pests$station <- as.factor(pests$station)
pests$trt<- as.factor(pests$trt)

ggplot(pests, aes(x=trt,y=N, col=trt))+
    geom_point()+
    geom_smooth(method = 'lm')+
    facet_wrap(~station)



hist(pests$N)

mod.pests1 <- glmer(N~trt+(1|station),pests,family = poisson())

sum(resid(mod.pests1,type='pearson')^2/df.residual(mod.pests1))
#结果是over dispersion
#把Nseeds来log一下，因为hist检查是非正态，再检查还是over
#try OLRE，就是新建一个随机效应（正态），为pest#indx

mod.pests2 <- glmer(N~trt+(1+trt|station),pests,family = poisson())
summary(mod.pests2)
anova(mod.pests2)
sum(resid(mod.pests2,type='pearson')^2/df.residual(mod.pests2))
#结果是非over dispersion

#问题7
seedling <- read.csv("exercises/plantdamage2.csv")

seedling$shadehouse <- as.factor(seedling$shadehouse)
seedling$light<- as.factor(seedling$light)


mod.seedling <- glmer(cbind(survs, 4-survs)~damage*light+(1|shadehouse),
                data=seedling,family=binomial)
summary(mod.seedling)  
#因为结果都不显著但是lightL的p值相对较小所以，
mod.seedling2 <- glmer(cbind(survs, 4-survs)~light+(1|shadehouse),
                data=seedling,family=binomial)
summary(mod.seedling2)  
发现light L 的p值非常小了，可以做一下bootsrap发现不夸零，是显著有影响的

library(lme4)
library(ggplot2)
library(lattice)
library(arm)
library(utils)
library(car)

preddat <- expand.grid(light=c('D', 'L'),
                       damage=c(0, 0.1, 0.25),
                       shadehouse =1:10)
summary(preddat)
dim(preddat)
preddat$shadehouse <-as.factor(preddat$shadehouse)
                                      
preddat$pred.fix <- 
  predict(mod.seedling, newdata=preddat, re.form=~0) 
subset(preddat, light=='D' & damage==0) 
preddat$pred.sh <- predict(mod.seedling, newdata=preddat, re.form=NULL) 
subset(preddat, light=='D' & damage==0) ## different between shadehouses
ranef(mod.seedling)$shadehouse[1:10,1]  ## compare this with model random BLUPS


## first, for fixed effects only.

## get out the model formula
form <- formula(mod.seedling, fixed.only=T)
form
form <- update(form, NULL~.) ## removes the response part
form

#we already have a prediction data matrix
head(preddat)

## make the model matrix (we will use our preddat matrix levels)
predmat <- model.matrix(form, data=preddat)
head(predmat)
dim(predmat)


preddat$pred.hand <- as.vector(predmat %*% fixef(mod.seedling))

ranefmat <- model.matrix(~0 + shadehouse, preddat)
dim(ranefmat) ## need 10 columns!
colnames(ranefmat)
head(ranefmat,20)

preddat$pred.hr <- predmat %*% fixef(mod.seedling) +
  ranefmat %*% (ranef(mod.seedling)$shadehouse[,1])

head(preddat,n=10) ## manual and automated versions agree
dim(preddat)


vcv <- vcov(mod.seedling)
vcv
## notice that the standard errors of the fixed effects are
## equal to the sqrt of the diagonal of the vcv
sqrt(diag(vcov(mod.seedling)))
summary(mod.seedling)$coef

semod <-  sqrt(diag(predmat%*%vcv%*%t(predmat)))
semod
dim(seedling)

preddat$ucl <- preddat$pred.fix + semod*1.98
preddat$lcl <- preddat$pred.fix - semod*1.98

ggplot(data=preddat, aes(x=damage, y=pred.fix, 
ymin=lcl, ymax=ucl, colour=light)) +
geom_smooth(stat='identity') + theme_classic()

ggplot(data=preddat, aes(x=damage, y=pred.fix, ymin=lcl, ymax=ucl, colour=light)) +geom_smooth(stat='identity') + facet_wrap(~shadehouse)

preddat$ucl.sh <- preddat$pred.sh + semod*1.98
preddat$lcl.sh <- preddat$pred.sh - semod*1.98
ggplot(data=preddat, aes(x=damage, y=pred.sh, ymin=lcl.sh, ymax=ucl.sh, colour=light)) +geom_smooth(stat='identity') +facet_wrap(~shadehouse) +theme_bw()


VarCorr(mod.seedling)  
summary(mod.seedling)

## extract the variance using
vc1 <- as.data.frame(VarCorr(mod.seedling))
vc1
var.sh <- vc1[vc1$grp=='shadehouse', 'vcov']
var.sh
## now recalculate the std. errors
semod <-  sqrt(diag(predmat%*%vcv%*%t(predmat)) + var.sh)

semod

preddat$ucl.fix <- preddat$pred.fix + semod*1.98
preddat$lcl.fix <- preddat$pred.fix - semod*1.98

head(preddat,n=15)   #sequence repeats across greenhouses


## now we can plot the correct predictons and confidence intervals 

ggplot(data=preddat, aes(x=damage, y=pred.fix, ymin=lcl.fix, ymax=ucl.fix, colour=light)) +geom_smooth(stat='identity') +theme_classic()


## compare to earlier result for known greenhouses

quartz() ##windows()  in microsoft OS computer
ggplot(data=preddat, aes(x=damage, y=pred.fix, ymin=lcl, ymax=ucl, colour=light)) +geom_smooth(stat='identity') + theme_classic()


















  














