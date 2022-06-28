
getwd()
pollsdata_day2 <- read.csv('pollsdata_day2.csv', h=T)
summary(pollsdata_day2)
head(pollsdata_day2)
dim(pollsdata_day2)
#define as factor
pollsdata_day2$black<-as.factor(pollsdata_day2$black)
pollsdata_day2$female<-as.factor(pollsdata_day2$female)
pollsdata_day2$edu<-as.factor(pollsdata_day2$edu)
pollsdata_day2$state<-as.factor(pollsdata_day2$state)
#plot th data
plot(bush~ female + black, data=pollsdata_day2)

#Model
mod <- glm(bush~ female+black,
                data=pollsdata_day2,family=binomial(link=logit))
summary(mod) 
str(mod)
library(ggplot2)

#dev.off()
dev.off()
#z value is due to non-normal data
#Always use  cbind for (success and failures) (n.obs-n.response = no. of fails)
mod1 <- glmer(bush~female+black+ (1|state), data=pollsdata_day2,family=binomial(link=logit))
summary(mod1)
#ggplot here
ggplot(data=pollsdata_day2, aes(x=female, y=bush, color=black, group=state)) + geom_point()+
  geom_smooth(method="lm") +
  facet_wrap(~state)

mod2 <- glmer(bush~female+black+ (1+female|state), data=pollsdata_day2,family=binomial(link=logit))
summary(mod2)

#confidence interval
confint(mod2, methods="boot", oldNames=FALSE, nsim=499)

#Pest excercies
Pesticide <- read.csv('Pesticide.csv', h=T)
summary(Pesticide)
head(Pesticide)
dim(Pesticide)

Pesticide$station<-as.factor(Pesticide$station)
Pesticide$trt<-as.factor(Pesticide$trt)
summary(Pesticide)
head(Pesticide)
dim(Pesticide)

mod.pest.glmm1 <- glmer(Nseeds~trt + (1|station), family=poisson, data=Pesticide)

confint(mod.pest.glmm1, method='boot',oldNames=FALSE, nsim=200) ## takes a little time

summary(mod.pest.glmm1)
Pesticide$indx <- 1:nrow(Pesticide)
summary(Pesticide)

#now test the poisson model again with the OLRE
mod.pest.glmm3 <- glmer(N~trt + (1|station)+ (1|indx), data=Pesticide,family=poisson(link=log))
summary(mod.pest.glmm3)
logLik(mod.pest.glmm3)

Pesticide$fitted.poisson <- (fitted(mod.pest.glmm3))

summary(Pesticide)

#now test a gamma model
#i tried with the OLRE but the model was weird, so i dropped it
mod.pest.glmm4 <- glmer(Nseeds~trt + (1|station), data=Pesticide,family=Gamma(link=log))
summary(mod.pest.glmm4)
logLik(mod.pest.glmm4)

Pesticide$fitted.gamma <- (fitted(mod.pest.glmm4))

summary(Pesticide)

#now do model evaluation

library(MLmetrics)
library(MuMIn)
library(ggplot2)

summary(Pesticide)

#R2 estimates for each model
r.squaredGLMM(mod.pest.glmm3)
r.squaredGLMM(mod.pest.glmm4)

#RMSPE tests
RMSPE(y_pred=Pesticide$fitted.poisson,y_true=Pesticide$Nseeds)
RMSPE(y_pred=Pesticide$fitted.gamma,y_true=Pesticide$Nseeds)
#so here the 

#plot the models

ggplot(data=Pesticide, aes(x=log(Nseeds), y=log(fitted.poisson))) + geom_point(col="red")+
  geom_point(data=Pesticide,aes(x=log(Nseeds), y=log(fitted.gamma)),
             stat='identity')


par(mfrow=c(1,2))  
plot(log(fitted.poisson)~log(Nseeds),data=Pesticide,col="red")
abline(a=0,b=1,col="black")
plot(log(fitted.gamma)~log(Nseeds),data=Pesticide,col="blue")
abline(a=0,b=1,col="black")
