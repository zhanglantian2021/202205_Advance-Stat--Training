### comparing performance of models with different underlying distributions


library(lme4)

pest <- read.csv('exercises/Pesticide.csv')
summary(pest)
head(pest, 20)
dim(pest)

pest$station <- as.factor(pest$station)
pest$trt <- as.factor(pest$trt)

mod.pest.glmm1 <- glmer(N~trt + (1|station), family=poisson, data=pest)

#now lets use bootstraps to check confidence intervals for the fixed and random effects

## calculating confidence intervals
confint(mod.pest.glmm1, method='boot',oldNames=FALSE, nsim=200) ## takes a little time

summary(mod.pest.glmm1)


## can add an "observation level random effect"
## (this is a common fix in overdispersed poisson models)
pest$indx <- 1:nrow(pest)
summary(pest)

#now test the poisson model again with the OLRE
mod.pest.glmm3 <- glmer(N~trt + (1|station)+ (1|indx), data=pest,family=poisson(link=log))
summary(mod.pest.glmm3)
logLik(mod.pest.glmm3)

pest$fitted.poisson <- (fitted(mod.pest.glmm3))

summary(pest)

#now test a gamma model
#i tried with the OLRE but the model was weird, so i dropped it
mod.pest.glmm4 <- glmer(N~trt + (1|station), data=pest,family=Gamma(link=log))
summary(mod.pest.glmm4)
logLik(mod.pest.glmm4)

pest$fitted.gamma <- (fitted(mod.pest.glmm4))

summary(pest)




#now do model evaluation

library(MLmetrics)
library(MuMIn)
library(ggplot2)

summary(pest)

#R2 estimates for each model
r.squaredGLMM(mod.pest.glmm3)
r.squaredGLMM(mod.pest.glmm4)

#RMSPE tests，均方根误差（RMSE）”是回归模型两项主要性能指标中的一项。 它度量模型所预测的值与实际值之间的平均差值。 它用来估计模型预测目标值的性能（准确度）。 “均方根误差”的值越小，模型的质量越好
RMSPE(y_pred=pest$fitted.poisson,y_true=pest$N)
RMSPE(y_pred=pest$fitted.gamma,y_true=pest$N)
#so here the 

#plot the models

ggplot(data=pest, aes(x=log(N), y=log(fitted.poisson))) + geom_point(col="red")+
  geom_point(data=pest,aes(x=log(N), y=log(fitted.gamma)),
              stat='identity')


par(mfrow=c(1,2))  
plot(log(fitted.poisson)~log(N),data=pest,col="red")
abline(a=0,b=1,col="black")
plot(log(fitted.gamma)~log(N),data=pest,col="blue")
abline(a=0,b=1,col="black")






