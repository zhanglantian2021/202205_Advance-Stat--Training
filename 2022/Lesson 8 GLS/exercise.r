setwd("E:/academic_resources/Advance Stat Training/2022/Lesson 8 GLS")

plants <- read.csv('exampledata/plantdamage.csv')
plants$light <- as.factor(plants$light)
plants$shadehouse <- as.factor(plants$shadehouse)

library(lme4)
library(nlme)

mod1 <- lmer(growth~light*damage + (1|shadehouse),data=plants)
summary(mod1)
plot(resid(mod1))
plot(mod1,resid(.)~fitted(.)|light,abline=0,layout=c(2,1))


mod2 <- lme(growth~light*damage, random = ~1|shadehouse,data=plants)
summary(mod2)
plot(mod2,resid(.)~fitted(.)|light,abline=0,layout=c(2,1))

plot(mod3, abs(resid(., type='p'))~fitted(.), type=c('p', 'r')) ## better






trees <- read.csv('exampledata/trees.csv')

trees$sp <- as.factor(trees$sp)
hist(trees$RGR)
mod1 <- gls(RGR~sp*li,data=trees)

plot(mod1)

plot(Variogram(mod1, form=~x+y, resType='n')) 

mod2 <- update(mod1, correlation=corExp(form=~ x + y))

mod3 <- update(mod1, correlation=corExp(form=~ x + y,nugget=T))


plot(Variogram(mod3, form=~x+y, resType='n')) 






