setwd("E:/academic_resources/Advance Stat Training/2022/Lesson 4 LMM types")

biodepth <- read.csv('exercises/Biodepth.csv', h=TRUE)
summary(biodepth)
dim(biodepth)

biodepth$location <- as.factor(biodepth$location)
biodepth$block <- as.factor(biodepth$block)
biodepth$plot <- as.factor(biodepth$plot)

isNested(biodepth$block,biodepth$location)

table(biodepth$location,biodepth$block)


biodepth$block.loc <- paste0(biodepth$location,"_",biodepth$block)

head(biodepth)

isNested(biodepth$block.loc,biodepth$location)

table(biodepth$location,biodepth$block.loc)

#check Normality of response data

quartz()
par(mfrow=c(2,2))
hist(biodepth$biomass)
hist(log(biodepth$biomass+1))
hist((biodepth$biomass)^(1/2))
hist((biodepth$biomass)^(1/3))

biodepth$biomass.cubrt <- (biodepth$biomass)^(1/3)

summary(biodepth)
#make the model
lmm1 <- lmer(biomass.cubrt~log.diversity + (1|location/block.loc),biodepth)
summary(lmm1)
#block.loc explains nothing!

#check diagnostics (residual properties)

#residual normality
par(mfrow=c(1,2))
qqPlot(resid(lmm1))

#residual homogeneity

## plot the sqrt of the absolute residuals against fitted values
plot(sqrt(abs(resid(lmm1)))~ fitted(lmm1))
lines(lowess(sqrt(abs(resid(lmm1)))~
               fitted(lmm1)), col='red')

#random effect normality

qqPlot(ranef(lmm1)$location$'(Intercept)')
qqPlot(ranef(lmm1)$block.loc$'(Intercept)')


#Q2 is asking for a random slopes model
#we know that block.loc explains nothing so lets drop it and in random slopes

lmm2 <- lmer(biomass.cubrt~log.diversity + (1+log.diversity|location),biodepth)
summary(lmm2)

#random intercepts not very important