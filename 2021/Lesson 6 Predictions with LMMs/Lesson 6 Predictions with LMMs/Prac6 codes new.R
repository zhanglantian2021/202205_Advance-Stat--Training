
#############################
##Prediction with LMMs
#############################

## CODE from Robert Bagchi (University of Connecticut)


rm(list=ls())
options(digits=3, width=60)


## load the relevant libraries
library(lme4)
library(ggplot2)
library(lattice)
library(arm)
library(utils)
library(car)


#########################
## Code 6.1: Predictions from lm models
#########################

## let's use the irrigation data to do some predictions


irr <- read.csv("exampledata/irrigation.csv")

summary(irr)
head(irr)
irr$soil <- as.factor(irr$soil)


#First, let us fit a model where growth is a function of soil type and water added
mod.irr <- lm(growth ~ soil*water, data=irr)

par(mfrow=c(2,2)); plot(mod.irr)  #ok
summary(mod.irr) #everything significant


#now let us get out the predicted values (fitted values) for the data points that we used in the dataset
fit1 <- fitted(mod.irr)
head(fit1)
irr$fitted <- fitted(mod.irr)
head(irr)

plot(growth~fitted,irr)


## compare this to the observed data
par(mfrow=c(1,1))
plot(fit1, irr$growth)
abline(a=0,b=1, col=2) ## 1:1 line



## now let us make a prediction in a new situation

## Say we want to 
##predict growth 
## after adding 2 units of water
## when the soil is loam
newdat <- data.frame(soil=factor('sand'), water=2)
newdat

## Can predict at this new case with
predict(mod.irr, newdata=newdat) ## gives the expected growth rate

## we can now compute the standard error for this estimate with
predict(mod.irr, newdata=newdat, se.fit=TRUE)

#95% CI: x or y-hat +- 1.96*s.e.

0.388 - 1.96*0.049
0.388 + 1.96*0.049

## If we want confidence intervals (95%)
predict(mod.irr, newdata=newdat, interval='confidence', level=0.95)

## For prediction intervals
predict(mod.irr, newdata=newdat, interval='prediction', level=0.95)

#note that they are not the same!! The 95% prediction interval is more conservative (larger ranges than the confidence interval)

## you might like to read these: 
http://www.graphpad.com/support/faqid/1506/
http://robjhyndman.com/hyndsight/intervals/



##********* BACK TO POWERPOINT *********##


###############################
## Code 6.2 predictions with lmer
##############################


# 1. first, let's try using the predict function method

## use the plant damage data set
plants <- read.csv('exampledata/plantdamage.csv')
summary(plants)
head(plants)

## make 'shadehouse' and 'block' factors so that they can be treated as random effects
plants$light <- as.factor(plants$light)
plants$shadehouse <- as.factor(plants$shadehouse)
plants$block <- as.factor(plants$block)
summary(plants)


## fit a model with only shadehouse as a random effect

mod.pl.lmer1 <- lmer(growth~light*damage + (1|shadehouse),data=plants)
summary(mod.pl.lmer1)


##PREDICT USING THE predict FUNCTION##

## Now to predict some values, you first need to come up with the data that you want to make predictions for

summary(plants)
plants$shadehouse
## All combinations of light and damage names must be EXACTLY the same as for the model
preddat <- expand.grid(light=c('D', 'L'),
                       damage=c(0, 0.1, 0.25),
                       shadehouse =1:10)
summary(preddat)
dim(preddat)
preddat$shadehouse <-as.factor(preddat$shadehouse)
                                      
summary(preddat)
preddat                 
                       


## now compare the levels of factors in original data and
## in the prediction matrix we have constructed
levels(preddat$light)
levels(plants$light)

levels(preddat$shadehouse) 
levels(plants$shadehouse)  

#all looks good

##now we move to use the predict function

## prediction with NO random effects
preddat$pred.fix <- 
  predict(mod.pl.lmer1, newdata=preddat, re.form=~0) 

summary(mod.pl.lmer1 ) 
preddat
head(preddat)  
subset(preddat, light=='D' & damage==0) ## same in each shadehouse

## prediction WITH random effects for shadehouse (NOTE re.form)
preddat$pred.sh <- predict(mod.pl.lmer1, newdata=preddat, re.form=NULL) 
  
#or, if you have multiple random effects and only wish to consider one
#preddat$pred.sh <- predict(mod.pl.lmer1, newdata=preddat, re.form=~(1|shadehouse))
  
head(preddat,n=12)
  
subset(preddat, light=='D' & damage==0) ## different between shadehouses
ranef(mod.pl.lmer1)$shadehouse[1:10,1]  ## compare this with model random BLUPS
0.261 - 1.69  ## = to the BLUP for shadehouse 1


## what about standard errors?
?predict.merMod
# check Details


##########################################
##********* BACK TO POWERPOINT *********##
##########################################



##2. Making predictions using matrix algebra ('by hand')

#recall model form
summary(mod.pl.lmer1)

## first, for fixed effects only.

## get out the model formula
form <- formula(mod.pl.lmer1, fixed.only=T)
form
form <- update(form, NULL~.) ## removes the response part
form

#we already have a prediction data matrix
head(preddat)

## make the model matrix (we will use our preddat matrix levels)
predmat <- model.matrix(form, data=preddat)
head(predmat)
dim(predmat)

## remember that y -hat= X*beta
## prediction is
preddat$pred.hand <- as.vector(predmat %*% fixef(mod.pl.lmer1))
                                 
head(preddat,n=12) ## matrix pred identical to original modelled predictions

## what about the random effects?
## for a random intercept model add the appropriate effect
## we need to make a random effects matrix (Z matrix) first

ranefmat <- model.matrix(~0 + shadehouse, preddat)
dim(ranefmat) ## need 10 columns!
colnames(ranefmat)
head(ranefmat,20)

##now we can make predictions for random levels using model matrices
#y-hat (random effects) = XB +ZB
preddat$pred.hr <- predmat %*% fixef(mod.pl.lmer1) +
  ranefmat %*% (ranef(mod.pl.lmer1)$shadehouse[,1])

head(preddat,n=10) ## manual and automated versions agree
dim(preddat)



##********* BACK TO POWERPOINT *********##


##############################
## Code 6.3: Confidence intervals
##############################


## first we need to extract the variance-covariance matrix so that we can get the standard errors
summary(mod.pl.lmer1)
vcv <- vcov(mod.pl.lmer1)
vcv
## notice that the standard errors of the fixed effects are
## equal to the sqrt of the diagonal of the vcv
sqrt(diag(vcov(mod.pl.lmer1)))
summary(mod.pl.lmer1)$coef



## then calculate the standard errors for each prediction pt 
## using the model matrix
semod <-  sqrt(diag(predmat%*%vcv%*%t(predmat)))
semod


dim(plants)
#we have a very large number of sampling units, so we can probably assume that the t value is being estimated from a standard t-distribution where t = 1.96
# however fo the sake of learning, our model had 103 data pts, 4 fixed effects and 1 random effect
qt(p=0.025, df = 103-4-1)  #this tells us the t-value we could use here for constructing the confidence intervals
## then we can get the confidence intervals
preddat$ucl <- preddat$pred.fix + semod*1.98
preddat$lcl <- preddat$pred.fix - semod*1.98
head(preddat)


## Now we have our predictions!
## plot them
ggplot(data=preddat, aes(x=damage, y=pred.fix, 
ymin=lcl, ymax=ucl, colour=light)) +
geom_smooth(stat='identity') + theme_classic()

ggplot(data=preddat, aes(x=damage, y=pred.fix, ymin=lcl, ymax=ucl, colour=light)) +geom_smooth(stat='identity') + facet_wrap(~shadehouse)
## identical because we haven't dealt with the random effects
head(preddat)

## so redefine the confidence intervals
preddat$ucl.sh <- preddat$pred.sh + semod*1.98
preddat$lcl.sh <- preddat$pred.sh - semod*1.98

ggplot(data=preddat, aes(x=damage, y=pred.sh, ymin=lcl.sh, ymax=ucl.sh, colour=light)) +geom_smooth(stat='identity') +facet_wrap(~shadehouse) +theme_bw()

## difference among shadehouses can be seen.

## the confidence intervals above only account for uncertainty in the FIXED EFFECTS


## If we want to predict for an UNKNOWN SHADEHOUSE, then there should be more uncertainty...

## how do we add this in?
##  the variation among shadehouses is given by
VarCorr(mod.pl.lmer1)  
summary(mod.pl.lmer1)

## extract the variance using
vc1 <- as.data.frame(VarCorr(mod.pl.lmer1))
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

#quartz() ##
windows()  #in microsoft OS computer
ggplot(data=preddat, aes(x=damage, y=pred.fix, ymin=lcl, ymax=ucl, colour=light)) +geom_smooth(stat='identity') + theme_classic()


## the confidence limits are different; wider for the newer results because we have added additional uncertainty into the problem




##********* BACK TO POWERPOINT *********##
#Excercise 6.1
biodepth <- read.csv("exampledata/biodepth.csv")

summary(biodepth)
head(biodepth)
biodepth$location <- as.factor(biodepth$location)
biodepth$block <- as.factor(biodepth$block)
biodepth$plot <- as.factor(biodepth$plot)

#Model
mod.irr <- lm(growth ~ soil*water, data=irr)