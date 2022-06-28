# PRACTICAL 2: GENERALIZED LINEAR MODELS


#in this practical you will learn how to conduct generalised linear models in R

setwd("E:/academic_resources/Advance Stat Training/2022/Lesson 2 Generalised linear models")



# PLEASE NOTE:
# The main functions we use here are found in the basic package of R installed on your computer
# That said, for a number of things we shall use additional packages that you will need to install. These are mentioned as they arise.

library(ggplot2)


## EXAMPLE 2.1. LOGISTIC REGRESSION WITH BINOMIAL DATA 1: THE SEAL RESPONSE DATA
# A behavioural ecologist  approached 50 nursing seals several times over a 2 day period and recorded the number of aggressive responses. She wanted to check if the age of cubs influences the probability of an aggressive response.

#upload the sealData1.csv file
#and look at its structure

seal <- read.csv('sealData1.csv', h=T)
summary(seal)
head(seal)
dim(seal)


#IMPORTANT: the present data  shows the total number of experimented for a given pupage in 'n.obs' and reports the number of "successes" in 'n.response' , so we are not looking at the raw 0-1 data that was collected, but the actual sums of those individual "experiments" categorised by age. This is not a problem in R. As long as we know the total number of experiments ('n.obs'), and the number of successes ('n.response'), R will will be able to make the calculation, as you will see below.


#Now let's plot the proportion of successes against age

plot(n.response/n.obs~pupage,data=seal)

#here we can see the data is bound between 0 and 1 and appears to decline as pupge increases


#so lets run the generalised linear model to test whether age affects aggression
#for generalised linear models we use the glm() function
#we must specify:
#1. the regression relationship we are interested in (cbind(n.response, n.obs-n.response)~pupage)
#(cbind allows us to use the totals in the analysis. Essentially we are binding the number of successes ("n.response") with the number of failures ("n.obs-n.response") in each age category

#2. the residual distribution of the response (family=binomial)
#3. the source of the data for the analysis (seal)

#we can also specify the link function (logit), but this is not necessary
#because the glm will choose the canonical (standard) functions associated with the distribution
#of the response
#compare the two statements below; they should produce the same results

## fit the binomial model from the first day
mod.seal <- glm(cbind(n.response, n.obs-n.response)~pupage,
                data=seal,family=binomial(link=logit))
summary(mod.seal)            
#cbind？
?cbind
#上下比较
mod.seal1 <- glm(cbind(n.response, n.obs-n.response)~pupage,
                 data=seal,family=binomial)
summary(mod.seal1)



#####################################################
############### RETURN TO LECTURE ###################
#####################################################

#PREDICTION

mod.seal <- glm(cbind(n.response, n.obs-n.response)~pupage,
                data=seal,family=binomial(link=logit))
summary(mod.seal)  


#now we will make a output figures using the derived linear generalised linear model

summary(seal)

## build new data frame to house predictions
##it should cover the same range as the values you have recorded in oberevations
newdata <- data.frame(pupage=1:30)
newdata
summary(mod.seal)
## make the predictions - notice the type='link' argument
preds <- predict(mod.seal, newdata=newdata, se.fit=T)

preds #the command produces a fitted value with standard error for each predictor value in newdata
summary(seal)
dim(seal)
## add CI predictions to new data frame in the transformed range 
# as n = 50 for original data, using z=1.96 to get the 95% CIs is reasonable

newdata$fit_T <- preds$fit
newdata$upr_T <- preds$fit+preds$se.fit*1.96
newdata$lwr_T <- preds$fit-preds$se.fit*1.96

head(newdata)

## plot with prediction in the transformed range
plot(newdata$pupage,newdata$fit_T,type="l")
lines(newdata$pupage,newdata$upr_T,lty=2)
lines(newdata$pupage,newdata$lwr_T,lty=2)


#back-transform them with plogis to get response range
?plogis
head(newdata)
newdata$fit <- plogis(preds$fit)
newdata$upr <- plogis(preds$fit+preds$se.fit*1.96) 
newdata$lwr <- plogis(preds$fit-preds$se.fit*1.96)


head(newdata) #note the extra columns in the dataframe
str(newdata)


## plot raw data and model prediction in the response range
#for this we will use ggplot
ggplot(data=seal, aes(x=pupage, y=n.response/n.obs)) + geom_point()+
geom_smooth(data=newdata,aes(x=pupage, y=fit, ymin=lwr, ymax=upr),
              stat='identity')





#####################################################
############### RETURN TO LECTURE ###################
#####################################################

#INFERENCE

#returning to the seal problem, lets check whether the tested model explains a significant amount of deviance

## fit the tested model: aggression ~ pupage           
mod.seal1 <- glm(cbind(n.response, n.obs-n.response)~pupage, data=seal,family=binomial)
summary(mod.seal1)
# residual and null deviance estimates provided at the bottom for tested and null models respectively
logLik(mod.seal1)

       
       
## fit the null model: aggression ~ 1
mod.seal0 <- glm(cbind(n.response, n.obs-n.response)~1, data=seal,family=binomial)
summary(mod.seal0)
# residual and null deviances are identical
logLik(mod.seal0)


#how to get a saturated model? A saturated model has one parameter for each data point
#we can make a dummy to carry this
seal $x1 <- as.factor(1:nrow(seal))
summary(seal)
mod.sealS <- glm(cbind(n.response, n.obs-n.response)~x1, data=seal,family=binomial)
summary(mod.sealS)
#now you can see the residual deviance for the saturated model is basically zero
logLik(mod.sealS)

#so residual deviance for tested model is:
2*(logLik(mod.sealS)-logLik(mod.seal1))

#so residual deviance for null model is:
2*(logLik(mod.sealS)-logLik(mod.seal0))

#compare with output from summary(mod.seal1)
summary(mod.seal1)

#Now compare the residual and null deviances (ask whether there has been a significant reduction in the deviance caused by adding in the predictor pupage)
anova(mod.seal1,test='Chisq')
#the result shows that there is a very significant effect of age on agression


#if you want more detail, look at this explanation:
# https://stats.stackexchange.com/questions/316763/log-likelihood-function-in-poisson-regression
# https://en.wikipedia.org/wiki/Logistic_regression


##Please note: typically you would do this evaluation step before moving to prediction. However I introduce it the other way around because in that way you can see how the glm model prediction relates to the raw data. 



#####################################################
############### RETURN TO LECTURE ###################
#####################################################


#Example 2.2: Poisson data with overdispersion: Aphid data set

aphid2 <- read.csv('/Users/kyletomlinson/Dropbox/Teaching/XTBG advanced stats 2022/Lectures/Lesson 2 Generalised linear models/AphidData2.csv', h=T)

#check out the data
summary(aphid2) #two treatments and counts of aphids. typical poisson data
dim(aphid2)
head(aphid2)
str(aphid2)

aphid2$trt <- as.factor(aphid2$trt)


mod1 <- glm(n.aphids~trt, data=aphid2, family=poisson)
summary(mod1)

par(mfrow=c(2,2)); plot(mod1) 
#look at the size of the residuals in the bottom left plot. some are huge for the one group (2 is about an acceptable cutoff)



#check for overdispersion formally
chisq <- sum(resid(mod1, type='pearson')^2)
chisq/df.residual(mod1) ## much greater than 1

## significantly so?
1-pchisq(chisq, df.residual(mod1)) ## Very significant


#include library to run the dispersiontest
library(AER)
#from the chisq residuals estimate, we know its overdispersion, so set up the test to  check for this
dispersiontest(mod1,alternative = "greater")  
#significantly overdispersed


#so we need to try something else to impove the s.e.'s (quasi models) or account for the overdispersion using a more complex model


#####################################################

# 1. we could try quasi-poisson 

mod.qp <- glm(n.aphids~trt, data=aphid2, family=quasipoisson(link=log))

coef(mod.qp)
coef(mod1)
#coefficients the same

quartz(); par(mfrow=c(2,2)); plot(mod1) 
#windows people, please use "windows()" in place of "quartz()"
quartz(); par(mfrow=c(2,2)); plot(mod.qp) 
## looks similar to mod1
## only thing that changes is the scaling (bottom right plot)
summary(mod.qp)
summary(mod1)
anova(mod.qp, test='Chisq') 

sum(resid(mod1, type='pearson')^2)/df.residual(mod1)
sum(resid(mod.qp, type='pearson')^2)/df.residual(mod.qp) ## the same

summary(mod1)$dispersion
summary(mod.qp)$dispersion ## different between models

#but look what has happened to the standard errors of the model outputs
coef(summary(mod1))
coef(summary(mod.qp))
#they are much larger for the quasipoisson model, which is good becuase now the significance of the estimates are more conservative

#you can also test the significance of the quasipoisson model using ANOVA
anova(mod.qp, test='F')

#a rather unsatisfactory explanation of why quasipoisson uses F-tests and not Chi-sq tests is provided here:
#https://stat.ethz.ch/R-manual/R-devel/library/stats/html/anova.glm.html

#ok lets make some predictions quickly..
fit=fitted.values(mod.qp)
yrep<-as.data.frame(fit)
yrep$trt<-aphid2$trt
library(ggplot2)
ggplot()+
  geom_boxplot(aes(x = aphid2$trt, y=log(aphid2$n.aphids)))+
  geom_point(aes(x = aphid2$trt, y=log(aphid2$n.aphids)))


#####################################################

#2. we could try to fit a more complex model that accounts for the overdispersion. we will use the negative binomial model and the Maxwell-Conway-Poisson model here

#need special packages for this

#negative binomial
library(MASS)
mod.nb <- glm.nb(n.aphids~trt, data=aphid2)
summary(mod.nb)
quartz(); par(mfrow=c(2,2));plot(mod.nb)
#note the size of standardised residuals in the bottom left plot


#CMP
#library(mpcmp)
library(COMPoissonReg)
#https://cran.nyuad.nyu.edu/web/packages/COMPoissonReg/COMPoissonReg.pdf
mod.CMP <- glm.cmp(n.aphids~trt, data=aphid2)
summary(mod.CMP)

#compare coefficients of models
coef(mod1)
coef(mod.qp)
coef(mod.nb)
coef(mod.CMP)


#compare AIC values of models
AIC(mod1); AIC(mod.qp);AIC(mod.nb); AIC(mod.CMP)



#####################################################
############### RETURN TO LECTURE ###################
#####################################################


#EXERCISE 2.3. BINOMIAL DATA 2: THE GRADUATE ADMISSIONS PROBLEM
# problem taken from: http://www.ats.ucla.edu/stat/r/dae/logit.htm

#A researcher is interested in how variables, such as GRE (Graduate Record Exam scores), GPA (grade point average) and rank (prestige of the undergraduate institution) affect admission into graduate school. The response variable, "admit" is a binary variable (admit/don't admit).

#load the data and look at its structure
mydata <- read.csv("binary.csv",header=T,sep=" ")

## view the first few rows of the data
head(mydata)
summary(mydata)
str(mydata)

#This dataset has a binary response variable called admit. Note that in this problem we are dealing with the individual responses (0,1 data) and not the summed proportions. There are three predictor variables: gre, gpa and rank. We will treat the variables gre and gpa as continuous. The variable rank is a factor that takes on the values 1 through 4. Institutions with a rank of 1 have the highest prestige, while those with a rank of 4 have the lowest. We can get basic descriptives for the entire data set by using summary. 

#NOTE: R is reading Rank as a numeric. We need to convert it to a factor. (In general, ALWAYS check the summary() statement or str() statment of imported data to make sure that your data columns are of the correct type.) 
mydata$rank <- as.factor(mydata$rank)
#recheck the data
summary(mydata)
str(mydata)

#now we run the glm with all parameters included
#please note that now we are not dealing with summed successes and failures but rather with the the individual data, so the response is written directly, rather than through using the cbind() function
#admit是0和1，所以是binoaml data
mylogit <- glm(admit ~ gre+gpa+rank, data = mydata, family = "binomial")
summary(mylogit)

anova(mylogit,test='Chisq')


#now check the residual diagnostics

par(mfrow=c(2,2)); plot(mylogit)
# OK now you can see how crazy residuals can look. recall that now we are dealing with individual values of 0 or 1 not proportions (as in the previous example). so these residuals arise because of differences between the fitted sigmoidal line and the transformed zeros and ones.
#Rather what's important here is whether there is overdispersion or not

# perform an overdispersion test (recall we compare residual chi-sq to the residual degrees of freedom)
chisq <- sum(resid(mylogit, type='pearson')^2)
chisq/df.residual(mylogit) ## close to 1; no overdispersion problem


#evaluating overdispersion in binomial models is non-trivial because of the structure of the data (0's and 1's). one suggested solution is to make a binned plot. a binned plot cuts up the fitted vs residual plot into sections call 'bins'. It then calculates a range of 2 standard errors around 0 for the values in each bin using the residuals in each bin. This should contain 95% of the data points and are represented as shaded areas. The dots represent the average residual size in each bin.  Bins with overdispersed values can then be seen lying outside these shaded areas. Overdispersion is only considered a problem if you have many values lying outside of the range.  

library(arm)
?binnedplot
 x <- predict(mylogit)
 y <- resid(mylogit)
 binnedplot(x,y)


#check whether the fitted model explains a significant amount of variation (deviance) in the response variable

null.logit <- glm(admit ~ 1, data = mydata, family = "binomial")
summary(null.logit)
dim(mydata)

anova(null.logit,mylogit,test="Chisq")

#now check the deviance explained by each parameter in the model
anova(mylogit,test="Chisq")
summary(mylogit)





#####################################################
############### RETURN TO LECTURE ###################
#####################################################

#EXERCISE 2.3: Admission continued

mydata <- read.csv("binary.csv",header=T,sep=" ")
summary(mydata)
mydata$rank <- as.factor(mydata$rank)

#a lof models.. lets cut to the additive model

# 1. we could use likelihood ratio tests following backward model selection from the most complex model to the lest complex
#rank means how  is the school
#full model
glm1 <- glm(admit ~ gre+gpa+rank, data = mydata, family = "binomial")
summary(glm1)

#subset model
glm2 <- glm(admit ~ gpa+rank, data = mydata, family = "binomial")
summary(glm2)

anova(glm1,glm2,test="Chisq")

AIC(glm1,glm2)

#this causes a significant change in the deviance, so keep it. 
#then  proceed to check all of the other terms as well

# 2. all subset models using part of the information theoretic approach

library(MuMIn)

options(na.action = "na.fail") 
#supply same maximum model
glm1 <- glm(admit ~ gre+gpa+rank, data = mydata, family = "binomial")
#run the dredge function, which gets all subset models and ranks them highest to lowest
dredge(glm1)
#models are compared using AICc and are ranked from best (lowest AICc value) to worst (highest AICc value)
#output suggests that the model including all 3 predictors is best


#NOTE: for very simple models (such as the one tested here), backward selection using LRTs may be OK, but for complex models a more thorough comparison of the subset models is required


#####################################################
############### RETURN TO LECTURE ###################
#####################################################


#Optional Exercise

#Exercise 2.4. POISSON REGRESSION: THE CHILDREN EVER BORN DATA
# problem taken from: http://data.princeton.edu/wws509/datasets

#These are the data from Fiji on children ever born. Reference: Little, R. J. A. (1978). Generalized Linear Models for Cross-Classified Data from the WFS. World Fertility Survey Technical Bulletins, Number 5.

#The dataset has 70 rows representing grouped individual data. Each row has entries for:
# 		The cell number (1 to 71, cell 68 has no observations),
# 		marriage duration (1=0-4, 2=5-9, 3=10-14, 4=15-19, 5=20-24, 6=25-29),
#		residence (1=Suva, 2=Urban, 3=Rural),
#		education (1=none, 2=lower primary, 3=upper primary, 4=secondary+),
#		mean number of children ever born (e.g. 0.50),
#		variance of children ever born (e.g. 1.14), and
#		number of women in the cell (e.g. 8).


mydata <- read.csv("ceb.csv",header=T)

mydata
summary(mydata)
head(mydata)
mydata$res <- as.factor(mydata$res)
mydata$dur <- as.factor(mydata$dur)
mydata$educ <- as.factor(mydata$educ)

#So now suppose we think education is a good predictor of the number of children a woman would have over the duration of her life. We can plot the counts for the number of children born to each woman for each education group, using ggplot().
library("ggplot2")
ggplot(mydata, aes(mean, fill = educ)) + geom_histogram(binwidth = 0.5, position = "dodge")

# The graph actually shows you a series of poisson distributions (1 for each education category). Looking at the graph it does seem as though women with secondary education have fewer children than women with less education

#we can test this idea more formally using Poisson regression
#Proceed to analyse the data according to the instructions you were given.

mylogit <- glm(mean ~ dur+res+educ, data = mydata)




############### END OF PRACTICAL ################### 

