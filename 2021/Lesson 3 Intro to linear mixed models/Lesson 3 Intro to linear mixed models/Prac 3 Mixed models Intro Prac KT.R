# PRACTICAL 3: Intro to Mixed models

#here I shall very basically introduce you to mixed models

#you will need to install and load the package lme4, a package designed specifically to handle mixed effects problems

library(lme4)


# we shall start with a problem that we have already handled in ANOVA to show that the mixed model function lmer() produces the same result as the aov() function that you have already seen. In this example we will also see that lmer() is able to cope with missing values whereas aov() does not.

# then we shall move onto a problem with model selection





#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#EXAMPLE 3.1. REPEATED MEASURES PROBLEM for Dogs

#recall the problem:
#16 dogs (DogID) used to consider effects of two Drug treatments (morphine vs trimethaphan) (Drug) on blood histamine levels (Yhist) in two types of dogs (dogs with high vs low histamine levels) (DogType)
#histamine levels measured on four occasions (0,1,3,5 minutes after drug administration) (Tme)

#therefore total replication = 16 dogs * 4 time measurements = 64


#so here we make the same dummy dataset

DogID = c(rep(1:16, each = 4))
Drug = c(rep(rep(1:2,each = 4),8))
DogType = c(rep(1:2,each = 32))
Time = c(rep(1:4,16))
Yhist = c(12,3,45,34,45,67,78,34,56,45,67,89,43,52,65,76,98,13,24,16,78,75,48,24,34,12,46,23,46,18,43,52,65,76,98,13,24,16,78,12,3,45,34,45,67,78,34,56,45,67,89,43,52,65,76,98,13,24,16,78,75,48,24,34)


df1 <- data.frame(DogID,Drug,DogType,Time,Yhist)


df1$DogID <- as.factor(df1$DogID)
df1$Drug <- as.factor(df1$Drug)
df1$DogType <- as.factor(df1$DogType)
df1$Time <- as.factor(df1$Time)

summary(df1)

#we can handle this with ANOVA
mod1 <- aov(Yhist~Drug*DogType*Time+Error(DogID/Time),df1)
summary(mod1)


#it can also be handled in lme4. Here we simply need to specify the structure of the random effects in the lmer() function using the +(1|random) specification
#in the present problem, this is simply DogID. lmer() will figure out from the datastructure that time is nested inside DogID whereas the other two variables are not

MM1 <- lmer(Yhist~Drug*DogType*Time + (1|DogID),df1)
summary(MM1)

anova(MM1)

summary(mod1)

#check the F-values (not the probabilities) for each estimate under each model run. You will see that these are near identical for both models. So they produce the same results



#NOW lets see what happens when we introduce some missing values into the response data
#please note that I introduce missing values using "NA", which basically means 'not available'
Yhist = c(12,3,45,34,45,67,78,34,56,45,67,89,43,52,65,76,98,13,24,16,78,75,48,56,34,12,46,23,46,18,43,52,65,76,98,13,24,16,78,12,3,45,NA,45,67,78,34,56,45,67,89,43,52,65,76,98,13,24,16,78,75,48,24,34)


df2 <- data.frame(DogID,Drug,DogType,Time,Yhist)

df2$DogID <- as.factor(df2$DogID)
df2$Drug <- as.factor(df2$Drug)
df2$DogType <- as.factor(df2$DogType)
df2$Time <- as.factor(df2$Time)

summary(df2)

#try ANOVA:

mod2 <- aov(Yhist~Drug*DogType*Time+Error(DogID/Time),df2)

# ANOVA breaks down because of one missing value!!!!
#in the old days we were forced to estimate the missing value to complete the data. now we side-step this problem using mixed models

#try lmer()

MM2 <- lmer(Yhist~Drug*DogType*Time + (1|DogID),df2)
summary(MM2)

#it works. and the coefficients havent changed that much
anova(MM2)

anova(MM1)

#lmer is using numerical approximation to find estimates for the maximum likelihood by trying different parameter values (betas and sigmas)



#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# END OF PRACTICAL
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





