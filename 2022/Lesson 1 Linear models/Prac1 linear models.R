##Practical 1: Linear regression models


setwd("E:/academic_resources/Advance Stat Training/2022/Lesson 1 Linear models")

getwd()



#Example 1.1: Comparing a t-test with a linear model for groups

#READ IN THE DATA
# import the data "fishspeed.csv" into R and attach the data
fishspeed <-read.table(file="fishspeed.csv",header=T,row.names=NULL,sep=",")
summary(fishspeed)
dim(fishspeed)

fishspeed$Species <- as.factor(fishspeed$Species)
summary(fishspeed)

summary(lm(Speed~Temperature,data= fishspeed))
summary(lm(Speed~Species,data= fishspeed))
#
t.test(Speed~Species,data= fishspeed,var.equal = TRUE)


#####################################################
############### RETURN TO LECTURE ###################
#####################################################


# Example 1.2: FINDING BETAS USING MATRICES

# we will do this for the fish problem using only temperature as a predictor
# recall that we need to add a matrix tha completes these 
# we need to make a matrix that has data for all our coefficients, including B0, the intercept

dim(fishspeed)[1]

# we need a column of 1's for the intercept of same length as the number of rows in the fishspeed data
# to get intercept repeat values of 1 n times


b0 <- c(rep(1,dim(fishspeed)[1]))
#告诉R，β0一整列是1
b1 <- fishspeed$Temperature

X <- data.frame(b0,b1)
X <- as.matrix(X)
X
y <- fishspeed$Speed

#recall B-cap = (X'X)-1 X'y
#we need to make all the pieces

#create transpose
Xt <- t(X)
#这个是把X矩阵倒置了
dim(X);dim(Xt)

#multiply Xt by X using MATRIX MULTIPLICATION
XX <- Xt %*% X
dim(XX)

#find the inverse (X'X)-1
XXinv <- solve(XX)
XXinv

#find X'y
Xy <- Xt%*%y

#now we have all the pieces
B <- XXinv %*% Xy
B

#compare this result to the lm() result

summary(lm(Speed~Temperature,data= fishspeed))
coef(lm(Speed~Temperature,data= fishspeed))

# Et Voila!!




#####################################################

# now do it using Species as the predictor
summary(lm(Speed~Species,data= fishspeed))
#intercept is species of galxies, and the sum of int and speciestrout is species pf trout
# vector for b0 does not change
b0 <- c(rep(1,dim(fishspeed)[1]))

# but for a factor you need to use a dummy variate to represent the different cases of fish species
#so lets check how many we have of each fish species
fishspeed$Species
table(fishspeed$Species)
#得到每个物种的样本有多少
#so first 10 values are trout and next 8 are galaxias
b1 <- c(rep(1,10),rep(0,8))
#这是指把1重复10遍代表galaxis有8个样本，同理，tourt有10个样本


X <- data.frame(b0,b1)
X <- as.matrix(X)
X
y <- fishspeed$Speed

#recall B-cap = (X'X)-1 X'y
#we need to make all the pieces

#create transpose
Xt <- t(X)

#multiply Xt by X using MATRIX MULTIPLICATION
XX <- Xt %*% X

#find the inverse (X'X)-1
XXinv <- solve(XX)

#find X'y
Xy <- Xt%*%y

#now we have all the pieces
B <- XXinv %*% Xy
B

#compare this result to the lm() result

summary(lm(Speed~Species,data= fishspeed))

# Et Voila!!


#####################################################


#NOW, please calculate the residual variance for the model:
#    Speed ~ Species
# CAREFUL: you need the first X matrix you made (this changes according to the formulation of the model)
#calculate (y'y - B'Xy)*(1/(n-p))

yt <- t(y)
yy <- yt %*% y
Bt <- t(B)

BXy <- Bt %*% Xy
n= dim(fishspeed)[1]
p=2

res.var <- (yy - BXy)*(1/(n-p))
res.var
sqrt(res.var)   # residual standard error

#compare with lm() result; look at the residual standard error
summary(lm(Speed~Species,data= fishspeed))





#####################################################
############### RETURN TO LECTURE ###################
#####################################################




#Example 1.3: multi-model with factor and variate predictors
 
#A researcher measures the maximum swimming speed of 10 brown trout and 10 Canterbury galaxias at a range of temperatures: 
#Response (Y) = swimming speed
#Factor = species 
#Variate (X1) = temperature
#Dummy variable (D) [not necessary to use, but illustrating how the data is constructed]


#plot the data with separate colours for each fish type
dev.off()
par(mfrow=c(1,1))
plot(fishspeed $Temperature, fishspeed $Speed,col= c("red", "blue")[as.numeric(fishspeed $Species)])

#certainly looks like they behave differently
#test this using the linear model

lm1 <- lm(Speed~Temperature*Species,fishspeed)
#species在后面是因为基于species不同所以温度和速度之间的关系模式有变化，而且species是factor。
summary(lm1)
#everything significant; most important result is interaction is significant.
#two categories = 2 lines
# write out the full model and write out the equations for each species

S = b1 + b2*T + b3[Tro] + b4*T[Tro]

G: S = b1 + b2*T
Tro: S = b1 + b2*T + b3[Tro]  + b4*T[Tro] 

#check anova as well
anova(lm1)

# plot the line of fitted values
plot(fishspeed $Temperature,lm1$fitted,ylab="Fish speed",xlab="Temperature")
# add the scatter plot points
points(fishspeed $Temperature,fishspeed $Speed,pch=3,col="blue")




#####################################################
############### RETURN TO LECTURE ###################
#####################################################



#Exercise 1.4: Problem with three factor levels
#Here we have the same problem as before but now we need to compare the swimming speed of 3 fish species: trout, galaxias, tenebrias


#Please do the following
#1. read in the data (fishspeed2.csv)
#2. plot the data for the 3 species
#3. write out the full linear model with dummies into your solution file 
#4. use lm() to evaluate the model: 
#check significance levels and write out the linear models for each species
#Please note!!!! you need to run lm() twice with different default levels of the fish. so please look up the function relevel(), which you will need to change the default level
#SMALL CLUE: EVERYTHING GETS COMPARED TO THE DEFAULT CASE!!!


#####################################################


#Exercise 1.5: 
# Solve the coefficients in the 3-species model using matrix algebra ;-)
# NOTE: you will need to reset "Species' back to three levels; a simple way is just to read in the data again
# your X matrix will need columns for the interaction terms seen in Exercise 1.4



#####################################################
############### RETURN TO LECTURE ###################
#####################################################

# Example 1.6: checking diagnostics on lm() objects

#use the model you built with fishspeed2 in EXERCISE 1E

fishspeed2 <-read.table(file="fishspeed2.csv",header=T,row.names=NULL,sep=",")
summary(fishspeed2)

lm3 <- lm(Speed~Temperature*Species,fishspeed2)
summary(lm3)

#running diagnostics is simple
#first partition

par(mfrow=c(2,2))
plot(lm3)


## SOME FORMAL TESTS FOR THE RESIDUALS

library(car)
###Assumptions###
#Normality
resids <- resid(lm3, type='pearson')

## qq plot with 95% CIs
par(mfrow=c(1,1))
qqPlot(resids)

#Shapiro test to check for non-normality (a significant test result means non-normal)
shapiro.test(residuals(lm3))
#not significant so residuals are normal

#Homoscedasticity
par(mfrow=c(1,1))
plot(sqrt(abs(resids))~lm3$fitted); abline(a = 1.96, b = 0, col = 2)

#Levene test for homogeneity of variance across groups
leveneTest(residuals(lm3), as.factor(fishspeed2$Species))
#not significant, so groups variances are not significantly different


#more fancy plots
library(tidyverse)
library(performance)
library(see)
library(qqplotr)
pp_check(lm3) %>% plot()
check_model(lm3)



############### END OF PRACTICAL ################### 





