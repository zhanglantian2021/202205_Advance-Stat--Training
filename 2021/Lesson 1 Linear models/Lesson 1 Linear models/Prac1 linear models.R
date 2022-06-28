##Practical 1: Linear regression models

setwd("C:\Users\SHRISTEE\Desktop\Advance Stat Training\Lesson 1 Linear models\Lesson 1 Linear models")

getwd()


#Example 1.1: Comparing a t-test with a linear model for groups

#READ IN THE DATA
# import the data "fishspeed.csv" into R and attach the data
fishspeed <-read.table(file="fishspeed.csv",header=T,row.names=NULL,sep=",")
summary(fishspeed)
dim(fishspeed)

fishspeed$Species <- as.factor(fishspeed$Species)
# run summary again
summary(lm(Speed~Temperature,data= fishspeed))
summary(lm(Speed~Species,data= fishspeed))
t.test(Speed~Species,data= fishspeed,var.equal = TRUE)


#####################################################
############### RETURN TO LECTURE ###################
#####################################################


# Example 1.2: FINDING BETAS USING MATRICES

# we will do this for the fish problem using only temperature as a predictor
# recall that we need to add a matrix tha completes these 
# we need to make a matrix that has data for all our coefficients, including B0, the intercept

dim(fishspeed)[1] #always row comes first
dim(fishspeed)
# we need a column of 1's for the intercept of same length as the number of rows in the fishspeed data
# to get intercept repeat values of 1 n times


b0 <- c(rep(1,dim(fishspeed)[1]))
b1 <- fishspeed$Temperature

X <- data.frame(b0,b1)
X <- as.matrix(X)
X
y <- fishspeed$Speed

#recall B-cap = (X'X)-1 X'y
#we need to make all the pieces

#create transpose
Xt <- t(X)
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

# Et Voila!!




#####################################################

# now do it using Species as the predictor #species is factor so we need to use dummy variable

# vector for b0 does not change
b0 <- c(rep(1,dim(fishspeed)[1]))

# but for a factor you need to use a dummy variate to represent the different cases of fish species
#so lets check how many we have of each fish species
fishspeed$Species
#so first 10 values are trout and next 8 are galaxias
b1 <- c(rep(1,10),rep(0,8))
#defaul 0 for galaxy and 1 for trout

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
B
#find X'y
Xy <- Xt%*%y

#now we have all the pieces
B <- XXinv %*% Xy
B

#compare this result to the lm() result

summary(lm(Speed~Species,data= fishspeed))

fishspeed$Species <- relevel(fishspeed$Species,ref="Trout")

# Et Voila!!


#####################################################


#NOW, please calculate the residual variance for the model:
#    Speed ~ species
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
summary(lm1)
#everything significant; most important result is interaction is significant.
#two categories = 2 lines
#please write out the equations for each species!


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

fishspeed2 <-read.table(file="fishspeed2.csv",header=T,row.names=NULL,sep=",")
summary(fishspeed2)
fishspeed2$Species <- as.factor(fishspeed2$Species)
lm1 <- lm(Speed~Temperature*Species,fishspeed2)
summary(lm1)

fishspeed2$Species <- relevel(fishspeed2$Species,ref="Galaxias")
#full linear model
S = b1 + b2*T + b3[Ten] + b4[Tro] + b5*T[Ten] + b6*T[Tro]

G: S = b1 + b2*T
Ten: S = b1 + b2*T + b3[Ten]  + b5*T[Ten] 
Tro: S = b1 + b2*T + b4[Tro]  + b6*T[Tro]

#rerun the model
#....
#following this it would make sense to combine values of Trout and Tenebrias to get a better common estimate using the combined data
#check the levels order:
levels(fishspeed2$Species)
#change the names by reassignment according to that order!
levels(fishspeed2$Species) <- c("Galaxias", "TT", "TT")
#check summary(fishspeed2) and then rerun the model
fishspeed2$Species <- relevel(fishspeed2$Species, ref = "Galaxias")
S = b1 + b2*T + b3[Ten] + b4[Tro] + b5*T[Ten] + b6*T[Tro] #ull model after change

summary(lm(breaks ~ wool + tension, data = fishspeed2$Species))
#####################################################


#Exercise 1.5: 
# Solve the coefficients in the 3-species model using matrix algebra ;-)
# NOTE: you will need to reset "Species' back to three levels; a simple way is just to read in the data again
# your X matrix will need columns for the interaction terms seen in Exercise 1.4


#our equation is:

levels(fishspeed2$Species)

S = b0 + b1*T + b2[Ten] + b3[Tro] + b4*T[Ten] + b5*T[Tro]

# we will do this for the fish problem using only temperature as a predictor
# recall that we need to add a matrix that completes these 
# we need to make a matrix that has data for all our coefficients, including B0, the intercept

dim(fishspeed)[1]

# we need a column of 1's for the intercept of same length as the number of rows in the fishspeed data
# to get intercept repeat values of 1 n times
summary(fishspeed2)
fishspeed2$Species
table(fishspeed2$Species)
b0 <- c(rep(1,dim(fishspeed2)[1]))
b1 <- fishspeed2$Temperature
D1 <- c(rep(1,10),rep(0,8),rep(0,11)) #trout
D2 <- c(rep(0,10),rep(0,8),rep(1,11))  #tenebri
b2 <- D1
b3 <- D2
b4 <- b1*D1 #trout
b5 <- b1*D2 #tenebri


X <- data.frame(b0,b1,b2,b3,b4,b5)
X <- as.matrix(X)
X
y <- fishspeed2$Speed

dim(X)
length(y)
#recall B-cap = (X'X)-1 X'y
#we need to make all the pieces

#create transpose
Xt <- t(X)
dim(X);dim(Xt)

#multiply Xt by X using MATRIX MULTIPLICATION
XX <- Xt %*% X
dim(XX)

#find the inverse (X'X)-1
XXinv <- solve(XX)
XXinv
dim(XXinv)

#find X'y
Xy <- Xt%*%y

#now we have all the pieces
B <- XXinv %*% Xy
B

#compare this result to the lm() result
fishspeed2$Species <- relevel(fishspeed2, ref = )
summary(lm(Speed~Temperature*Species,data= fishspeed2))






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
plot(sqrt(abs(resids))~lm3$fitted); abline(a = 1.96, b = 0, col = 2)

#Levene test for homogeneity of variance across groups
leveneTest(residuals(lm3), as.factor(fishspeed2$Species))
#not significant, so groups variances are not significantly different






############### END OF PRACTICAL ################### 





