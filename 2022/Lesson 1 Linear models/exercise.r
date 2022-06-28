setwd("E:/academic_resources/Advance Stat Training/2022/Lesson 1 Linear models")

fishspeed <-read.table(file="fishspeed2.csv",header=T,row.names=NULL,sep=",")

fishspeed$Species <- as.factor(fishspeed$Species)
summary(fishspeed)
dim(fishspeed)

par(mfrow=c(1,1))
plot(fishspeed $Temperature, fishspeed $Speed,col= c("red", "blue","green")[as.numeric(fishspeed $Species)])

lm1 <- lm(Speed~Temperature*Species,fishspeed)
summary(lm1)

#y(G)=22.44+x；y(Te)=(22.44+21.58)+(1+1.37)x;y(tr)=(22.44+18.35)+(1+1.63)x

#改变default case，本来默认是galaxes，也就是把它编码为1，其他编码为0，是以它为基准的。基准线改变，但是不会改变每个方程的系数，可是会改变结果显示是否显著，因为比较的基准线变了，是指另外两天线分别和基准线的截距和斜率相比，另外两条线的截距和斜率与基准相比的变化的部分是否有显著差异。
?relevel
fishspeed$Species <- relevel(fishspeed$Species, ref = "Tenebrias")
lm2 <- lm(Speed~Temperature*Species,fishspeed)

#full model 应该这么写:S = b0 + b1*T + b2[Ten] + b3[Tro]  + b4*T[Ten]  + b5*T[Tro]



#以矩阵的方式写,是要以上面full model 的形式写的，
#以Tenebrias为基准
b0 <- c(rep(1,dim(fishspeed)[1]))
table(fishspeed$Species)
#在只有两个物种的时候没有把温度赋予到b1是
b1 <- fishspeed$Temperature
b2 <- c(rep(0,8),rep(1,11),rep(0,10))
b3 <- c(rep(0,8),rep(0,11),rep(1,10))

#这是关于相互作用
b4 <- b2*b1
b5 <- b3*b1
X <- data.frame(b0,b1,b2,b3,b4,b5)
X <- as.matrix(X)
X
y <- fishspeed$Speed
Xt <- t(X)
XX <- Xt %*% X
XXinv <- solve(XX)
Xy <- Xt%*%y
B <- XXinv %*% Xy
B

#kyle答案
#Exercise 1.4: Problem with three factor levels
#Here we have the same problem as before but now we need to compare the swimming speed of 3 fish species: trout, galaxias, tenebrias


#Please do the following
#1. read in the data (fishspeed2.csv)
#2. plot the data for the 3 species
#3. write out the full linear model with dummies into your solution file 
#4. use lm() to evaluate the model: 
#check significance levels and write out the linear models for each species
#Please note!!!! you need to run lm() twice with different default levels of the fish. so please look up the function relevel(), which you will need to change the default level
#REMEMBER: EVERYTHING GETS COMPARED TO THE DEFAULT CASE!!!


fishspeed2 <-read.table(file="fishspeed2.csv",header=T,row.names=NULL,sep=",")
summary(fishspeed2)
dim(fishspeed2)

fishspeed2$Species <- as.factor(fishspeed2$Species)

lm1 <- lm(Speed~Temperature*Species,fishspeed2)
summary(lm1)

quartz()
plot(fishspeed2 $Temperature, fishspeed2 $Speed,col= c("red","green", "blue")[as.numeric(fishspeed2 $Species)])

#change default

fishspeed2 $Species <- relevel(fishspeed2 $Species,ref="Tenebrias")

lm1 <- lm(Speed~Temperature*Species,fishspeed2)
summary(lm1)

#with galaxias as default



#####################################################


#Exercise 1.5: 
# Solve the coefficients in the 3-species model using matrix algebra ;-)
# NOTE: you will need to reset "Species' back to three levels; a simple way is just to read in the data again
# your X matrix will need columns for the interaction terms seen in Exercise 1.4

fishspeed2 $Species <- relevel(fishspeed2 $Species,ref="Galaxias")

S = b0 + b1*T + b2[Ten] + b3[Tro]  + b4*T[Ten]  + b5*T[Tro]

B ~ (X'X)-1 X'y

b0 <- c(rep(1,dim(fishspeed2)[1]))
b1 <- fishspeed2$Temperature

#for the factor intercepts
fishspeed2$Species
table(fishspeed2$Species)
b2 <- c(rep(0,10),rep(0,8),rep(1,11))
b3 <- c(rep(1,10),rep(0,8),rep(0,11))

#for the interactions
b4 <- b2*b1 #tenebrias
b5 <- b3*b1 #trout

X <- data.frame(b0,b1,b2,b3,b4,b5)
X <- as.matrix(X)
X
y <- fishspeed2$Speed

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

coef(lm(Speed~Temperature*Species,fishspeed2))