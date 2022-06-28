#Optional Exercise

#Exercise 2.4. POISSON REGRESSION: THE CHILDREN EVER BORN DATA
# problem taken from: http://data.princeton.edu/wws509/datasets

#These are the data from Fiji on children ever born. Reference: Little, R. J. A. (1978). Generalized Linear Models for Cross-Classified Data from the WFS. World Fertility Survey Technical Bulletins, Number 5.

#The dataset has 70 rows representing grouped individual data. Each row has entries for:
#   The cell number (1 to 71, cell 68 has no observations),
#   marriage duration (1=0-4, 2=5-9, 3=10-14, 4=15-19, 5=20-24, 6=25-29),
#  residence (1=Suva, 2=Urban, 3=Rural),
#  education (1=none, 2=lower primary, 3=upper primary, 4=secondary+),
#  mean number of children ever born (e.g. 0.50),
#  variance of children ever born (e.g. 1.14), and
#  number of women in the cell (e.g. 8).


mydata <- read.csv("ceb.csv",header=T)

mydata
summary(mydata)
head(mydata)

mydata$dur <- as.factor(mydata$dur)
mydata$res <- as.factor(mydata$res)
mydata$educ <- as.factor(mydata$educ)

#So now suppose we think education is a good predictor of the number of children a woman would have over the duration of her life. We can plot the counts for the number of children born to each woman for each education group, using ggplot().
library(ggplot2)
ggplot(mydata, aes(mean, fill = educ)) + geom_histogram(binwidth = 0.5, position = "dodge")

# The graph actually shows you a series of poisson distributions (1 for each education category). Looking at the graph it does seem as though women with secondary education have fewer children than women with less education

#we can test this idea more formally using Poisson regression
#Proceed to analyse the data according to the instructions you were given.


#a poisson model?
hist(mydata$mean)
p1 <- glm(mean~ dur + res + educ,mydata,family="poisson")
summary(p1)
anova(p1,test="Chisq")
logLik(p1)

g1 <- glm(mean~ dur + res + educ,mydata,family="Gamma")
summary(g1)
anova(g1,test="Chisq")
logLik(g1)

# perform an overdispersion test (recall we compare residual chi-sq to the residual degrees of freedom)
chisq <- sum(resid(g1, type='pearson')^2)
chisq/df.residual(g1) ## underdispersed

#model selection: all subsets

library(MuMIn)

options(na.action = "na.fail") 
#supply same maximum model
p1 <- glm(mean~ dur + res + educ,mydata,family="poisson")
#run the dredge function, which gets all subset models and ranks them highest to lowest
dredge(p1)
#problem!!

g1 <- glm(mean~ dur + res + educ,mydata,family="Gamma")
#run the dredge function, which gets all subset models and ranks them highest to lowest
dredge(g1)