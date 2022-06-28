rm(list=ls())
options(digits=3, width=60)


##########################################
## Anthony Ives & Kyle Tomlinson
## 21/3/2019
##########################################

## Prac 9: Phylogenetic regression

#you will need a bunch of R packages to do this prac

library(ape)
library(phylolm)
library(lme4)
library(nlme)
library(phytools)
#library(phyr)
#library(picante)
library(caper)


####################################################
##  Code 9.1 - random evolution of a trait up a tree
####################################################

#make a phylogeny with 15 tips (i.e. 15 species)
n <- 15
phy <- compute.brlen(rtree(n=n), method = "Grafen", power = 1)

par(mfrow=c(1,1))
plot(phy)

# This function simulates some trait data for the species in the phylogeny evolving according to a given model of evolution (in this case Brownian motion - random drift model) (we will discuss brownian motion further later on)
x <- rTraitCont(phy, model = "BM", sigma = 1)


par(mfrow=c(1,2))
contMap(phy, x=x)
par(mai=c(.8,.1,.1,.1))
plot(x,1:length(x))

# it is clear that trait values can evolve differently across the tree purely by chance
# it is clear that closely related species end up with similar trait values because of their shared evolutionary history



####################################################
##  Code 9.2 - Turning a phylogeny into a covariance matrix
####################################################

#recall what our phylogeny, 'phy', looks like:

par(mfrow=c(1,2))
plot(phy)

# construct the covariance matrix of shared path lengths
TreeCovar <- vcv.phylo(phy) 
TreeCovar




####################################################
##  Code 9.3 - Modifying phylogenetic signal
####################################################

# Phylogenetic signal and branch-length transformations
# Pagel's branch-length transformation

# lets modify the branch lengths in 'phy' to impose a different evolutionary history
# we start using Pagel's lambda, which adjusts branch lengths
## Pagel's lambda is bound in [0,1] where lambda ->1 implies strong phylo signal and lambda -> 0 means no phylogenetic signal
# we will make 3 trees

phy.lam.1.0 <- transf.branch.lengths(phy=phy, model="lambda", parameters=list(lambda = 1.0))$tree

phy.lam.0.5 <- transf.branch.lengths(phy=phy, model="lambda", parameters=list(lambda = 0.5))$tree

phy.lam.0.0 <- transf.branch.lengths(phy=phy, model="lambda", parameters=list(lambda = 0.0))$tree

# plot them
par(mfrow=c(1,4))
plot(phy); plot(phy.lam.1.0); plot(phy.lam.0.5); plot(phy.lam.0.0)

# Pagel's lambda shifts the phylogenetic signal from the full phylogenetic tree structure (lambda=1) towards deeper and deeper separation and finally no covariance at all (lamda = 0) as indicated by zero shared path length. This last phylogeny we call the star phylogeny. It indicates no phylogenetic signal

# Now look at the covariance matrices constructed from these trees
# notice what happens to the off-diagonal elements

TreeCovar <- vcv.phylo(phy)
TreeCovar.1.0 <- vcv.phylo(phy.lam.1.0)
TreeCovar.0.5 <- vcv.phylo(phy.lam.0.5)
TreeCovar.0.0 <- vcv.phylo(phy.lam.0.0)

TreeCovar
TreeCovar.1.0
TreeCovar.0.5
TreeCovar.0.0   # this is an Identity matrix



####################################################
##  Code 9.4 - running gls()
####################################################

#recall what our phylogeny, 'phy', looks like, and also the associated trait we made in the first part:

# just a reminder of what this looks like
quartz()
par(mfrow=c(1,2))
contMap(phy, x=x)
par(mai=c(.8,.1,.1,.1))
plot(x,1:length(x))

#now make a dataframe of the trait data (this is not actually necessary..)
d <- data.frame(x)

#  now model x using ordinary lm and gls without a covariance structure

lm1 <- lm(x~1,data=d)
summary(lm1)

gls1 <- gls(x~1, d)
summary(gls1) # identical

# now add in the phylogeny with pagel's lambda = 0.5 
gls2 <- gls(x~1, d, correlation = corPagel(0.5 , phy = phy,fixed=TRUE))
summary(gls2)

# ok lets impose the covariance matrix 'by hand'

#first, an identity matrix
CovarI <- diag(15)

gls3 <- gls(x~1, d, correlation = corSymm(CovarI[lower.tri(CovarI)], fixed=TRUE))
summary(gls3)  #compare to lm1 -> identical!

#second, the covariance matrix for the tree we made with Pagel's lambda = 0.0 (TreeCovar.0.0)
gls4 <- gls(x~1, d, correlation = corSymm(TreeCovar.0.0[lower.tri(TreeCovar.0.0)],fixed=TRUE)) 
summary(gls4)  #compare to lm1 -> identical!

#third, the covariance matrix for the tree we made with Pagel's lambda = 0.5 (TreeCovar.0.5)
gls5 <- gls(x~1, d, correlation = corSymm(TreeCovar.0.5[lower.tri(TreeCovar.0.5)],fixed=TRUE)) 
summary(gls5)  #compare to gls2 -> identical
# Aside: the command fixed=TRUE is actually essential to ensure this; otherwise corPagel is set up to search for a 'best' value for lambda




####################################################
##  Code 9.5 - Hierarchical data: lmm vs gls
####################################################


# first lets make hierarchical data
d <- data.frame(site=rep(1:4, each=6), plot=rep(1:6, times=4), x=0)
dim(d)
d
d$site <- as.factor(d$site)

#here it should be obvious that the data is grouped into 4 sites (random effect)

#now lets make some modelled response values
b0 <- 0
sd.b <- 1   #across sites variation
sd.e <- .5  #within sites variation

for(i in d$site){
	dd <- d[d$site == i,]
	nn <- nrow(dd)
	#introduce between site variation
	b1.site <- rnorm(n=1, mean=0, sd=sd.b) 
	#introduce within site variation 
	d$x[d$site == i] <- b0 + b1.site + rnorm(n=nn, sd=sd.e)  
}

summary(d)
dim(d)
d
par(mfrow=c(1,1))
plot(x ~ site, d, xlab="site", ylab="x")
#here it is obvious that the different sites have different mean values, so there is significant structure in the data caused by different sites, hence significant random effect

#but we can also set up a covariance matrix that represents these relations
#which can then be turned into a phylogeny

#first define the grouping structure
vcv.1 <- kronecker(sd.b^2*diag(nrow=4), matrix(1, nrow=6, ncol=6))
#second impose within sites covariance as the maximum covariance across the main  diagonal (the covariance between an site and itself)
vcv <- vcv.1 + sd.e^2*diag(dim(d)[1])  
# this is our covariance matrix describing the hierarchy in the data

#now turn this into a phylogeny 
phy.lm <- vcv2phylo(vcv, tolerance = 1e-3)
phy.lm$tip.label <- 1:24
phy.lm <- multi2di(phy.lm) #turns multichotomies into dichotomies with zero branch lengths (necessary for gls to run..)

par(mfrow=c(1,2))
plot(phy.lm)
plot(d$x, 1:length(d$x), xlab="site", ylab="x")


# analyse of LMM and PGLS

z.lm <- lm(x ~ 1, REML=F, data=d)
summary(z.lm)
AIC(z.lm)


z.lmm <- lmer(x ~ 1 + (1|site), REML=F, data=d)
summary(z.lmm)


z.gls <- gls(x ~ 1, correlation = corPagel(1.0 , phy = phy.lm) , data = d, method="ML") 
summary(z.gls)

AIC(z.lm,z.lmm,z.gls)
#importantly please note that lmm and gls have generated near identical AIC values






####################################################
##  Code 9.6 - Model of evolution: BM vs OU
####################################################


#make a phylogeny with 15 tips (i.e. 15 species)
n <- 15
phy <- compute.brlen(rtree(n=n), method = "Grafen", power = 1)

par(mfrow=c(1,2))
plot(phy)

# model brownian motion
x <- rTraitCont(phy, model = "BM", sigma = 1)
# model OU process (stabilizing selection)
y <- rTraitCont(phy, model = "OU", sigma = 1,alpha=5)
#note that stabilisation increases with increasing alpha; consequently alpha = 0 yields BM

par(mfrow=c(2,2))
contMap(phy, x=x)
par(mai=c(.8,.1,.1,.1))
plot(x,1:length(x))
contMap(phy, x=y)
par(mai=c(.8,.1,.1,.1))
plot(y,1:length(y))

#its clear from the image that trait variation is accumulating more slowly with OU than BM (compare the trait ranges)


# an alternative way of imposing an OU process is by transforming the branchlengths of the phylogeny and then imposing BM onto the modified tree

# OU branch-length transformation
alpha <- 5
phy.OU <- transf.branch.lengths(phy=phy, model="OUfixedRoot", parameters=list(alpha = alpha))$tree

par(mfrow=c(1,2))
plot(phy)
plot(phy.OU)

#please note what this implies: it says that an OU process reduces the importance of evolutionary history in explaining present day trait values

# what is useful here is that we can use a branchlength transform to impose a OU process onto a tree 




####################################################
##  Code 9.7 - Testing for phylogenetic signal
####################################################

# Tests for phylogenetic signal
n <- 50
phy <- compute.brlen(rtree(n=n), method = "Grafen", power = 1)
plot(phy)


# Simulate data with Pagel's lambda transform
lam <- .5   #the point being that we want phylo signal
phy.lam <- transf.branch.lengths(phy=phy, model="lambda", parameters=list(lambda = lam))$tree
plot(phy.lam)
x <- rTraitCont(phy.lam, model = "BM", sigma = 1)

# Estimate phylogenetic signal using maximum likelihood
# estimate separately using pagel's lambda and using an OU process
z.lam <- phylolm(x ~ 1, phy=phy, model = "lambda")
z.OU <- phylolm(x ~ 1, phy=phy, model = "OUfixedRoot")
summary(z.lam)
summary(z.OU)   
# of course these estimates will be unique to your simulated dataset
AIC(z.lam); AIC(z.OU)  #compare which method gives a better fit to the data




# Method 1. Likelihood Ratio Test (LRT)
z.0 <- lm(x ~ 1)  #null model with no phylo signal/ star phylogeny
AIC(z.0)

LR.lam <- 2*(z.lam$logLik - logLik(z.0)[1])
LR.OU <- 2*(z.OU$logLik - logLik(z.0)[1])

c(LR.lam, pchisq(LR.lam, df=1, lower.tail=F))
c(LR.OU, pchisq(LR.OU, df=1, lower.tail=F))



# Method 2. Parametric bootstrap over the fitted model
z.lam <- phylolm(x ~ 1, phy=phy, model = "lambda", boot=500, full.matrix = TRUE)
z.OU <- phylolm(x ~ 1, phy=phy, model = "OUfixedRoot", boot=500, full.matrix = TRUE)
summary(z.lam)
summary(z.OU)

par(mfrow=c(2,1), mai=c(1,1,.4,.4))
hist(z.lam$bootstrap[,3], main="lam", xlab="")
lines(z.lam$optpar * c(1,1), c(0,1000), col="red")
hist(z.OU$bootstrap[,3], main="OU", xlab="")
lines(z.OU$optpar * c(1,1), c(0,1000), col="red")



# Method 3: Parametric bootstrap over H0
z.0 <- lm(x ~ 1) #set up null model
summary(z.0)
nboot <- 500
boot <- data.frame(lam=rep(0,nboot), alpha=0)
head(boot,n=20)

x.sim <- x
for(i in 1:nboot){
	x.sim <- simulate(z.0)[[1]]
	names(x.sim) <- names(x)
	z.lam.sim <- phylolm(x.sim ~ 1, phy=phy, model = "lambda")
	z.OU.sim <- phylolm(x.sim ~ 1, phy=phy, model = "OUfixedRoot")
	boot$lam[i] <- z.lam.sim$optpar
	boot$alpha[i] <- z.OU.sim$optpar
}

boot

P.lam <- mean(boot$lam > z.lam$optpar)
P.OU <- mean(boot$alpha < z.OU$optpar)

par(mfrow=c(2,1), mai=c(1,1,.4,.4))
hist(boot$lam, main=paste("P =", round(P.lam, digits=3)), xlab="",xlim=c(0,1))
lines(z.lam$optpar * c(1,1), c(0,nboot), col="red")
hist(boot$alpha, main=paste("P =", round(P.OU, digits=3)), xlab="")
lines(z.OU$optpar * c(1,1), c(0,nboot), col="red")



# Method 4: Permutation test
nperm <- 500
perm <- data.frame(lam=rep(0,nperm), alpha=0)
for(i in 1:nperm){
	x.perm <- x[sample(1:length(x))]
	names(x.perm) <- names(x)
	z.lam.perm <- phylolm(x.perm ~ 1, phy=phy, model = "lambda")
	z.OU.perm <- phylolm(x.perm ~ 1, phy=phy, model = "OUfixedRoot")
	perm$lam[i] <- z.lam.perm$optpar
	perm$alpha[i] <- z.OU.perm$optpar
}

P.lam <- mean(perm$lam > z.lam$optpar)
P.OU <- mean(perm$alpha < z.OU$optpar)

par(mfrow=c(2,1), mai=c(1,1,.4,.4))
hist(perm$lam, main=paste("P =", round(P.lam, digits=3)), xlab="")
lines(z.lam$optpar * c(1,1), c(0,nperm), col="red")
hist(perm$alpha, main=paste("P =", round(P.OU, digits=3)), xlab="")
lines(z.OU$optpar * c(1,1), c(0,nperm), col="red")





####################################################
##  Code 9.9 - Phylogenetic signal in the independent variable
####################################################


# make a new tree
n <- 30
phy <- compute.brlen(rtree(n=n), method = "Grafen", power = 1)

#model phylogenetic signal in the predictor x and residuals e
b1 <- 0
lam.x <- 1
lam.e <- 1
phy.x <- transf.branch.lengths(phy=phy, model="lambda", parameters=list(lambda = lam.x))$tree
phy.e <- transf.branch.lengths(phy=phy, model="lambda", parameters=list(lambda = lam.e))$tree

par(mfrow=c(1,2))
plot(phy)
plot(phy.e, root.edge=T)

x <- rTraitCont(phy.x, model = "BM", sigma = 1)
e <- rTraitCont(phy.e, model = "BM", sigma = 1)

#now use these values to predict y
# it is obvious that y will have strong phylo signal given that both x and e have full phylogenetic signal
x <- x[match(names(e), names(x))]
y <- b1 * x + e
y <- array(y)
rownames(y) <- phy$tip.label
	
z <- phylolm(y ~ x, phy=phy, model = "lambda")
summary(z)

#now resimulate x with a different lambda but dont change y
# here i use the star phylogeny for x (lambda = 0)
phy.x2 <- transf.branch.lengths(phy=phy, model="lambda", parameters=list(lambda = 0.0))$tree
x2 <- rTraitCont(phy.x2, model = "BM", sigma = 1)

plot(x,x2)  #so x and x2 are quite unrelated, and you already know that they were built with different phylogenetic models

z2 <- phylolm(y ~ x2, phy=phy, model = "lambda")
summary(z2)
summary(z)
# note that although the predictor coefficients changed a lot, the estimated lambda (phylogenetic signal in y) is about the same







####################################################
##  Code 9.9 - The savanna tree seedling dataset
##              H1: RMF is greater among humid species
####################################################

# A. read in the data and phylogeny

meandata <- read.csv("exampledata/Savanna_seedlings.csv", header=TRUE, na.strings="*")
summary(meandata)
dim(meandata)
head(meandata)

tree1 <- read.nexus("exampledata/Savanna_seedling_nexus",tree.names = NULL)
tree1
#so our trait df has 52 species but our phylogeny has 51 species

plot(tree1)
#please note the structure of the tree. The brnches have different lengths
# this is called an Additive tree, which is based on the raw mutation rates estimated in each branch of the phylogeny by the tree construction method

#we can deal with additive trees, but typically pgls methods use an ultrametric tree where the branches have been scaled to be equal length
#the following function can do this using a parameter called lambda. THIS IS NOT PAGEL's LAMBDA
#deciding on which transformation is appropriate is a whole other discussion which i would like to avoid for now (this is longhand for i dont understand this well yet ;-))
#typically I also run models with additive trees
treeU <- chronopl(tree1, lambda=0.1)

par(mfrow=c(1,2), mai=c(.8,.1,.1,.1))
plot(tree1)
plot(treeU)

# B. reduce both data types to the minmum set of common species

#now link the data to the tree using the caper function 'comparative.data'
#the trick here is that the species names must be the same in the tree and the df and also the species column in the df must be named "Species"
combine<-comparative.data(treeU,meandata,Species,vcv=TRUE,na.omit=F)
combine
# note that now there are only 51 values, meaning that the command has removed the non-matched species

# NOTE: Here i will do a divergence
#I can actually analyse data using the comparative.data object in the function pgls()
#but i can also simply extract the reduced phylogeny and reduced dataframe from this object and use them independently in the way i did in previous analyses

phy.min <- combine$phy
phy.min
plot(phy.min)

meandata2 <- combine$data
summary(meandata2)
dim(meandata2)


####################################################
## FIRST in pgls()

# C. test for phylo signal 

# the trait i am interested in is root mass fraction
# here i am using pagel's lambda 
mod0 <- pgls(RMF  ~ 1, combine, lambda='ML')
summary(mod0)
#technically this model suggests that lambda is not different from zero (because it overlaps zero) 
#but we will proceed using lambda for the sake of completeness

# D. run full model with specified fixed effects

#these are plant traits on seedlings, which tend to change a lot wth plant size, so log plant size needs to be included as a covariate in the model (why log??)

mod1 <- pgls(RMF  ~ log(Mass20) + Continent + Climate, combine, lambda='ML')
summary(mod1)


#then you could test for significant effects in the usual ways


####################################################
## SECOND using the separate reduced phylogeny and dataframe
#rownames(meandata2) <- phy.min$tip.label
meandata2

#C. test for phylogenetic signal

#ignore the warning
z.lam <- phylolm(meandata2$RMF ~ 1, phy=phy.min, model = "lambda")
summary(z.lam)  #same result as above

z.OU <- phylolm(meandata2$RMF ~ 1, phy=phy.min, model = "OUfixedRoot")
summary(z.OU) 


#if you dont trust the ML solution you can then proceed to run permutation tests


x<-meandata2$RMF
phy <- phy.min
nperm <- 500
perm <- data.frame(lam=rep(0,nperm), alpha=0)
for(i in 1:nperm){
	x.perm <- x[sample(1:length(x))]
	names(x.perm) <- names(x)
	z.lam.perm <- phylolm(x.perm ~ 1, phy=phy, model = "lambda")
	z.OU.perm <- phylolm(x.perm ~ 1, phy=phy, model = "OUfixedRoot")
	perm$lam[i] <- z.lam.perm$optpar
	perm$alpha[i] <- z.OU.perm$optpar
}

P.lam <- mean(perm$lam > z.lam$optpar)
P.OU <- mean(perm$alpha < z.OU$optpar)

par(mfrow=c(2,1), mai=c(1,1,.4,.4))
hist(perm$lam, main=paste("P =", round(P.lam, digits=3)), xlab="")
lines(z.lam$optpar * c(1,1), c(0,nperm), col="red")
hist(perm$alpha, main=paste("P =", round(P.OU, digits=3)), xlab="")
lines(z.OU$optpar * c(1,1), c(0,nperm), col="red")

#this indicates that lambda is significant whereas the OU transform is not significant
#as a matter of caution proceed assuming phylogenetic signal


# D. run full model with specified fixed effects

# Now add in the full fixed effects model to test (bX)
z.lam.1 <- phylolm(RMF ~ log(Mass20) + Continent + Climate, data = meandata2, phy=phy.min, model = "lambda")
summary(z.lam.1)  #same result as above; RMF is borderline significant


#parametric bootstrap
z.lam.1.boot <- phylolm(RMF ~ log(Mass20) + Continent + Climate, data = meandata2, phy=phy.min, model = "lambda",boot=500)

confint(z.lam.1.boot)
#so RMF lower in SA species

#what are the coefficients?
coef(z.lam.1)
#coefficient for ClimateSA
coef(z.lam.1)[5]

#here i will permute, but only for lambda models (where significant phylogenetic signal was detected)


x<-meandata2$RMF
phy <- phy.min
nperm <- 500
# we will record permuted coefficients for all predictors and for lambda
perm <- data.frame(lam=rep(0,nperm), Intercept=0,lnMass20=0,ContinentAUS=0,ContinentSAM=0,ClimateSA=0)
for(i in 1:nperm){
	x.perm <- x[sample(1:length(x))]
	names(x.perm) <- names(x)
	z.lam.perm <- phylolm(x.perm ~ log(Mass20) + Continent + Climate, data = meandata2, phy=phy, model = "lambda")
	
	perm$lam[i] <- z.lam.perm$optpar
	perm$Intercept[i] <- coef(z.lam.perm)[1]
	perm$lnMass20[i] <- coef(z.lam.perm)[2]
	perm$ContinentAUS[i] <- coef(z.lam.perm)[3]
	perm$ContinentSAM[i] <- coef(z.lam.perm)[4]
	perm$ClimateSA[i] <- coef(z.lam.perm)[5]
}


#so lets just check root mass climate type
#careful now: coef(z.lam.1)[5] was negative so you want to find probability of getting a value smaller than this number
P.ClimateSA <- mean(perm$ClimateSA < coef(z.lam.1)[5])

hist(perm$ClimateSA, main=paste("P =", round(P.ClimateSA, digits=3)), xlab="")
lines(coef(z.lam.1)[5]* c(1,1),c(0,nperm),  col="red")

#to really understand a permutation result you should model a range of different permutation numbers and only trust the result once the probabilities have stabilised
#if you use your noggin you can work out a loop to do this ;-)


#######  END OF PRAC ######










