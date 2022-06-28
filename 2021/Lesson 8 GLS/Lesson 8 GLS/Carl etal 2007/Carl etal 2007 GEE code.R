#%%%%%%%%%%
# Helper functions for the spatial usage of gee and geepack
# Code written by Gudrun Carl, 2005-2007

# for more details on GEE see also: 
# Carl, G. & Kühn, I. 2007. Analyzing spatial autocorrelation in species distributions using Gaussian and logit models. Ecological Modelling, in press.
##############################################################

#########################################################################
dat.nn<-function(data,n){
#########################################################################
# A function to generate clusters and order variables and
# to produce a data frame with response, predictors, coordinates, and 
# 3 new parameters:
# o for order
# id for identifying clusters and 
# waves for identifying members of clusters
#
# Arguments
# data      a data frame with response and predictors and 
#           in the last columns with ordered cartesian coordinates
# n         for maximal cluster size  n*n 
#########################################################################
  
l<-dim(data)[2]
OST<-data[,l-1]
NORD<-data[,l]
ko<-OST-min(OST)
idx<-(ko-(ko%%(n)))/n+1
ks<-NORD-min(NORD)
idy<-(ks-(ks%%(n)))/n+1
ie<-(idy-1)*max(idx)+idx
idwx<-ko%%(n)+1
idwy<-ks%%(n)+1
wav<-(idwy-1)*n+idwx
data<-as.matrix(data)
o<-order(ie,wav)
id<-ie[o]
waves<-wav[o]
dat.new1<-data[o,]
dat.new2<-cbind(dat.new1,o,id,waves)
dat.new<-as.data.frame(dat.new2)
}
 

#########################################################################
a.gee<-function(mgee,n,type="glm",corstr="independence",quad=T) {
#########################################################################
# A function to order correlation parameters of Generalized Estimating 
# Equation Models
# Arguments
# mgee       matrix or vector of correlation parameters according to model
# n          for maximal cluster size n*n
# type       type of model 
#            "glm", "gee", "geese" are allowed
# corstr     correlation structure
#            "independence", "exchangeable", "userdefined" are allowed
# quad       by default quadratic correlation structure
#            for model "geese" and "userdefined" correlation only
#########################################################################

if(n==2)n3<-6
if(n==3)n3<-36
if(n==4)n3<-120
a<-rep(0,n3)
if(type=="glm") a<-a
if(type=="gee"){ 
  if(corstr=="exchangeable") a[c(1:n3)]<-mgee[1,2]
  if(corstr=="independence") a<-a
}
a<-as.vector(a)

if(type=="geese") {
 if(corstr=="userdefined"){
if(quad) {
if(n==2)  {
 a<-rep(0,6)
 a[c(1,2,5,6)]<-mgee[1]
 a[c(3,4)]<-mgee[2]
}
if(n==3)  {
 a<-rep(0,36)
 a[c(1,3,9,11,18,22,24,27,29,33,34,36)]<-mgee[1]
 a[c(2,6,14,21,23,35)]<-mgee[2]
 a[c(4,10,12,17,25,28,30,32)]<-mgee[3]
 a[c(5,7,13,15,16,20,26,31)]<-mgee[4]	
 a[c(8,19)]<-mgee[5]	
}
if(n==4)  {
 a<-rep(0,120)
 a[c(1,4,16,19,30,33,46,55,58,66,69,76,79,88,93,96,100,103,106,109,
114,115,118,120)]<-mgee[1]
 a[c(2,8,17,23,37,50,56,62,67,73,83,92,94,101,116,119)]<-mgee[2]
 a[c(3,12,27,41,54,57,95,117)]<-mgee[3]
 a[c(5,18,20,32,34,45,59,68,70,78,80,87,97,102,104,108,110,113)]<-mgee[4]
 a[c(6,9,21,22,24,31,36,38,44,49,60,63,71,72,74,77,82,84,86,91,98,
105,107,112)]<-mgee[5]
 a[c(7,13,26,28,40,42,43,53,61,85,99,111)]<-mgee[6]
 a[c(10,25,35,48,64,75,81,90)]<-mgee[7]
 a[c(11,14,29,39,47,52,65,89)]<-mgee[8]
 a[c(15,51)]<-mgee[9]
}}
if(!quad) a<-mgee 
}
if(corstr=="exchangeable") a[c(1:n3)]<-mgee
if(corstr=="independence") a<-a
 }
a<-as.vector(a)
}

#########################################################################
clus.sz<-function(id){
#########################################################################
# A function to calculate sizes of clusters
# Argument
# id     vector which identifies the clusters
#########################################################################

clus<-rep(0,length(id))
k0<-0
k1<-1
for(i in 2:length(id)) { i1<-i-1
if(id[i]==id[i1]) {k1<-k1+1
if(i==length(id)) {k0<-k0+1
                  clus[k0]<-k1}}
if(id[i]!=id[i1]) {k0<-k0+1
                  clus[k0]<-k1
                  k1<-1
if(i==length(id)) {k0<-k0+1
                  clus[k0]<-k1 }}}
clusz<-clus[clus>0]
}

#########################################################################  
zcor.quad<-function(zcor,n,quad=TRUE) {
#########################################################################
# A function to create a quadratic correlation structure 
# zcor    an object of class "genZcor" (see: geepack)
# n       for maximal cluster size n*n
# quad    by default quadratic correlation structure
#########################################################################

if(quad) {
if(n==2)  {
 zcorn<-matrix(0,dim(zcor)[1],2)
 zcorn[,1]<-zcor[,1]+zcor[,2]+zcor[,5]+zcor[,6]
 zcorn[,2]<-zcor[,3]+zcor[,4]
}
if(n==3)  {
 zcorn<-matrix(0,dim(zcor)[1],5)
 zcorn[,1]<-zcor[,1]+zcor[,3]+zcor[,9]+zcor[,11]+zcor[,18]+zcor[,22]+
 zcor[,24]+zcor[,27]+zcor[,29]+zcor[,33]+zcor[,34]+zcor[,36]
 zcorn[,2]<-zcor[,2]+zcor[,6]+zcor[,14]+zcor[,21]+zcor[,23]+zcor[,35]
 zcorn[,3]<-zcor[,4]+zcor[,10]+zcor[,12]+zcor[,17]+zcor[,25]+zcor[,28]+
 zcor[,30]+zcor[,32]
 zcorn[,4]<-zcor[,5]+zcor[,7]+zcor[,13]+zcor[,15]+zcor[,16]+zcor[,20]+
 zcor[,26]+zcor[,31]
 zcorn[,5]<-zcor[,8]+zcor[,19]
}
if(n==4)  {
 zcorn<-matrix(0,dim(zcor)[1],9)
  zcorn[,1]<-zcor[,1]+zcor[,4]+zcor[,16]+zcor[,19]+zcor[,30]+zcor[,33]+
zcor[,46]+zcor[,55]+zcor[,58]+zcor[,66]+zcor[,69]+zcor[,76]+
zcor[,79]+zcor[,88]+zcor[,93]+zcor[,96]+zcor[,100]+zcor[,103]+
zcor[,106]+zcor[,109]+zcor[,114]+zcor[,115]+zcor[,118]+zcor[,120]
  zcorn[,2]<-zcor[,2]+zcor[,8]+zcor[,17]+zcor[,23]+zcor[,37]+zcor[,50]+
zcor[,56]+zcor[,62]+zcor[,67]+zcor[,73]+zcor[,83]+zcor[,92]+
zcor[,94]+zcor[,101]+zcor[,116]+zcor[,119]
  zcorn[,3]<-zcor[,3]+zcor[,12]+zcor[,27]+zcor[,41]+zcor[,54]+zcor[,57]+
zcor[,95]+zcor[,117]
  zcorn[,4]<-zcor[,5]+zcor[,18]+zcor[,20]+zcor[,32]+zcor[,34]+zcor[,45]+
zcor[,59]+zcor[,68]+zcor[,70]+zcor[,78]+zcor[,80]+zcor[,87]+
zcor[,97]+zcor[,102]+zcor[,104]+zcor[,108]+zcor[,110]+zcor[,113]
  zcorn[,5]<-zcor[,6]+zcor[,9]+zcor[,21]+zcor[,22]+zcor[,24]+zcor[,31]+
zcor[,36]+zcor[,38]+zcor[,44]+zcor[,49]+zcor[,60]+zcor[,63]+
zcor[,71]+zcor[,72]+zcor[,74]+zcor[,77]+zcor[,82]+zcor[,84]+
zcor[,86]+zcor[,91]+zcor[,98]+zcor[,105]+zcor[,107]+zcor[,112]
  zcorn[,6]<-zcor[,7]+zcor[,13]+zcor[,26]+zcor[,28]+zcor[,40]+zcor[,42]+
zcor[,43]+zcor[,53]+zcor[,61]+zcor[,85]+zcor[,99]+zcor[,111]
  zcorn[,7]<-zcor[,10]+zcor[,25]+zcor[,35]+zcor[,48]+zcor[,64]+zcor[,75]+
zcor[,81]+zcor[,90]
  zcorn[,8]<-zcor[,11]+zcor[,14]+zcor[,29]+zcor[,39]+zcor[,47]+zcor[,52]+
zcor[,65]+zcor[,89]
  zcorn[,9]<-zcor[,15]+zcor[,51]
}
} 
if(!quad) zcorn<-zcor
zcorn<-as.matrix(zcorn)
}

#########################################################################
res.gee<-fu