###################################################
# Author: Vitor Sagae
# Creating design matrices
###################################################

############### Environments (E) ############################
rm(list=ls()) 

setwd("./data/")

Y<-read.csv("./Y.csv",sep=";")
Y<-Y[Y$year!="2011",]
#Y<-Y[!Y$environ=="NE_2012",]
head(Y)

# We assign the column that links phenotypes with environments
colIDy <- 4
IDs <- Y[,colIDy]
Y[,colIDy]<-factor(Y[,colIDy])  

# Must create a matrix 
Z<-as.matrix(model.matrix(~Y[,colIDy]-1)) 
d<-colSums(Z)
V<-Z

for(i in 1:ncol(Z)){ V[,i]<-V[,i]/sqrt(d[i]) } 
EVD.E <- list(vectors=V,values=d)
G.E <- tcrossprod(Z)

save(G.E,file='./G.E.rda')
save(EVD.E,file='./EVD.E.rda')

############### W Matrix (WP) ############################
rm(list=ls())


W<-read.csv("./WP.csv",sep=',',na.strings = "NA")
W<-W[,-1]
W<-as.matrix(W)
Y<-read.csv("./Y.csv",sep=";")
Y<-Y[Y$year!=2011,]
#Y<-Y[!Y$environ=="NE_2012",]

head(Y)
n<-dim(W)[1]
p<-dim(W)[2]

## Imputing, centering, standarizing. 
# Parameters for centering and standarizing

ctr <- TRUE
std <- TRUE

S<-0
for(i in 1:ncol(W)){
  meanWi <- mean(W[,i],na.rm=TRUE)
  # naive imputation
  W[,i] <- ifelse(is.na(W[,i]),meanWi, W[,i]) 
  if(ctr){ W[,i]<-W[,i]-meanWi } #centering
  if(std){ W[,i]<-W[,i]/sd(W[,i]) } # standarizing
  S<-S+1    # var(W[,i])
}

W<-W[,colSums(is.na(W))==0]
#save(W,file="W.matrix.rda")
G.W<-tcrossprod(W)/ncol(W)


save(G.W,file='./G.W.rda')

EVD.W<-eigen(G.W)
rownames(EVD.W$vectors)<-rownames(G.W)

save(EVD.W,file='./EVD.W.rda')
#image(G.W)

############### Genotypes (G) ############################
rm(list=ls()) 

Y<-read.csv("./Y.csv",sep=";")
Y<-Y[Y$year!="2011",]
#Y<-Y[!Y$environ=="NE_2012",]

X<-read.csv("./X.csv",sep=",")
rownames(X)<-X$X
X<-X[,-1]

gennames<-unique(Y$strain)

X<-X[rownames(X) %in% gennames,]

head(Y)  
n<-dim(X)[1]  

# We assign the column that links phenotypes with genotypes
colIDy <-4
IDs <- Y[,colIDy]

colNames <- colnames(X)
IDx <- rownames(X)

p<-length(colNames)

rownames(X)<-IDx
colnames(X)<-colNames

## Imputing, centering, standarizing. 
# Parameters for centering and standarizing
ctr <- TRUE
std <- TRUE

S<-0
for(i in 1:ncol(X)){
  meanXi <- mean(X[,i],na.rm=TRUE)
  # naive imputation
  X[,i] <- ifelse(is.na(X[,i]),meanXi, X[,i]) 
  if(ctr){ X[,i]<-X[,i]-meanXi } #centering
  if(std){ X[,i]<-X[,i]/sd(X[,i]) } # standarizing
  S<-S+1
}

X<-X[,colSums(is.na(X))==0]
G.G<-tcrossprod(as.matrix(X))/ncol(X)

IDs<-factor(IDs,levels=rownames(G.G))
Z<-as.matrix(model.matrix(~IDs-1))
G.G<-tcrossprod(tcrossprod(Z,G.G),Z)  #This is ZGZ'

save(G.G,file='./G.G.rda')

EVD.G<-eigen(G.G)
rownames(EVD.G$vectors)<-rownames(G.G)

save(EVD.G,file='./EVD.G.rda')

############### Genotype by Environment (GxE) ############################
rm(list=ls())

load('./G.E.rda')
load('./G.G.rda')

G.GxE <- G.G*G.E # Hadamard product
EVD.GxE <- eigen(G.GxE)

save(G.GxE,file='./G.GxE.rda')
save(EVD.GxE,file='./EVD.GxE.rda')
############### Genotype by Weather covariates (GxW) ############################
rm(list=ls())

load('./G.G.rda')
load('./G.W.rda')

G.GxW <- G.G*G.W # Hadamard product
EVD.GxW <- eigen(G.GxW)

save(G.GxW,file='./G.GxW.rda')
save(EVD.GxW,file='./EVD.GxW.rda')


############### W matrix avg (AVG) ############################
rm(list=ls())

W<-read.csv("./WCONV.csv",sep=',',na.strings = "NA")
W<-W[,-1]
W<-as.matrix(W)
Y<-read.csv("./Y.csv",sep=";")
Y<-Y[Y$year!="2011",]
#Y<-Y[!Y$environ=="NE_2012",]

head(Y)
n<-dim(W)[1]
p<-dim(W)[2]

## Imputing, centering, standarizing. 
# Parameters for centering and standarizing

ctr <- TRUE
std <- TRUE

S<-0
for(i in 1:ncol(W)){
  meanWi <- mean(W[,i],na.rm=TRUE)
  # naive imputation
  W[,i] <- ifelse(is.na(W[,i]),meanWi, W[,i]) 
  if(ctr){ W[,i]<-W[,i]-meanWi } #centering
  if(std){ W[,i]<-W[,i]/sd(W[,i]) } # standarizing
  S<-S+1    # var(W[,i])
}

W<-W[,colSums(is.na(W))==0]
G.W<-tcrossprod(W)/ncol(W)


save(G.W,file='./G.W_conv.rda')

EVD.W_conv<-eigen(G.W)
rownames(EVD.W_conv$vectors)<-rownames(G.W)

save(EVD.W_conv,file='./EVD.W_conv.rda')
#image(G.W)

############### Genotype by weather cov AVG (GxW_AVG) ############################
rm(list=ls())

load('./G.G.rda')
load('./G.W_conv.rda')

G.GxW <- G.G*G.W # Hadamard product
EVD.GxW_conv <- eigen(G.GxW)

save(G.GxW,file='./G.GxW_conv.rda')
save(EVD.GxW_conv,file='./EVD.GxW_conv.rda')

############# Genotype by Environment ################
remove(list=ls())

Y<-read.csv("./Y.csv",sep=";")
Y<-Y[Y$year!="2011",]
Y$environ<-as.factor(Y$environ)
nEnv<-length(unique(Y$environ))
Ze<-model.matrix(~Y$environ-1)

load('./G.G.rda')

ZEZE<-tcrossprod(Ze)
GxE<-G.G*ZEZE

EVD.GxE<-eigen(GxE)
rownames(EVD.GxE$vectors)<-rownames(GxE)

save(GxE,file="./GxE.rda")
save(EVD.GxE,file="./EVD.GxE.rda")

####################### W matrix FILT ################
rm(list=ls())

W<-read.csv("./WFILT.csv",sep=',',na.strings = "NA")
W<-W[,-1]
W<-as.matrix(W)
Y<-read.csv("./Y.csv",sep=";")
#Y<-Y[!Y$environ=="NE_2012",]

head(Y)
n<-dim(W)[1]
p<-dim(W)[2]

## Imputing, centering, standarizing. 
# Parameters for centering and standarizing

ctr <- TRUE
std <- TRUE

S<-0
for(i in 1:ncol(W)){
  meanWi <- mean(W[,i],na.rm=TRUE)
  # naive imputation
  W[,i] <- ifelse(is.na(W[,i]),meanWi, W[,i]) 
  if(ctr){ W[,i]<-W[,i]-meanWi } #centering
  if(std){ W[,i]<-W[,i]/sd(W[,i]) } # standarizing
  S<-S+1    # var(W[,i])
}

W<-W[,colSums(is.na(W))==0]
#save(W,file="W.matrix.rda")
G.W<-tcrossprod(W)/ncol(W)


save(G.W,file='./G.W_filt_all.rda')

EVD.W_filt_all<-eigen(G.W)
rownames(EVD.W_filt_all$vectors)<-rownames(G.W)

save(EVD.W_filt_all,file='./EVD.W_filt_all.rda')
#image(G.W)

########### Genotype by weather FILT ##########
rm(list=ls())
setwd( "E:\\UF\\Work") 
load('./G.G.rda')
load('./G.W_filt_all.rda')

G.GxW <- G.G*G.W # Hadamard product
EVD.GxW_filt_all <- eigen(G.GxW)

save(G.GxW,file='./G.GxW_filt_all.rda')
save(EVD.GxW_filt_all,file='./EVD.GxW_filt_all.rda')

####################### W matrix STG ################
rm(list=ls())

W<-read.csv("./WSTG.csv",sep=',',na.strings = "NA")
W<-W[,-1]

W<-as.matrix(W)
W<-W[,colSums(is.na(W))==0]
#save(W,file="W.matrix.rda")
G.W<-tcrossprod(W)/ncol(W)


save(G.W,file='./G.W_Stages.rda')

EVD.W_Stages<-eigen(G.W)
rownames(EVD.W_Stages$vectors)<-rownames(G.W)

save(EVD.W_Stages,file='./EVD.W_Stages.rda')
#image(G.W)

########### Genotype by weather considering physiological stages ##########
rm(list=ls())

load('./G.G.rda') # Data2
load('./G.W_Stages.rda')

G.GxW <- G.G*G.W # Hadamard product
EVD.GxW_Stages <- eigen(G.GxW)

save(G.GxW,file='./G.GxW_Stages.rda')
save(EVD.GxW_Stages,file='./EVD.GxW_Stages.rda')
