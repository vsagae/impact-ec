###################################################
# Author: Vitor Sagae
# Assign to cross-validation
###################################################
rm(list=ls())

# Reads data sets
Y<-read.csv('./data/Y.csv',sep=";",na.strings="NA")
Y<-Y[Y$year!="2011",]
# Specifications
# Vector of phenotypes
colPhen <-14
envID<-as.factor(Y[,3])
IDs<-as.character(unique(Y$strain))
IDy<-Y$strain

env<-unique(envID)

###Preparing CV1 e CV2 and CV0

#Getting the genotypes ID
colIDy<-4
IDs <- unique(Y[,colIDy])

## Assigns lines to folds (CV1) ###########################

length(IDs)

#Number of folds
set.seed(123)
folds = 5
partitions= 1
fold <- rep(1:folds,each=ceiling(length(IDs)/folds))[order(runif(length(IDs)))]

#CV1
Y$CV1 <- NA
for(i in 1:folds){
  tmp<-Y[,4]%in%IDs[fold==i]
  Y$CV1[tmp]<-i
}


#CV2
Y$CV2 <- NA
for(i in IDs){
  tmp<-which(Y[,4]==i)
  ni<-length(tmp)
  fold<-sample(1:folds,size=ni,replace=length(tmp)>folds)
  Y$CV2[tmp]<-fold
  #print(table(fold))
}


#CV0
#Getting the Environment ID
IDsEnv <- unique(Y[,3])

## Assigns Environment to folds (CV0) ###########################
length(IDsEnv)

#Number of folds
folds = length(IDsEnv)
fold <- sample(1:length(IDsEnv),length(IDsEnv))

#CV0
Y$CV0 <- NA
for(i in 1:folds){
  tmp<-Y[,3]%in%IDsEnv[fold==i]
  Y$CV0[tmp]<-i
}

write.csv(Y,file="./data/Y_CV.csv")
