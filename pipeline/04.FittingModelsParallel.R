###################################################
# Author: Vitor Sagae
# Fitting models
###################################################
rm(list=ls())

library(BGLR)
library(doSNOW)
library(foreach)
library(doParallel)

# Reads data sets
Y<-read.csv('./data/Y_CV.csv',sep=",",na.strings="NA")
Y<-Y[,-1]
load('./data/EVD.E.rda')
load('./data/EVD.G.rda')
load('./data/EVD.W.rda')
load('./data/EVD.W_conv.rda')
load('./data/EVD.W_filt_uc.rda')
load('./data/EVD.W_filt_all.rda')
load('./data/EVD.W_Stages.rda')
load('./data/EVD.GxW.rda')
load('./data/EVD.GxW_conv.rda')
load('./data/EVD.GxW_filt_uc.rda')
load('./data/EVD.GxW_filt_all.rda')
load('./data/EVD.GxW_Stages.rda')

nEnv<-length(unique(Y$environ))

ETA.new <- list()
ETA.new[[1]]<-list(V=EVD.E$vectors,d=EVD.E$values,model='RKHS')
ETA.new[[2]]<-list(V=EVD.G$vectors,d=EVD.G$values,model='RKHS')
ETA.new[[3]]<-list(V=EVD.W$vectors,d=EVD.W$values,model='RKHS')
ETA.new[[4]]<-list(V=EVD.GxW$vectors,d=EVD.GxW$values,model='RKHS')

ETA.conv <- list()
ETA.conv[[1]]<-list(V=EVD.E$vectors,d=EVD.E$values,model='RKHS')
ETA.conv[[2]]<-list(V=EVD.G$vectors,d=EVD.G$values,model='RKHS')
ETA.conv[[3]]<-list(V=EVD.W_conv$vectors,d=EVD.W_conv$values,model='RKHS')
ETA.conv[[4]]<-list(V=EVD.GxW_conv$vectors,d=EVD.GxW_conv$values,model='RKHS')


ETA.filt.all <- list()
ETA.filt.all[[1]]<-list(V=EVD.E$vectors,d=EVD.E$values,model='RKHS')
ETA.filt.all[[2]]<-list(V=EVD.G$vectors,d=EVD.G$values,model='RKHS')
ETA.filt.all[[3]]<-list(V=EVD.W_filt_all$vectors,d=EVD.W_filt_all$values,model='RKHS')
ETA.filt.all[[4]]<-list(V=EVD.GxW_filt_all$vectors,d=EVD.GxW_filt_all$values,model='RKHS')

ETA.stages <- list()
ETA.stages[[1]]<-list(V=EVD.E$vectors,d=EVD.E$values,model='RKHS')
ETA.stages[[2]]<-list(V=EVD.G$vectors,d=EVD.G$values,model='RKHS')
ETA.stages[[3]]<-list(V=EVD.W_Stages$vectors,d=EVD.W_Stages$values,model='RKHS')
ETA.stages[[4]]<-list(V=EVD.GxW_Stages$vectors,d=EVD.GxW_Stages$values,model='RKHS')

ETA.g <- list()
ETA.g[[1]]<-list(V=EVD.E$vectors,d=EVD.E$values,model='RKHS')
ETA.g[[2]]<-list(V=EVD.G$vectors,d=EVD.G$values,model='RKHS')

# Specifications
# Vector of phenotypes
colPhen <-14
envID<-as.factor(Y[,3])
IDs<-as.character(unique(Y$strain))
IDy<-Y$strain

env<-unique(envID)


# number of iterations
nIter <- 12000
# burnin period
burnIn <- 2000

#############################################################################
# Fitting models


# MODEL Y = E + G + W + GxW + e with the W matrix in a alternative format
# MODEL Y = E + G + W + GxW + e with the W matrix in a conventional format by the mean
# MODEL Y = E + G + W + GxW + e with the W matrix filtered and using usefulnes criterion
# MODEL Y = E + G + W + GxW + e with the W matrix filtered alternative format

## Assigns lines to folds (CV1) ###########################

length(IDs)

#Number of folds
set.seed(123)
folds = 5
partitions= 10
fold <- rep(1:folds,each=ceiling(length(IDs)/folds))[order(runif(length(IDs)))]

namescol<-expand.grid(1:folds,1:partitions)
nameslist<-list()
for (g in 1:nrow(namescol)){
  nameslist[[g]]<-paste("fold.",namescol[g,1],".rep.",namescol[g,2],sep="")
}
namescol<-do.call(c,nameslist)
#Aleatoriza dentro de cada environment e posteriormente repete a mesma coisa nos outros (mesmo id em todos ambientes)

#CV1
Y$CV1 <- NA
for(i in 1:folds){
  tmp<-Y[,4]%in%IDs[fold==i]
  Y$CV1[tmp]<-i
}

# Fitting CV1
library(foreach)
library(doSNOW)

my.cluster <- makeCluster(6)
doParallel::registerDoParallel(cl = my.cluster)
registerDoSNOW(my.cluster)


# Fitting CV1
COR<-matrix(0,nrow=nEnv,ncol=(folds*1))
COR.conv<-matrix(0,nrow=nEnv,ncol=(folds*1))
COR.filt.all<-matrix(0,nrow=nEnv,ncol=(folds*partitions),dimnames = list(1:nEnv,namescol))
COR.stages<-matrix(0,nrow=nEnv,ncol=(folds*1))
COR.g<-matrix(0,nrow=nEnv,ncol=(folds*1))

corCV1<-foreach(s=1:partitions,.combine="cbind",.packages=c("dplyr","BGLR"),.multicombine=T)%dopar%{
  for (i in 1:folds){
    NA.y<-Y
    NA.y[NA.y$CV1==i,colPhen]<-NA
    A<-as.matrix(NA.y[,colPhen])
    fm<-BGLR(y=A,ETA=ETA.new,nIter=nIter,burnIn=burnIn)
    fm.conv<-BGLR(y=A,ETA=ETA.conv,nIter=nIter,burnIn=burnIn)
    fm.filt.all<-BGLR(y=A,ETA=ETA.filt.all,nIter=nIter,burnIn=burnIn)
    fm.stages<-BGLR(y=A,ETA=ETA.stages,nIter=nIter,burnIn=burnIn)
    fm.g<-BGLR(y=A,ETA=ETA.g,nIter=nIter,burnIn=burnIn)
    
    for (g in 1:length(env)){
      tst1<-which(is.na(NA.y[Y$environ==env[g],colPhen]))
      COR[g,((1-1)*folds+i)]<-cor(Y[Y$environ==env[g],colPhen][tst1],fm$yHat[Y$environ==env[g]][tst1])
      COR.conv[g,((1-1)*folds+i)]<-cor(Y[Y$environ==env[g],colPhen][tst1],fm.conv$yHat[Y$environ==env[g]][tst1])
      COR.filt.all[g,((s-1)*folds+i)]<-cor(Y[Y$environ==env[g],colPhen][tst1],fm.filt.all$yHat[Y$environ==env[g]][tst1])
      COR.stages[g,((1-1)*folds+i)]<-cor(Y[Y$environ==env[g],colPhen][tst1],fm.stages$yHat[Y$environ==env[g]][tst1])
      COR.g[g,((1-1)*folds+i)]<-cor(Y[Y$environ==env[g],colPhen][tst1],fm.g$yHat[Y$environ==env[g]][tst1])
      
    }
  }
  return(rbind(cbind("COR",COR),cbind("cor.mean",COR.conv),cbind("cor.stg",COR.stages),cbind("cor.g",COR.g)))
}





dir.create(paste(".\\Output_CV1\\",sep=""))
filename<-paste(".\\Output_CV1\\",sep="")
write.csv(corCV1,file=paste(file.path(filename),"CORCV1.csv",sep=""),row.names=F)

# ----------------------------------------------------------
# Fitting CV2
COR<-matrix(0,nrow=nEnv,ncol=(folds*1))
COR.conv<-matrix(0,nrow=nEnv,ncol=(folds*1))
COR.filt.all<-matrix(0,nrow=nEnv,ncol=(folds*partitions),dimnames = list(1:nEnv,namescol))
COR.stages<-matrix(0,nrow=nEnv,ncol=(folds*1))
COR.g<-matrix(0,nrow=nEnv,ncol=(folds*1))

corCV2<-foreach(s=1:partitions,.combine="cbind",.packages=c("dplyr","BGLR"),.multicombine=T)%dopar%{
  
  for (i in 1:folds){
    NA.y<-Y
    NA.y[NA.y$CV2==i,colPhen]<-NA
    A<-as.matrix(NA.y[,colPhen])
    fm<-BGLR(y=A,ETA=ETA.new,nIter=nIter,burnIn=burnIn)
    fm.conv<-BGLR(y=A,ETA=ETA.conv,nIter=nIter,burnIn=burnIn)
    fm.filt.all<-BGLR(y=A,ETA=ETA.filt.all,nIter=nIter,burnIn=burnIn)
    fm.stages<-BGLR(y=A,ETA=ETA.stages,nIter=nIter,burnIn=burnIn)
    fm.g<-BGLR(y=A,ETA=ETA.g,nIter=nIter,burnIn=burnIn)
    
    for (g in 1:length(env)){
      tst1<-which(is.na(NA.y[Y$environ==env[g],colPhen]))
      COR[g,((s-1)*folds+i)]<-cor(Y[Y$environ==env[g],colPhen][tst1],fm$yHat[Y$environ==env[g]][tst1])
      COR.conv[g,((s-1)*folds+i)]<-cor(Y[Y$environ==env[g],colPhen][tst1],fm.conv$yHat[Y$environ==env[g]][tst1])
      COR.filt.all[g,((s-1)*folds+i)]<-cor(Y[Y$environ==env[g],colPhen][tst1],fm.filt.all$yHat[Y$environ==env[g]][tst1])
      COR.stages[g,((s-1)*folds+i)]<-cor(Y[Y$environ==env[g],colPhen][tst1],fm.stages$yHat[Y$environ==env[g]][tst1])
      COR.g[g,((s-1)*folds+i)]<-cor(Y[Y$environ==env[g],colPhen][tst1],fm.g$yHat[Y$environ==env[g]][tst1])
      
    }
  }
  return(rbind(cbind("COR",COR),cbind("cor.mean",COR.conv),cbind("cor.stg",COR.stages),cbind("cor.g",COR.g)))
}

dir.create(paste(".\\Output_CV2\\",sep=""))
filename<-paste(".\\Output_CV2\\",sep="")
write.csv(corCV2,file=paste(file.path(filename),"CORCV2.csv",sep=""),row.names=F)

# -------------------------------------------------------------------------------------------

# Fitting CV0
COR<-matrix(0,nrow=nEnv,ncol=(folds*1))
COR.conv<-matrix(0,nrow=nEnv,ncol=(folds*1))
COR.filt.all<-matrix(0,nrow=nEnv,ncol=(folds*partitions),dimnames = list(1:nEnv,namescol))
COR.stages<-matrix(0,nrow=nEnv,ncol=(folds*1))
COR.g<-matrix(0,nrow=nEnv,ncol=(folds*1))

corCV0<-foreach(s=1:partitions,.combine="cbind",.packages=c("dplyr","BGLR"),.multicombine=T)%dopar%{
  
  for (i in 1:length(env)){
    NA.y<-Y
    NA.y[NA.y$CV0==i,colPhen]<-NA
    A<-as.matrix(NA.y[,colPhen])
    fm<-BGLR(y=A,ETA=ETA.new,nIter=nIter,burnIn=burnIn)
    fm.conv<-BGLR(y=A,ETA=ETA.conv,nIter=nIter,burnIn=burnIn)
    fm.filt.all<-BGLR(y=A,ETA=ETA.filt.all,nIter=nIter,burnIn=burnIn)
    fm.stages<-BGLR(y=A,ETA=ETA.stages,nIter=nIter,burnIn=burnIn)
    fm.g<-BGLR(y=A,ETA=ETA.g,nIter=nIter,burnIn=burnIn)
    
    #  
    tst1<-fm$whichNa
    COR[i,partitions]<-cor(Y[tst1,colPhen],fm$yHat[tst1])
    COR.conv[i,partitions]<-cor(Y[tst1,colPhen],fm.conv$yHat[tst1])
    COR.filt.all[i,partitions]<-cor(Y[tst1,colPhen],fm.filt.all$yHat[tst1])
    COR.stages[i,partitions]<-cor(Y[tst1,colPhen],fm.stages$yHat[tst1])
    COR.g[i,partitions]<-cor(Y[tst1,colPhen],fm.g$yHat[tst1])
    
}
  return(rbind(cbind("COR",COR),cbind("cor.mean",COR.conv),cbind("cor.stg",COR.stages),cbind("cor.stg",COR.stages)))
}

dir.create(paste(".\\Output_CV0\\",sep=""))
filename<-paste(".\\Output_CV0\\",sep="")
write.csv(corCV0,file=paste(file.path(filename),"CORCV0.csv",sep=""),row.names=F)

close(pb)
parallel::stopCluster(cl = my.cluster)
