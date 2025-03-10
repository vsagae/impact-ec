###################################################
# Author: Vitor Sagae
###################################################

rm(list=ls())

library(BGLR)
library(doSNOW)
library(foreach)
library(doParallel)
library(foreach)
library(doSNOW)
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

ETA.filt <- list()
ETA.filt[[1]]<-list(V=EVD.E$vectors,d=EVD.E$values,model='RKHS')
ETA.filt[[2]]<-list(V=EVD.G$vectors,d=EVD.G$values,model='RKHS')
ETA.filt[[3]]<-list(V=EVD.W_filt_uc$vectors,d=EVD.W_filt_uc$values,model='RKHS')
ETA.filt[[4]]<-list(V=EVD.GxW_filt_uc$vectors,d=EVD.GxW_filt_uc$values,model='RKHS')

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

# FULL MODELS
NA.y<-Y
A<-as.matrix(NA.y[,colPhen])
fm<-BGLR(y=A,ETA=ETA.new,nIter=nIter,burnIn=burnIn)
fm.conv<-BGLR(y=A,ETA=ETA.conv,nIter=nIter,burnIn=burnIn)
fm.filt.uc<-BGLR(y=A,ETA=ETA.filt,nIter=nIter,burnIn=burnIn)
fm.stages<-BGLR(y=A,ETA=ETA.stages,nIter=nIter,burnIn=burnIn)
fm.g<-BGLR(y=A,ETA=ETA.g,nIter=nIter,burnIn=burnIn)

save(fm,file="fm.rda")
save(fm.conv,file="fm.avg.rda")
save(fm.filt.uc,file="fm.filt.rda")
save(fm.stages,file="fm.stg.rda")
save(fm.g,file="fm.g.rda")