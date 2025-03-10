###################################################
# Author: Vitor Sagae
# Formating environmental data
###################################################

remove(list=ls())

library(tidyverse)

#Carregando dados fenotipicos
dados<-read.csv("./data/Y.csv",sep=";",na.strings = "NA")
dados<-dados[dados$year!="2011",] # Selecting just 2012
#dados<-dados[!dados$environ=="NE_2012",] # removing a environment

# Transformando colunas para fatores
dados$location<-as.factor(dados$location)
dados$strain<-as.factor(dados$strain)
dados$year<-as.factor(dados$year)
dados$environ<-as.factor(dados$environ)
dados<-dados |> droplevels()

m_Env <- tapply(dados$yield, dados$environ,mean)
m_Env
#save(m_Env,file="./data/m_Env.rda")
# -------------------------------------------

# WP
We_data<-read.csv("./data/env.cov.csv",sep=",")
We_data<-We_data[,-1]
#We_data<-We_data[!We_data$ENV=="NE_2012",] # Removing a environment

We_data2<-reshape(data=We_data,
                  idvar="ENV",
                  v.names=c("T2M_MAX","T2M_MIN","RH2M","T2M","PRECTOTCORR","WS2M","QV2M"),
                  timevar="DOY",
                  direction="wide")
We_data2<-We_data2[,-c(3,4,5)]
colnames(We_data2)[2]<-"environ"

dados.to.join<-dados[,c(1:4)]
We.data.joined<-full_join(dados.to.join,We_data2[,-1],by="environ")

W<-We.data.joined[,-c(1:5)]

WP<-W
write.csv(W,"./data/W.csv")
write.csv(WP,"./data/WP.csv")

# -------------------------------------------
# AVG
We_data<-read.csv("./data/env.cov.csv",sep=",")
#We_data<-We_data[!We_data$ENV=="NE_2012",] # Removing a environment
We_conv<-We_data |> group_by(ENV) |> summarise_at(vars(T2M_MAX,T2M_MIN,RH2M,T2M,PRECTOTCORR,WS2M,QV2M),mean)

colnames(We_conv)[1]<-"environ"
dados.to.join<-dados[,c(1:4)]
We.conv.joined<-full_join(dados.to.join,We_conv,by="environ")

W<-We.conv.joined[,-c(1:4)]
WCONV<-W
write.csv(WCONV,"./data/WCONV.csv")


# -------------------------------------------
# Filt
We_data<-read.csv("./data/env.cov.csv",sep=",")
We_data<-We_data[,-1]
We_data$DOY<-We_data$DOY-130

filtro<-c(14:22,41:51,56:65,96:111,142:171,100:111,139:147,35:51,96:110,101:132,143:175,20:35,37:57,78:94,126:143,152:189)
We_data<-We_data |> filter(DOY%in%filtro)

We_data2<-reshape(data=We_data,
                  idvar="ENV",
                  v.names=c("T2M_MAX","T2M_MIN","RH2M","T2M","PRECTOTCORR","WS2M","QV2M"),
                  timevar="DOY",
                  direction="wide")
We_data2<-We_data2[,-c(3,4,5)]
colnames(We_data2)[2]<-"environ"

dados.to.join<-dados[,c(1:4)]
We.data.joined<-full_join(dados.to.join,We_data2[,-1],by="environ")

W<-We.data.joined[,-c(1:5)]
WFILT<-W

write.csv(WFILT,"./data/WFILT.csv")

# --------------------------------------------------------------------

# STG
We_data<-read.csv("./data/W_Stages.csv",sep=",")
We_data<-We_data[,-1]

We_data<-cbind(We_data2[,2],We_data)
colnames(We_data)[1]<-"environ"

dados.to.join<-dados[,c(1:4)]
We.data.joined<-full_join(dados.to.join,We_data,by="environ")

W<-We.data.joined[,-c(1:4)]
WSTG<-W
write.csv(WSTG,"./data/WSTG.csv")

        