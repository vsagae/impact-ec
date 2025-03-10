###################################################
# Author: Vitor Sagae
# Extracting environmental data using nasapower
###################################################

rm(list = ls())

library(tidyverse)
# Location
LongLat1 <- c(-93.734651,42.015252)  #IOWA 2012
LongLat3 <- c(-88.221499,40.054879) #ILLINOIS 2012
LongLat4 <- c(-86.992231,40.471148) #INDIANA 2012
LongLat6<- c(-98.137824,40.575256) #NEBRASKA 2012

lonlat <- cbind(LongLat1,LongLat3,LongLat4,LongLat6 )

Y<-read.csv("./data/Y.csv",sep=";")
Y<-Y[Y$year!="2011",]
locnames<-unique(Y$environ)

# Year
date_Inteval_2 <- c("2012-05-10", "2012-11-14")

# Removing date_interval1
datas <- cbind("2012"=date_Inteval_2)
a = seq(from = as.Date(min(datas[,1])), to = as.Date(max(datas[,1])), by = 'day')
nDays <- length(a)

#Env data
#install.packages("nasapower")
library(nasapower)

NEnvVar <- 8
NEnv <- 4
#Matrix to receive the data


Env <- 4
Year <- 1
YearNames<-c("2012")

res<-list()
Env_M_data<-list()

for(i in 1:Year){
  for(j in 1:Env){
    daily_single_data <- get_power(
      community = "ag",
      lonlat = t(lonlat[,j]),
      pars = c("T2M_MAX", "T2M_MIN","RH2M", "T2M", "PRECTOTCORR", "WS2M","QV2M"),
      dates = t(datas[,i]),
      temporal_api = "daily")
   daily_single_data<-data.frame(ENV=paste(locnames[j],sep=""),daily_single_data)
   
   res[[j]]<-daily_single_data
   
  }
  Env_M_data[[i]]<-do.call(rbind,res)
} 

Env_M_data<-do.call(rbind,Env_M_data)
Env_M_data<-Env_M_data[,-c(2,3)]
Env_M_data<-data.frame(YearFilt=sub("^...", "", Env_M_data[,1]),Env_M_data)
Env_M_data_filt<-subset(Env_M_data,(as.character(YearFilt)==as.character(YEAR)))

write.csv(Env_M_data_filt,"./data/env.cov.csv")


# Splitting into stages
stg<-data.frame(ENV=c("IA_2012","IL_2012","IN_2012","NE_2012"),
                plant.date = as.Date(c("2012-05-21","2012-05-10","2012-05-16","2012-05-10")),
                harv.date = as.Date(c("2012-09-20","2012-09-20","2012-09-22","2012-09-18")))

time.windows=c(0,20,35,45,60,69,78,93,111,132) #Planting->V2->V5->R1->R2->R3->R4->R5->R6->R7

sel.Env_M_data<-Env_M_data_filt %>% left_join(stg,by="ENV") %>% filter(YYYYMMDD>=plant.date & cicle.length<=harv.date)
  

sel.Env_M_data <- sel.Env_M_data %>%
  group_by(ENV) %>%
  arrange(YYYYMMDD) %>% # Garante que as datas estão ordenadas
  mutate(day = as.numeric(YYYYMMDD - first(YYYYMMDD))) %>% 
  ungroup()

sel.Env_M_data <- sel.Env_M_data %>%
  mutate(across(c("T2M_MAX", "T2M_MIN","RH2M", "T2M", "PRECTOTCORR", "WS2M","QV2M"), ~as.numeric(scale(.x)))) %>% 
  ungroup()

sel.Env_M_data <- sel.Env_M_data %>%
  mutate(windows = cut(day, breaks = time.windows, include.lowest = FALSE, right = FALSE,labels=FALSE))


w<-sel.Env_M_data %>% group_by(ENV,windows) %>% summarise(across(c("T2M_MAX", "T2M_MIN","RH2M", "T2M", "PRECTOTCORR", "WS2M","QV2M"),mean,na.rm=T))

w$windows<-as.integer(rep(time.windows,times=4))
w$ENV<-as.character(w$ENV)

w<-as.data.frame(w)
w_stg<-reshape(data=w,
                  idvar="ENV",
                  v.names=c("T2M_MAX","T2M_MIN","RH2M","T2M","PRECTOTCORR","WS2M","QV2M"),
                  timevar="windows",
                  direction="wide")

write.csv(w_stg,"./data/W_stages.csv")
