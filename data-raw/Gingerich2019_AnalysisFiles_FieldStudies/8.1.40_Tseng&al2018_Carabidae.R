##Philip D. Gingerich RATES OF EVOLUTION (Cambridge University Press 2019) 
##8.1.40_Tseng&al2018_Carabidae
#==========================================================================
cat(rep("\n",50))							#clear console
print (date()) #Ctr-s to save, Ctr-a to select all, Ctr-r to run
rm(list=ls(all=TRUE))#remove/clear all prev. variables
assign('last.warning',NULL,envir=baseenv())
ptm<-proc.time()
##======================================================================##
##library("devtools");library(roxygen2)
##setwd("c:/R_aaPackages/RATES");document();setwd("..");install("RATES") 
#=====================================================================Setup
library(RATES);library(MASS)
#--------------------------------------------------------------------------
list.files("C://R_aaRateBook//RB_Chapter08",pattern=".csv")
file1<-"C://R_aaROEVchapt08//8.1.40_Tseng&al2018_Carabidae.csv" 
T<-read.csv(file1,header=TRUE,sep=",",quote="\"",dec=".",fill=TRUE)
T[1:4,]
nrow(T)
minn=10;minyr=min(T[,3]);maxyr=max(T[,3])
#------------------------------------------------------Logged matrix
LT<-T[,c(2,3,5,7)]
  colnames(T)=c('spec','year','region','ln_ely_mm');LT[1:4,]
for (i in 1:nrow(T)){
  LT[i,4]=log(T[i,7])
};LT[1:4,]
#------------------------------------------------------ Pterostichus algidus
algidus<-matrix(nrow=0,ncol=4);colnames(algidus)=c('yr','N','ln_ely','s_ely');algidus
for (i in minyr:maxyr){
  temp<-matrix(nrow=0,ncol=4);colnames(temp)=c('spec','year','region','ln_ely_mm')
  nc=0
  for(j in 1:nrow(LT)){
    if(LT[j,1]=='algidus'&&LT[j,2]==i){
      nc=nc+1
      temp=rbind(temp,LT[j,])
    }
  }
  if(nc>minn){
    al1=i;al2=nrow(temp);al3=mean(temp[,4]);al4=sd(temp[,4])
    algidus=rbind(algidus,c(al1,al2,al3,al4))
  }
}
algidus
#------------------------------------------------------ Scaphinotus angusticollis
angust<-matrix(nrow=0,ncol=4);colnames(angust)=c('yr','N','ln_ely','s_ely');angust
for (i in minyr:maxyr){
  temp<-matrix(nrow=0,ncol=4);colnames(temp)=c('spec','year','region','ln_ely_mm')
  nc=0
  for(j in 1:nrow(LT)){
    if(LT[j,1]=='angusticollis'&&LT[j,2]==i){
      nc=nc+1
      temp=rbind(temp,LT[j,])
    }
  }
  if(nc>minn){
    al1=i;al2=nrow(temp);al3=mean(temp[,4]);al4=sd(temp[,4])
    angust=rbind(angust,c(al1,al2,al3,al4))
  }
}
angust
#------------------------------------------------------ Pterostichus melanarius
melana<-matrix(nrow=0,ncol=4);colnames(melana)=c('yr','N','ln_ely','s_ely');melana
for (i in minyr:maxyr){
  temp<-matrix(nrow=0,ncol=4);colnames(temp)=c('spec','year','region','ln_ely_mm')
  nc=0
  for(j in 1:nrow(LT)){
    if(LT[j,1]=='melanarius'&&LT[j,2]==i){
      nc=nc+1
      temp=rbind(temp,LT[j,])
    }
  }
  if(nc>minn){
    al1=i;al2=nrow(temp);al3=mean(temp[,4]);al4=sd(temp[,4])
    melana=rbind(melana,c(al1,al2,al3,al4))
  }
}
melana
#------------------------------------------------------ Carabus nemoralis
nemora<-matrix(nrow=0,ncol=4);colnames(nemora)=c('yr','N','ln_ely','s_ely');nemora
for (i in minyr:maxyr){
  temp<-matrix(nrow=0,ncol=4);colnames(temp)=c('spec','year','region','ln_ely_mm')
  nc=0
  for(j in 1:nrow(LT)){
    if(LT[j,1]=='nemoralis'&&LT[j,2]==i){
      nc=nc+1
      temp=rbind(temp,LT[j,])
    }
  }
  if(nc>minn){
    al1=i;al2=nrow(temp);al3=mean(temp[,4]);al4=sd(temp[,4])
    nemora=rbind(nemora,c(al1,al2,al3,al4))
  }
}
nemora
#------------------------------------------------------ Harpalus fraternus
frater<-matrix(nrow=0,ncol=4);colnames(frater)=c('yr','N','ln_ely','s_ely');frater
for (i in minyr:maxyr){
  temp<-matrix(nrow=0,ncol=4);colnames(temp)=c('spec','year','region','ln_ely_mm')
  nc=0
  for(j in 1:nrow(LT)){
    if(LT[j,1]=='fraternus'&&LT[j,2]==i){
      nc=nc+1
      temp=rbind(temp,LT[j,])
    }
  }
  if(nc>minn){
    al1=i;al2=nrow(temp);al3=mean(temp[,4]);al4=sd(temp[,4])
    frater=rbind(frater,c(al1,al2,al3,al4))
  }
}
frater
#------------------------------------------------------ Euryderus grossus
grossus<-matrix(nrow=0,ncol=4);colnames(grossus)=c('yr','N','ln_ely','s_ely');grossus
for (i in minyr:maxyr){
  temp<-matrix(nrow=0,ncol=4);colnames(temp)=c('spec','year','region','ln_ely_mm')
  nc=0
  for(j in 1:nrow(LT)){
    if(LT[j,1]=='grossus'&&LT[j,2]==i){
      nc=nc+1
      temp=rbind(temp,LT[j,])
    }
  }
  if(nc>minn){
    al1=i;al2=nrow(temp);al3=mean(temp[,4]);al4=sd(temp[,4])
    grossus=rbind(grossus,c(al1,al2,al3,al4))
  }
}
grossus
#------------------------------------------------------ Cymindis planipennis
planip<-matrix(nrow=0,ncol=4);colnames(planip)=c('yr','N','ln_ely','s_ely');planip
for (i in minyr:maxyr){
  temp<-matrix(nrow=0,ncol=4);colnames(temp)=c('spec','year','region','ln_ely_mm')
  nc=0
  for(j in 1:nrow(LT)){
    if(LT[j,1]=='planipennis'&&LT[j,2]==i){
      nc=nc+1
      temp=rbind(temp,LT[j,])
    }
  }
  if(nc>minn){
    al1=i;al2=nrow(temp);al3=mean(temp[,4]);al4=sd(temp[,4])
    planip=rbind(planip,c(al1,al2,al3,al4))
  }
}
planip
#------------------------------------------------------ Amara quenseli
quenseli<-matrix(nrow=0,ncol=4);colnames(quenseli)=c('yr','N','ln_ely','s_ely');quenseli
for (i in minyr:maxyr){
  temp<-matrix(nrow=0,ncol=4);colnames(temp)=c('spec','year','region','ln_ely_mm')
  nc=0
  for(j in 1:nrow(LT)){
    if(LT[j,1]=='quenseli'&&LT[j,2]==i){
      nc=nc+1
      temp=rbind(temp,LT[j,])
    }
  }
  if(nc>minn){
    al1=i;al2=nrow(temp);al3=mean(temp[,4]);al4=sd(temp[,4])
    quenseli=rbind(quenseli,c(al1,al2,al3,al4))
  }
}
quenseli
#================================================== Calculate rates 
gt=1
idr=matrix(nrow=0,ncol=9)
  colnames(idr)=c('int','diff.sd','rate.sd.g','log.i','log|d|','log|r|',
    'sbn','wgt','sp.no.')
for(sn in 1:8){								#species number
  if(sn==1){test=algidus}
  if(sn==2){test=angust}
  if(sn==3){test=melana}
  if(sn==4){test=nemora}
  if(sn==5){test=frater}
  if(sn==6){test=grossus}
  if(sn==7){test=planip}
  if(sn==8){test=quenseli}
  for (sr in 1:(nrow(test)-1)){ 					#start row
    for (i in 1:(nrow(test)-sr)){					#interval steps
      id<-numeric(length=9)
      k=test[sr+i,1]-test[sr,1]					#interval in yrs
      if(k<=gt){id[1]=1}else{id[1]=k/gt}				#interval in gen
      meandiff=test[(sr+i),3]-test[sr,3]				#mean diff.
      psd=PoolSD(test[sr,2],test[sr+i,2],test[sr,4],test[sr+i,4])#n1,n2,sd1,sd2
      id[2]=meandiff/psd							#diff.sd
      id[3]=abs(id[2])/id[1]						#rate.sd.gen
      id[4]=log10(id[1])							#log.i
      id[5]=log10(abs(id[2]))						#log.d
      id[6]=log10(id[3])							#log.r
      if(i==1){
        if(id[1]<=1){id[7]=1}else{id[7]=2}				#sbn
      }else{id[7]=3}
      id[8]=1/id[1]							#wgt
      id[9]=sn #id[4]+id[6]						#sum
      idr=rbind(idr,id)
    }
  }
}
nrow(idr)
plot(idr[,4],idr[,6])
min(idr[,9]);max(idr[,9])

#--------------------------------------------------------------------------
idrx=idr[!rowSums(!is.finite(idr)),]  #remove rows that have -Inf
idrxC=idrx[,1:8]
idrc=idrx[idrx[,9]>=1,]  #keep only rows where interval >= gentime
idrcC=idrc[,1:8]
nrow(idr);nrow(idrx);nrow(idrxC);nrow(idrcC)
writefile<-"C://R_aaROEVchapt08//8.1.40_Tseng&al2018_Carabidae_out.csv"
write.csv(idrx,file=writefile,na="NA",row.names=FALSE) 


