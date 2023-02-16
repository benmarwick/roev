##Philip D. Gingerich RATES OF EVOLUTION (Cambridge University Press 2019) 
##8.2.04_Stearns1983_Gambusia
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
#============================================================ Load files
file1<-"C://R_aaROEVchapt08//8.2.04_Stearns1983_Gambusia.csv"
S<-read.csv(file1,header=TRUE,sep=",",quote="\"",dec=".",fill=TRUE)
S
#----------------------------------------- female age at maturity (days)
FA<-matrix(nrow=6,ncol=6)
  colnames(FA)=c('Kay','Twin','31','33','40','81')
  rownames(FA)=c('N','Mean','2SE','Stdev','Ln.M','Ln.S')
for (i in 1:6){
  FA[1,i]=S[1,3*i-1];FA[2,i]=S[1,3*i];FA[3,i]=S[1,3*i+1]
};FA=t(FA);FA
FA[,4]=.5*FA[,3]*sqrt(FA[,1]);FA[,5]=log(FA[,2]);FA[,6]=FA[,4]/FA[,2];FA
#----------------------------------------- female length at maturity (mm)
FL<-matrix(nrow=6,ncol=6)
  colnames(FL)=c('Kay','Twin','31','33','40','81')
  rownames(FL)=c('N','Mean','2SE','Stdev','Ln.M','Ln.S')
for (i in 1:6){
  FL[1,i]=S[2,3*i-1];FL[2,i]=S[2,3*i];FL[3,i]=S[2,3*i+1]
};FL=t(FL);FL
FL[,4]=.5*FL[,3]*sqrt(FL[,1]);FL[,5]=log(FL[,2]);FL[,6]=FL[,4]/FL[,2];FL
#----------------------------------------- female weight of offspring (mg)
FW<-matrix(nrow=6,ncol=6)
  colnames(FW)=c('Kay','Twin','31','33','40','81')
  rownames(FW)=c('N','Mean','2SE','Stdev','Ln.M','Ln.S')
for (i in 1:6){
  FW[1,i]=S[3,3*i-1];FW[2,i]=S[3,3*i];FW[3,i]=S[3,3*i+1]
};FW=t(FW);FW
FW[,4]=.5*FW[,3]*sqrt(FW[,1]);FW[,5]=log(FW[,2]);FW[,6]=FW[,4]/FW[,2];FW
#----------------------------------------- male age at maturity (days)
MA<-matrix(nrow=6,ncol=6)
  colnames(MA)=c('Kay','Twin','31','33','40','81')
  rownames(MA)=c('N','Mean','2SE','Stdev','Ln.M','Ln.S')
for (i in 1:6){
  MA[1,i]=S[4,3*i-1];MA[2,i]=S[4,3*i];MA[3,i]=S[4,3*i+1]
};MA=t(MA);MA
MA[,4]=.5*MA[,3]*sqrt(MA[,1]);MA[,5]=log(MA[,2]);MA[,6]=MA[,4]/MA[,2];MA
#----------------------------------------- male length at maturity (mm)
ML<-matrix(nrow=6,ncol=6)
  colnames(ML)=c('Kay','Twin','31','33','40','81')
  rownames(ML)=c('N','Mean','2SE','Stdev','Ln.M','Ln.S')
for (i in 1:6){
  ML[1,i]=S[5,3*i-1];ML[2,i]=S[5,3*i];ML[3,i]=S[5,3*i+1]
};ML=t(ML);ML
ML[,4]=.5*ML[,3]*sqrt(ML[,1]);ML[,5]=log(ML[,2]);ML[,6]=ML[,4]/ML[,2];ML
#=================================================Calculate rates
RFA<-matrix(nrow=15,ncol=8)
  colnames(RFA)=c('int','diff.sd','rate.sd.g','log.i','log|d|','log|r|',
    'S-B-N','wgt')
nrow=0
for (i in 1:5){
  for (j in (i+1):6){
    nrow=nrow+1
    RFA[nrow,1]=280
    poolsd=sqrt((FA[i,1]*FA[i,6]^2+FA[j,1]*FA[j,6]^2)/(FA[i,1]+FA[j,1]))
    RFA[nrow,2]=(FA[j,5]-FA[i,5])/poolsd
    RFA[nrow,3]=abs(RFA[nrow,2])/RFA[nrow,1]
    RFA[nrow,4]=log10(RFA[nrow,1]);RFA[nrow,5]=log10(abs(RFA[nrow,2]))
    RFA[nrow,6]=log10(RFA[nrow,3]);RFA[nrow,7]=2
    RFA[nrow,8]=1/RFA[nrow,1]
  }
};RFA
#----------------------------------------------------------RFL
RFL<-matrix(nrow=15,ncol=8)
  colnames(RFA)=c('int','diff.sd','rate.sd.g','log.i','log|d|','log|r|',
    'S-B-N','wgt')
nrow=0
for (i in 1:5){
  for (j in (i+1):6){
    nrow=nrow+1
    RFL[nrow,1]=280
    poolsd=sqrt((FL[i,1]*FL[i,6]^2+FL[j,1]*FL[j,6]^2)/(FL[i,1]+FL[j,1]))
    RFL[nrow,2]=(FL[j,5]-FL[i,5])/poolsd
    RFL[nrow,3]=abs(RFL[nrow,2])/RFL[nrow,1]
    RFL[nrow,4]=log10(RFL[nrow,1]);RFL[nrow,5]=log10(abs(RFL[nrow,2]))
    RFL[nrow,6]=log10(RFL[nrow,3]);RFL[nrow,7]=2
    RFL[nrow,8]=1/RFL[nrow,1]
  }
};RFL
#----------------------------------------------------------RFW
RFW<-matrix(nrow=15,ncol=8)
  colnames(RFA)=c('int','diff.sd','rate.sd.g','log.i','log|d|','log|r|',
    'S-B-N','wgt')
nrow=0
for (i in 1:5){
  for (j in (i+1):6){
    nrow=nrow+1
    RFW[nrow,1]=280
    poolsd=sqrt((FW[i,1]*FW[i,6]^2+FW[j,1]*FW[j,6]^2)/(FW[i,1]+FW[j,1]))
    RFW[nrow,2]=(FW[j,5]-FW[i,5])/poolsd
    RFW[nrow,3]=abs(RFW[nrow,2])/RFW[nrow,1]
    RFW[nrow,4]=log10(RFW[nrow,1]);RFW[nrow,5]=log10(abs(RFW[nrow,2]))
    RFW[nrow,6]=log10(RFW[nrow,3]);RFW[nrow,7]=2
    RFW[nrow,8]=1/RFW[nrow,1]
  }
};RFW
#----------------------------------------------------------RMA
RMA<-matrix(nrow=15,ncol=8)
  colnames(RFA)=c('int','diff.sd','rate.sd.g','log.i','log|d|','log|r|',
    'S-B-N','wgt')
nrow=0
for (i in 1:5){
  for (j in (i+1):6){
    nrow=nrow+1
    RMA[nrow,1]=280
    poolsd=sqrt((MA[i,1]*MA[i,6]^2+MA[j,1]*MA[j,6]^2)/(MA[i,1]+MA[j,1]))
    RMA[nrow,2]=(MA[j,5]-MA[i,5])/poolsd
    RMA[nrow,3]=abs(RMA[nrow,2])/RMA[nrow,1]
    RMA[nrow,4]=log10(RMA[nrow,1]);RMA[nrow,5]=log10(abs(RMA[nrow,2]))
    RMA[nrow,6]=log10(RMA[nrow,3]);RMA[nrow,7]=2
    RMA[nrow,8]=1/RMA[nrow,1]
  }
};RMA
#----------------------------------------------------------RML
RML<-matrix(nrow=15,ncol=8)
  colnames(RFA)=c('int','diff.sd','rate.sd.g','log.i','log|d|','log|r|',
    'S-B-N','wgt')
nrow=0
for (i in 1:5){
  for (j in (i+1):6){
    nrow=nrow+1
    RML[nrow,1]=280
    poolsd=sqrt((ML[i,1]*ML[i,6]^2+ML[j,1]*ML[j,6]^2)/(ML[i,1]+ML[j,1]))
    RML[nrow,2]=(ML[j,5]-ML[i,5])/poolsd
    RML[nrow,3]=abs(RML[nrow,2])/RML[nrow,1]
    RML[nrow,4]=log10(RML[nrow,1]);RML[nrow,5]=log10(abs(RML[nrow,2]))
    RML[nrow,6]=log10(RML[nrow,3]);RML[nrow,7]=2
    RML[nrow,8]=1/RML[nrow,1]
  }
};RML
#================================================Bind in one idr matrix
#idr=RFA;idr=rbind(idr,RFL);idr=rbind(idr,RFW)
#  idr=rbind(idr,RMA);idr=rbind(idr,RFL);idr
idr=rbind(RFA,RFL,RFW,RMA,RML)
idrx=idr[!rowSums(!is.finite(idr)),]  #remove rows that have -Inf
writefile<-"C://R_aaROEVchapt08//8.2.04_Stearns1983_Gambusia_out.csv"
write.csv(idrx,file=writefile,na="NA",row.names=FALSE) 

