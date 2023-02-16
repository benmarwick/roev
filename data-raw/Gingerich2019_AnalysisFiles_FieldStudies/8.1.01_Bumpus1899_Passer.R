##Philip D. Gingerich RATES OF EVOLUTION (Cambridge University Press 2019)
##8.1.01_Bumpus1899_Passer
#==========================================================================
cat(rep("\n",50))							#clear console
print (date())
rm(list=ls(all=TRUE))	#remove/clear all prev. variables
assign('last.warning',NULL,envir=baseenv())
ptm<-proc.time()
##======================================================================##
##library("devtools");library(roxygen2)
##setwd("c:/R_aaPackages/RATES");document();setwd("..");install("RATES")
#=====================================================================Setup
library(RATES);library(MASS)
#============================================================ Load file
  file1<-"../Gingerich2019_AnalysisFiles_FieldStudies/8.1.01_Bumpus1899_Passer.csv"
B<-read.csv(file1,header=TRUE,sep=",",quote="\"",dec=".",fill=TRUE)
B[1:3,]
#-----------------------------------------------------------------------
Bmean<-numeric(length=length(B[1,]));Bmean[1:length(B[1,])]=NA
Bsd<-numeric(length=length(B[1,]));Bsd[1:length(B[1,])]=NA
for (i in 3:length(B[1,])){
  Bmean[i]=mean(B[,i])
  Bsd[i]=sd(B[,i])
}
Bmean
Bsd
Bsd/Bmean		#coefficients of variation
BS<-B[B$SP=='S',];BS[1:3,]	#--------------------------------
BSmean<-numeric(length=length(BS[1,]));BSmean[1:length(B[1,])]=NA
BSsd<-numeric(length=length(BS[1,]));BSsd[1:length(B[1,])]=NA
BSN=length(BS[,1]);BSN
for (i in 3:length(BS[1,])){
  BSmean[i]=mean(BS[,i])
  BSsd[i]=sd(BS[,i])
}
BSmean
BSsd
BSsd/BSmean		#coefficients of variation
BSmean-Bmean
BSsd-Bsd

LB<-B;LB[,3:11]=log(B[,3:11]);LB[1:3,]	#--------------------------
LBN=length(LB[,1]);LBN
LBmean<-numeric(length=length(LB[1,]));LBmean[1:length(LB[1,])]=NA
LBsd<-numeric(length=length(LB[1,]));LBsd[1:length(LB[1,])]=NA
for (i in 3:length(LB[1,])){
  LBmean[i]=mean(LB[,i])
  LBsd[i]=sd(LB[,i])
}
LBS<-LB[LB$SP=='S',];LBS[1:3,]	#--------------------------------
LBSN=length(LBS[,1]);LBSN
LBSmean<-numeric(length=length(LBS[1,]));LBSmean[1:length(LBS[1,])]=NA
LBSsd<-numeric(length=length(LBS[1,]));LBSsd[1:length(LBS[1,])]=NA
for (i in 3:length(LBS[1,])){
  LBSmean[i]=mean(LBS[,i])
  LBSsd[i]=sd(LBS[,i])
}
c(LBmean,LBSmean)
c(LBsd,LBSsd)
#---------------------------------------------------------Rate calc.
LBSdiff<-numeric(length=length(LBS[1,]));LBSdiff[1:length(LBS[1,])]=NA
Psd<-numeric(length=length(LB[1,]));Psd[1:length(LB[1,])]=NA
Rate.sd<-numeric(length=length(LB[1,]));Rate.sd[1:length(LB[1,])]=NA
for (i in 1:length(LB[1,])){
  LBSdiff[i]=LBSmean[i]-LBmean[i]
  Psd[i]=sqrt(.5*(LBSsd[i]^2+LBsd[i]^2))
  Rate.sd[i]=LBSdiff[i]/Psd[i]
};LBSdiff;Psd;Rate.sd
#------------------------------------------------------------------
LBN
LBmean
LBsd
LBSN
LBSmean
LBSsd
LBSdiff
Psd
Rate.sd
log10(abs(Rate.sd))
median(abs(Rate.sd),na.rm=TRUE)
median(log10(abs(Rate.sd)),na.rm=TRUE)
10^median(log10(abs(Rate.sd)),na.rm=TRUE)
#-------------------------------------------------------Write archive file
idr<-matrix(nrow=length(LB[1,])-2,ncol=7)
  colnames(idr)=c('int','diff.sd','rate.sd','log.i','log.|d|','log.|r|','wgt')
idr[1:9,1]=c(1);idr[1:9,2]=Rate.sd[3:11];idr[1:9,3]=Rate.sd[3:11]
idr[1:9,4]=log10(idr[1:9,1]);idr[1:9,5]=log10(abs(idr[1:9,2]))
idr[1:9,6]=log10(abs(idr[1:9,3]));idr[1:9,7]=1/idr[1:9,1]
idr
idrx=idr[!rowSums(!is.finite(idr)),]  #remove rows that have -Inf
writefile<-"C://R_aaROEVchapt08//8.1.01_Bumpus1899_Passer_out.csv"
write.csv(idrx,file=writefile,na="NA",row.names=FALSE)

