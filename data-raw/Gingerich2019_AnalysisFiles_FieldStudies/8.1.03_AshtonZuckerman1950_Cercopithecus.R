##Philip D. Gingerich RATES OF EVOLUTION (Cambridge University Press 2019) 
##8.1.03_AshtonZuckerman1950_Cercopithecus
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
#============================================================ Load file
file1<-"C://R_aaROEVchapt08//8.1.03_AshtonZuckerman1950_Cercopithecus.csv"
A<-read.csv(file1,header=TRUE,sep=",",quote="\"",dec=".",fill=TRUE)
A[1:3,]
#-----------------------------------------------------------------------
LA<-A[,1:8]
  colnames(LA)=c('Tooth','Meas.','LGmean','LGsd','Gn','LKmean','LKsd','Kn')
LA[,3]=log(A[,3]);LA[,6]=log(A[,6])		#transform means
LA[,4]=A[,4]*sqrt(A[,5])/A[,3];LA[,7]=A[,7]*sqrt(A[,8])/A[,6]
LA[1:3,]
idr<-matrix(nrow=length(A[,1]),ncol=7)
  colnames(idr)=c('int','diff.sd','rate.sd','log.i','log.|d|','log.|r|','wgt')
psd<-numeric(length=length(A[,1]))		#pooled standard deviation
gen=30						#interval in generations
for (i in 1:length(A[,1])){
  idr[i,1]=gen
  psd[i]=sqrt(.5*(LA[i,4]^2+LA[i,7]^2))
  idr[i,2]=(LA[i,6]-LA[i,3])/psd[i]
  idr[i,3]=abs(idr[i,2])/idr[i,1]
  idr[i,4]=log10(abs(idr[i,1]));idr[i,5]=log10(abs(idr[i,2]))
  idr[i,6]=log10(idr[i,3]);idr[i,7]=1/idr[i,1]
}
idr[1:3,]
idrx=idr[!rowSums(!is.finite(idr)),]  #remove rows that have -Inf
writefile<-"C://R_aaROEVchapt08//8.1.03_AshtonZuckerman1950_Cercopithecus_out.csv"
write.csv(idrx,file=writefile,na="NA",row.names=FALSE) 
