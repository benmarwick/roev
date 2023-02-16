##Philip D. Gingerich RATES OF EVOLUTION (Cambridge University Press 2019) 
##8.2.01_JohnstonSelander1971_Passer
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
filename<-vector(length=6)
filename[1]="8.2.01_JohnstonSelander1971_Passer_fig04m"
filename[2]="8.2.01_JohnstonSelander1971_Passer_fig04f"
filename[3]="8.2.01_JohnstonSelander1971_Passer_fig05m"
filename[4]="8.2.01_JohnstonSelander1971_Passer_fig05f"
filename[5]="8.2.01_JohnstonSelander1971_Passer_fig06m"
filename[6]="8.2.01_JohnstonSelander1971_Passer_fig06f"
for (run in 1:6){
  file1<-paste("C://R_aaROEVchapt08//",filename[run],".csv",sep="")
  X<-read.csv(file1,header=TRUE,sep=",",quote="\"",dec=".",fill=TRUE);X
  Xjs=X[X$Keep=="Keep",];Xjs
  JS<-matrix(nrow=length(Xjs[,1]),ncol=3)
    colnames(JS)=c('Loc.','Lmean','Lsd');JS
  JS[,1]=Xjs[,1];JS[,2]=log(Xjs[,2]);JS[,3]=Xjs[,3];JS[1:3,]
#-------------------------------------------------------- Fill temp file
  idum<-matrix(nrow=0,ncol=3); colnames(idum)=c('int','diff.sd','rate.sd')
  for (i in 1:(length(JS[,1])-1)){
    for (j in (i+1):length(JS[,1])){
      dum<-numeric(length=3)
      dum[1]=111
      dum[2]=(JS[j,2]-JS[i,2])/JS[i,3]
      dum[3]=dum[2]/dum[1]
      idum=rbind(idum,dum)
    }
  }
  idum[1:3,];length(idum[,1])

#---------------------------------------------------------------------
  idr<-matrix(nrow=length(idum[,1]),ncol=8)
    colnames(idr)=c('int','diff.sd','rate.sd','log.i','log.|d|','log.|r|',
      'sbn','wgt')
  for (i in 1:length(idum[,1])){
    idr[i,1]=111
    idr[i,2]=idum[i,2]
    idr[i,3]=abs(idum[i,3])
    idr[i,4]=log10(idr[i,1]);idr[i,5]=log10(abs(idr[i,2]))
    idr[i,6]=log10(idr[i,3]);idr[i,7]=2;idr[i,8]=1/idr[i,1]
  }
  idr[1:3,]
  idrx=idr[!rowSums(!is.finite(idr)),]  #remove rows that have -Inf
  print(c(nrow(idr),nrow(idrx)))
  writefile<-paste("C://R_aaROEVchapt08//",
    filename[run],"_out.csv",sep="")
  write.csv(idrx,file=writefile,na="NA",row.names=FALSE) 
}
