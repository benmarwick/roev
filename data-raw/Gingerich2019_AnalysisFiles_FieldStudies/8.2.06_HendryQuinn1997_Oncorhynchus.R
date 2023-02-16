##Philip D. Gingerich RATES OF EVOLUTION (Cambridge University Press 2019) 
##8.2.06_HendryQuinn1997_Oncorhynchus
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
file1<-"C://R_aaROEVchapt08//8.2.06_HendryQuinn1997_Oncorhynchus.csv" 
H<-read.csv(file1,header=TRUE,sep=",",quote="\"",dec=".",fill=TRUE)
H[1:4,]
#------------------------------------------------------Logged matrix
LH<-matrix(nrow=nrow(H),ncol=8)
  colnames(LH)=c('I92','N92','M92','S92','I93','N93','M93','S93');LH[1,]
for (i in 1:nrow(H)){
  LH[i,1]=(H[i,5]-H[i,4])/H[i,3]
  LH[i,2]=H[i,6]
  LH[i,3]=log(H[i,7])
  LH[i,4]=(sqrt(H[i,6])*H[i,8])/H[i,7]

  LH[i,5]=(H[i,9]-H[i,4])/H[i,3]
  LH[i,6]=H[i,10]
  LH[i,7]=log(H[i,11])
  LH[i,8]=(sqrt(H[i,10])*H[i,12])/H[i,11]
}
LH[1:4,]
#================================================== Calculate rates in triplets
idr=matrix(nrow=2*nrow(LH),ncol=9)
  colnames(idr)=c('int','diff.sd','rate.sd.g','log.i','log|d|','log|r|',
    'sbn','wgt','logI+logR')
nc=0
for (j in 1:5){		#no of triplets in 1947-1992 comparisons
  for (k in 1:3){
    if (k==1){m=1;n=2}
    if (k==2){m=2;n=3}
    if (k==3){m=1;n=3}
    nc=nc+1
    idr[nc,1]=2*LH[3*(j-1)+k,1]				#int.g (two way)
    diff=abs(LH[3*(j-1)+m,3]-LH[3*(j-1)+n,3])
    poolsd=PoolSD(LH[3*(j-1)+m,2],LH[3*(j-1)+n,2],	#n1,n2,sd1,sd2
    LH[3*(j-1)+m,4],LH[3*(j-1)+n,4])				
    idr[nc,2]=diff/poolsd					#diff.sd
    idr[nc,3]=idr[nc,2]/idr[nc,1]				#rate.sd.g
    idr[nc,4]=log10(idr[nc,1])				#log.I
    idr[nc,5]=log10(idr[nc,2])				#log.D
    idr[nc,6]=log10(idr[nc,3])				#log.R
    idr[nc,7]=2							#sbn
    idr[nc,8]=1/idr[nc,1]					#wgt
    idr[nc,9]=idr[nc,4]+idr[nc,6]				#Log.I+Log.R
  }
}
for (j in 1:5){		#no of triplets in 1947-1993 comparisons
  for (k in 1:3){
    if (k==1){m=1;n=2}
    if (k==2){m=2;n=3}
    if (k==3){m=1;n=3}
    nc=nc+1
    idr[nc,1]=2*LH[3*(j-1)+k,5]				#int.g (two way)
    diff=abs(LH[3*(j-1)+m,7]-LH[3*(j-1)+n,7])
    poolsd=PoolSD(LH[3*(j-1)+m,6],LH[3*(j-1)+n,6],	#n1,n2,sd1,sd2
    LH[3*(j-1)+m,8],LH[3*(j-1)+n,8])				
    idr[nc,2]=diff/poolsd					#diff.sd
    idr[nc,3]=idr[nc,2]/idr[nc,1]				#rate.sd.g
    idr[nc,4]=log10(idr[nc,1])				#log.I
    idr[nc,5]=log10(idr[nc,2])				#log.D
    idr[nc,6]=log10(idr[nc,3])				#log.R
    idr[nc,7]=2							#sbn
    idr[nc,8]=1/idr[nc,1]					#wgt
    idr[nc,9]=idr[nc,4]+idr[nc,6]				#Log.I+Log.R
  }
}
idr


#--------------------------------------------------------------------------
writefile<-"C://R_aaROEVchapt08//8.2.06_HendryQuinn1997_Oncorhynchus_out.csv"
write.csv(idr,file=writefile,na="NA",row.names=FALSE) 


