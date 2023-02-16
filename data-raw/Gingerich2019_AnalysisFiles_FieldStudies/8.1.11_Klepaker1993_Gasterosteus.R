##Philip D. Gingerich RATES OF EVOLUTION (Cambridge University Press 2019) 
##8.1.11_Klepaker1993_Gasterosteus
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
file1<-"C://R_aaROEVchapt08//8.1.11_Klepaker1993_Gasterosteus.csv" 
K<-read.csv(file1,header=TRUE,sep=",",quote="\"",dec=".",fill=TRUE)
K[1:4,]
#------------------------------------------------------Logged matrix
LK<-matrix(nrow=nrow(K),ncol=12)
  colnames(LK)=c('I60','N60','M60','S60','I87','N87','M87','S87',
    'I91','N91','M91','S91');LK[1,]
for (i in 1:nrow(K)){
  LK[i,1]=K[i,2]					#year
  LK[i,2]=K[i,3]					#N
  LK[i,3]=log(K[i,4])				#ln mean
  LK[i,4]=K[i,5]/K[i,4]				#ln stdev

  LK[i,5]=K[i,6]					#year
  LK[i,6]=K[i,7]					#N
  LK[i,7]=log(K[i,8])				#ln mean
  LK[i,8]=K[i,9]/K[i,8]				#ln stdev

  LK[i,9]=K[i,10]					#year
  LK[i,10]=K[i,11]					#N
  LK[i,11]=log(K[i,12])				#ln mean
  LK[i,12]=K[i,13]/K[i,12]				#ln stdev
}
LK[1:4,]
#================================================== Calculate rates in triplets
gentime=2
idr=matrix(nrow=3*nrow(LK),ncol=9)
  colnames(idr)=c('int','diff.sd','rate.sd.g','log.i','log|d|','log|r|',
    'sbn','wgt','logI+logR')
nc=0
for (i in 1:3){			#no of columnar comparisons
  if(i==1){c1=1;c2=2;c3=3;c4=4;c5=5;c6=6;c7=7;c8=8}
  if(i==2){c1=1;c2=2;c3=3;c4=4;c5=9;c6=10;c7=11;c8=12}
  if(i==3){c1=5;c2=6;c3=7;c4=8;c5=9;c6=10;c7=11;c8=12}
  for (j in 1:nrow(LK)){
    nc=nc+1
    idr[nc,1]=(LK[nc-(i-1)*24,c5]-LK[nc-(i-1)*24,c1])/gentime#int.g (one way)
    diff=LK[nc-(i-1)*24,c7]-LK[nc-(i-1)*24,c3]
    poolsd=PoolSD(LK[nc-(i-1)*24,c2],LK[nc-(i-1)*24,c6],	#n1,n2,sd1,sd2
      LK[nc-(i-1)*24,c4],LK[nc-(i-1)*24,c8])				
    idr[nc,2]=diff/poolsd						#diff.sd
    idr[nc,3]=abs(idr[nc,2])/idr[nc,1]				#rate.sd.g
    idr[nc,4]=log10(idr[nc,1])					#log.I
    idr[nc,5]=log10(abs(idr[nc,2]))					#log.D
    idr[nc,6]=log10(idr[nc,3])					#log.R
    idr[nc,7]=2								#sbn
    idr[nc,8]=1/idr[nc,1]						#wgt
    idr[nc,9]=idr[nc,4]+idr[nc,6]					#Log.I+Log.R
  }
}
idrx=idr[!rowSums(!is.finite(idr)),]  #remove rows that have -Inf
nrow(idr);nrow(idrx)
plot(idr[,4],idr[,6])

#--------------------------------------------------------------------------
writefile<-"C://R_aaROEVchapt08//8.1.11_Klepaker1993_Gasterosteus_out.csv"
write.csv(idrx,file=writefile,na="NA",row.names=FALSE) 


