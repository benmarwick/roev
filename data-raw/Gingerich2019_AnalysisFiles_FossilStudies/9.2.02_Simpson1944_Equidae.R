##Philip D. Gingerich RATES OF EVOLUTION (Cambridge University Press 2019) 
##9.2.02_Simpson1944_Equidae
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
file1<-"C://R_aaROEVchapt09//9.2.02_Simpson1944_Equidae.csv" 
S<-read.csv(file1,header=TRUE,sep=",",quote="\"",dec=".",fill=TRUE)
S;nrow(S)
#------------------------------------------------------Logged matrix
LS=S;colnames(LS)=c('line','genspec','McF','McF_age','wgt_kg','N',
  'LnPm','LnPs','LnEm','LnEs','LnHm','Lnhs');LS[1,]
LS[,7]=log(S[,7]);LS[,8]=S[,8]/S[,7]
LS[,9]=log(S[,9]);LS[,10]=S[,10]/S[,9]
LS[,11]=log(S[,11]);LS[,12]=S[,12]/S[,11]
S;LS
#========================================================================== 
idr=matrix(nrow=0,ncol=9)
  colnames(idr)=c('int.g','diff.sd','rate.sd.g','log.i','log|d|','log|r|',
    'sbn','wgt','li+lr')
#----------------------------------------intervals in line 1:paracone hgt
col=7					#ln paracone hgt
nr=4					#number of rows
id<-numeric(length=9)
nc=0					#n count
for (i in 1:(nr-1)){		#increment
  for (sr in 1:(nr-i)){		#start row
    nc=nc+1
    iyr=1000000*(LS[(sr+i),4]-LS[sr,4])	#interval in years
    w1g=1000*LS[sr,5];w2g=1000*LS[(sr+i),5];wg=exp(.5*(log(w1g)+log(w2g)))
    gt=10^(.266*log10(wg)-.553)		#Eq 3.2 generation time in years
    id[1]=iyr/gt					#interval in generations
    meandiff=LS[(sr+i),col]-LS[sr,col]	#mean diff.
    psd=PoolSD(LS[sr,6],LS[sr+i,6],LS[sr,col+1],LS[sr+i,col+1])#n1,n2,sd1,sd2
    id[2]=meandiff/psd				#diff.sd
    id[3]=abs(id[2])/id[1]			#rate.sd.gen
    id[4]=log10(id[1])				#log.i
    id[5]=log10(abs(id[2]))			#log.d
    id[6]=log10(id[3])				#log.r
    if(i==1){id[7]=2}else{id[7]=3}		#sbn
    id[8]=1/id[1]					#wgt
    id[9]=id[4]+id[6]				#sum
    idr=rbind(idr,id)
  }
}
#----------------------------------------intervals in line 1:ectoloph len.
col=9					#ln paracone hgt
nr=4					#number of rows
id<-numeric(length=9)
nc=0					#n count
for (i in 1:(nr-1)){		#increment
  for (sr in 1:(nr-i)){		#start row
    nc=nc+1
    iyr=1000000*(LS[(sr+i),4]-LS[sr,4])	#interval in years
    w1g=1000*LS[sr,5];w2g=1000*LS[(sr+i),5];wg=exp(.5*(log(w1g)+log(w2g)))
    gt=10^(.266*log10(wg)-.553)		#Eq 3.2 generation time in years
    id[1]=iyr/gt					#interval in generations
    meandiff=LS[(sr+i),col]-LS[sr,col]	#mean diff.
    psd=PoolSD(LS[sr,6],LS[sr+i,6],LS[sr,col+1],LS[sr+i,col+1])#n1,n2,sd1,sd2
    id[2]=meandiff/psd				#diff.sd
    id[3]=abs(id[2])/id[1]			#rate.sd.gen
    id[4]=log10(id[1])				#log.i
    id[5]=log10(abs(id[2]))			#log.d
    id[6]=log10(id[3])				#log.r
    if(i==1){id[7]=2}else{id[7]=3}		#sbn
    id[8]=1/id[1]					#wgt
    id[9]=id[4]+id[6]				#sum
    idr=rbind(idr,id)
  }
}
#----------------------------------------intervals in line 1:hypsodonty
col=11				#ln paracone hgt
nr=4					#number of rows
id<-numeric(length=9)
nc=0					#n count
for (i in 1:(nr-1)){		#increment
  for (sr in 1:(nr-i)){		#start row
    nc=nc+1
    iyr=1000000*(LS[(sr+i),4]-LS[sr,4])	#interval in years
    w1g=1000*LS[sr,5];w2g=1000*LS[(sr+i),5];wg=exp(.5*(log(w1g)+log(w2g)))
    gt=10^(.266*log10(wg)-.553)		#Eq 3.2 generation time in years
    id[1]=iyr/gt					#interval in generations
    meandiff=LS[(sr+i),col]-LS[sr,col]	#mean diff.
    psd=PoolSD(LS[sr,6],LS[sr+i,6],LS[sr,col+1],LS[sr+i,col+1])#n1,n2,sd1,sd2
    id[2]=meandiff/psd				#diff.sd
    id[3]=abs(id[2])/id[1]			#rate.sd.gen
    id[4]=log10(id[1])				#log.i
    id[5]=log10(abs(id[2]))			#log.d
    id[6]=log10(id[3])				#log.r
    if(i==1){id[7]=2}else{id[7]=3}		#sbn
    id[8]=1/id[1]					#wgt
    id[9]=id[4]+id[6]				#sum
    idr=rbind(idr,id)
  }
}
#----------------------------------------intervals in line 2:paracone hgt
col=7					#ln paracone hgt
nr=3					#number of rows
id<-numeric(length=9)
nc=0					#n count
for (i in 1:(nr-1)){		#increment
  for (sr in (1+4):((nr-i)+4)){		#start row
    nc=nc+1
    iyr=1000000*(LS[(sr+i),4]-LS[sr,4])	#interval in years
    w1g=1000*LS[sr,5];w2g=1000*LS[(sr+i),5];wg=exp(.5*(log(w1g)+log(w2g)))
    gt=10^(.266*log10(wg)-.553)		#Eq 3.2 generation time in years
    id[1]=iyr/gt					#interval in generations
    meandiff=LS[(sr+i),col]-LS[sr,col]	#mean diff.
    psd=PoolSD(LS[sr,6],LS[sr+i,6],LS[sr,col+1],LS[sr+i,col+1])#n1,n2,sd1,sd2
    id[2]=meandiff/psd				#diff.sd
    id[3]=abs(id[2])/id[1]			#rate.sd.gen
    id[4]=log10(id[1])				#log.i
    id[5]=log10(abs(id[2]))			#log.d
    id[6]=log10(id[3])				#log.r
    if(i==1){id[7]=2}else{id[7]=3}		#sbn
    id[8]=1/id[1]					#wgt
    id[9]=id[4]+id[6]				#sum
    if(nc>1){idr=rbind(idr,id)}
  }
}
#----------------------------------------intervals in line 2:ectoloph len.
col=9					#ln paracone hgt
nr=3					#number of rows
id<-numeric(length=9)
nc=0					#n count
for (i in 1:(nr-1)){		#increment
  for (sr in (1+4):((nr-i)+4)){		#start row
    nc=nc+1
    iyr=1000000*(LS[(sr+i),4]-LS[sr,4])	#interval in years
    w1g=1000*LS[sr,5];w2g=1000*LS[(sr+i),5];wg=exp(.5*(log(w1g)+log(w2g)))
    gt=10^(.266*log10(wg)-.553)		#Eq 3.2 generation time in years
    id[1]=iyr/gt					#interval in generations
    meandiff=LS[(sr+i),col]-LS[sr,col]	#mean diff.
    psd=PoolSD(LS[sr,6],LS[sr+i,6],LS[sr,col+1],LS[sr+i,col+1])#n1,n2,sd1,sd2
    id[2]=meandiff/psd				#diff.sd
    id[3]=abs(id[2])/id[1]			#rate.sd.gen
    id[4]=log10(id[1])				#log.i
    id[5]=log10(abs(id[2]))			#log.d
    id[6]=log10(id[3])				#log.r
    if(i==1){id[7]=2}else{id[7]=3}		#sbn
    id[8]=1/id[1]					#wgt
    id[9]=id[4]+id[6]				#sum
    if(nc>1){idr=rbind(idr,id)}
  }
}
#----------------------------------------intervals in line 2:hypsodonty
col=11					#ln paracone hgt
nr=3					#number of rows
id<-numeric(length=9)
nc=0					#n count
for (i in 1:(nr-1)){		#increment
  for (sr in (1+4):((nr-i)+4)){		#start row
    nc=nc+1
    iyr=1000000*(LS[(sr+i),4]-LS[sr,4])	#interval in years
    w1g=1000*LS[sr,5];w2g=1000*LS[(sr+i),5];wg=exp(.5*(log(w1g)+log(w2g)))
    gt=10^(.266*log10(wg)-.553)		#Eq 3.2 generation time in years
    id[1]=iyr/gt					#interval in generations
    meandiff=LS[(sr+i),col]-LS[sr,col]	#mean diff.
    psd=PoolSD(LS[sr,6],LS[sr+i,6],LS[sr,col+1],LS[sr+i,col+1])#n1,n2,sd1,sd2
    id[2]=meandiff/psd				#diff.sd
    id[3]=abs(id[2])/id[1]			#rate.sd.gen
    id[4]=log10(id[1])				#log.i
    id[5]=log10(abs(id[2]))			#log.d
    id[6]=log10(id[3])				#log.r
    if(i==1){id[7]=2}else{id[7]=3}		#sbn
    id[8]=1/id[1]					#wgt
    id[9]=id[4]+id[6]				#sum
    if(nc>1){idr=rbind(idr,id)}
  }
}
#---------------------------------------------------------------------
nrow(idr)
plot(idr[,4],idr[,6])
min(idr[,9]);max(idr[,9])
#--------------------------------------------------------------------------
idrx=idr[!rowSums(!is.finite(idr)),]  #remove rows that have -Inf
nrow(idr);nrow(idrx)
writefile<-"C://R_aaROEVchapt09//9.2.02_Simpson1944_Equidae_out.csv"
write.csv(idrx,file=writefile,na="NA",row.names=FALSE) 


