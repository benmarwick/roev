##Philip D. Gingerich RATES OF EVOLUTION (Cambridge University Press 2019) 
##9.2.22_Forsten1990_Equus
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
file1<-"C://R_aaROEVchapt09//9.2.22_Forsten1990_Equus_Equus1.csv" 
file2<-"C://R_aaROEVchapt09//9.2.22_Forsten1990_Equus_Equus2.csv" 
#========================================================================== 
idr=matrix(nrow=0,ncol=9)
  colnames(idr)=c('int.g','diff.sd','rate.sd.g','log.i','log|d|','log|r|',
    'sbn','wgt','li+lr')
#----------------------------------------------------------Equus1 rates
#                                   convert data.frame row to vector
F<-read.csv(file1,header=TRUE,sep=",",quote="\"",dec=".",fill=TRUE)
F;nrow(F)
age<-apply(as.matrix.noquote(F[2,]),2,as.numeric);age
#------------------------------------------------------Logged matrix
F[3,]
FF=apply(as.matrix.noquote(F[3:nrow(F),]),2,as.numeric);FF[1,]
LF<-matrix(nrow=ncol(FF),ncol=nrow(FF))
for (cn in 1:7){
  for (rn in seq(1,ncol(LF),3)){
    LF[cn,rn-2]=FF[rn-2,cn]
    LF[cn,rn-1]=FF[rn-1,cn]
    LF[cn,rn]=FF[rn,cn]
  }
}
#---------------------------------------------------------------- rates
gt=8	#generation time 8 years fide Berger 1986, fig. 5.2
nr=nrow(LF)					#number of rows
nn=.5*(nr-1)*nr;nn1=(ncol(LF)/3)*nn
nc=0
for (c in seq(3,ncol(LF),3)){		#columns to analyze
  for (i in 1:(nr-1)){			#increment
    for (sr in 1:(nr-i)){		#start row
      id<-numeric(length=9)
      nc=nc+1
      iyr=age[sr+i]-age[sr]	#interval in years
      id[1]=iyr/gt					#interval in generations
      meandiff=LF[(sr),c-1]-LF[sr+i,c-1]		#mean diff.
      psd=PoolSD(LF[sr+i,c-2],LF[sr,c-2],LF[sr+i,c],LF[sr,c])#n1,n2,sd1,sd2
      id[2]=meandiff/psd				#diff.sd
      id[3]=abs(id[2])/id[1]				#rate.sd.gen
      id[4]=log10(id[1])				#log.i
      id[5]=log10(abs(id[2]))				#log.d
      id[6]=log10(id[3])				#log.r
      if(i==1){id[7]=2}else{id[7]=3}		#sbn
      id[8]=1/id[1]					#wgt
      id[9]=id[4]+id[6]					#sum
    	idr=rbind(idr,id)
    }
  }
}
nn1;nrow(idr)
#----------------------------------------------------------Equus2 rates
#                                   convert data.frame row to vector
F<-read.csv(file2,header=TRUE,sep=",",quote="\"",dec=".",fill=TRUE)
F;nrow(F)
age<-apply(as.matrix.noquote(F[2,]),2,as.numeric);age
#------------------------------------------------------Logged matrix
F[3,]
FF=apply(as.matrix.noquote(F[3:nrow(F),]),2,as.numeric);FF[1,]
LF<-matrix(nrow=ncol(FF),ncol=nrow(FF))
for (cn in 1:7){
  for (rn in seq(1,ncol(LF),3)){
    LF[cn,rn-2]=FF[rn-2,cn]
    LF[cn,rn-1]=FF[rn-1,cn]
    LF[cn,rn]=FF[rn,cn]
  }
}
#---------------------------------------------------------------- rates
gt=8	#generation time 8 years fide Berger 1986, fig. 5.2
nr=nrow(LF)					#number of rows
nn=.5*(nr-1)*nr;nn2=(ncol(LF)/3)*nn
nc=0
for (c in seq(3,ncol(LF),3)){		#columns to analyze
  for (i in 1:(nr-1)){			#increment
    for (sr in 1:(nr-i)){		#start row
      id<-numeric(length=9)
      nc=nc+1
      iyr=age[sr+i]-age[sr]	#interval in years
      id[1]=iyr/gt					#interval in generations
      meandiff=LF[(sr),c-1]-LF[sr+i,c-1]		#mean diff.
      psd=PoolSD(LF[sr+i,c-2],LF[sr,c-2],LF[sr+i,c],LF[sr,c])#n1,n2,sd1,sd2
      id[2]=meandiff/psd				#diff.sd
      id[3]=abs(id[2])/id[1]				#rate.sd.gen
      id[4]=log10(id[1])				#log.i
      id[5]=log10(abs(id[2]))				#log.d
      id[6]=log10(id[3])				#log.r
      if(i==1){id[7]=2}else{id[7]=3}		#sbn
      id[8]=1/id[1]					#wgt
      id[9]=id[4]+id[6]					#sum
    	idr=rbind(idr,id)
    }
  }
}
#=========================================================================
nrow(idr)
plot(idr[,4],idr[,6])
#--------------------------------------------------------------------------
idrx=idr[!rowSums(!is.finite(idr)),]  #remove rows that have -Inf
nn1+nn2;nrow(idr);nrow(idrx)
min(idrx[,9]);max(idrx[,9])
writefile<-"C://R_aaROEVchapt09//9.2.22_Forsten1990_Equus_out.csv"
write.csv(idrx,file=writefile,na="NA",row.names=FALSE) 


