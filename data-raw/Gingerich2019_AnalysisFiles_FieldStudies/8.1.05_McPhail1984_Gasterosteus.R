##Philip D. Gingerich RATES OF EVOLUTION (Cambridge University Press 2019) 
##8.1.05_McPhail1984_Gasterosteus
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
file1<-"C://R_aaROEVchapt08//8.1.05_McPhail1984_Gasterosteus.csv" 
M<-read.csv(file1,header=TRUE,sep=",",quote="\"",dec=".",fill=TRUE)
M[1:4,]
ncol(M)
#------------------------------------------------------Logged matrix
LM<-matrix(nrow=nrow(M),ncol=ncol(M))
  colnames(LM)=colnames(M)
for (cn in 1:3){							#column no.
  for (rn in 1:nrow(M)){ 					#row no.
    LM[rn,cn]=M[rn,cn]	
  }
}
for (cn in seq(4,32,2)){
  for (rn in 1:nrow(M)){
    LM[rn,cn]=log(M[rn,cn])
    LM[rn,cn+1]=M[rn,cn+1]/M[rn,cn]
  }
}  
#================================================== Calculate rates 
gt=2									#generation time (yr)
idr=matrix(nrow=6*2*15,ncol=9)
  colnames(idr)=c('int','diff.sd','rate.sd.g','log.i','log|d|','log|r|',
    'sbn','wgt','logI+logR')
nc=0
for (cn in seq(4,32,2)){						#column no.
  for (rn in 1:3){ 							#row no.
    for (i in 1:(4-rn)){						#interval yrs
     nc=nc+1
      k=LM[rn+i,2]-LM[rn,2]
      if(k<=gt){idr[nc,1]=1}else{idr[nc,1]=k/gt}		#interval
      meandiff=LM[(rn+i),cn]-LM[rn,cn]				#mean diff.
      psd=PoolSD(LM[rn,3],LM[rn+i,3],LM[rn,cn+1],LM[rn+i,cn+1])#n1,n2,sd1,sd2
      idr[nc,2]=meandiff/psd						#diff.sd
      idr[nc,3]=abs(idr[nc,2])/idr[nc,1]				#rate.sd.gen
      idr[nc,4]=log10(idr[nc,1])					#log.i
      idr[nc,5]=log10(abs(idr[nc,2]))				#log.d
      idr[nc,6]=log10(idr[nc,3])					#log.r
      if(i==1){
        if(idr[nc,1]<=1){idr[nc,7]=1}else{idr[nc,7]=2}		#sbn
      }else{idr[nc,7]=3}
      idr[nc,8]=1/idr[nc,1]						#wgt
      idr[nc,9]=idr[nc,4]+idr[nc,6]					#sum
    }
  }
  for (rn in 5:7){ 							#row no.
    for (i in 1:(8-rn)){							#interval yrs
    nc=nc+1
      k=LM[rn+i,2]-LM[rn,2]
      if(k<=gt){idr[nc,1]=1}else{idr[nc,1]=k/gt}			#interval
      meandiff=LM[(rn+i),cn]-LM[rn,cn]				#mean diff.
      psd=PoolSD(LM[rn,3],LM[rn+i,3],LM[rn,cn+1],LM[rn+i,cn+1])#n1,n2,sd1,sd2
      idr[nc,2]=meandiff/psd						#diff.sd
      idr[nc,3]=abs(idr[nc,2])/idr[nc,1]				#rate.sd.gen
      idr[nc,4]=log10(idr[nc,1])					#log.i
      idr[nc,5]=log10(abs(idr[nc,2]))				#log.d
      idr[nc,6]=log10(idr[nc,3])					#log.r
      if(i==1){
        if(idr[nc,1]<=1){idr[nc,7]=1}else{idr[nc,7]=2}		#sbn
      }else{idr[nc,7]=3}
      idr[nc,8]=1/idr[nc,1]						#wgt
      idr[nc,9]=idr[nc,4]+idr[nc,6]					#sum
    }
  }
}
plot(idr[,4],idr[,6])
summary(idr[,9])

#--------------------------------------------------------------------------
idrx=idr[!rowSums(!is.finite(idr)),]  #remove rows that have -Inf
idrxM=idrx[,1:8]
idrc=idrx[idrx[,9]>=1,]  #keep only rows where interval >= gentime
idrcM=idrc[,1:8]
nrow(idr);nrow(idrx);nrow(idrxM);nrow(idrcM)
writefile<-"C://R_aaROEVchapt08//8.1.05_McPhail1984_Gasterosteus_out.csv"
write.csv(idrxM,file=writefile,na="NA",row.names=FALSE) 


