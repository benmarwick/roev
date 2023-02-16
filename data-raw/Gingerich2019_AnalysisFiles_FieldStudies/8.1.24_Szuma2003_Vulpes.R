##Philip D. Gingerich RATES OF EVOLUTION (Cambridge University Press 2019) 
##8.1.24_Szuma2003_Vulpes
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
file1<-"C://R_aaROEVchapt08//8.1.24_Szuma2003_Vulpes.csv" 
S<-read.csv(file1,header=TRUE,sep=",",quote="\"",dec=".",fill=TRUE)
S[1:4,]
nrow(S)
#------------------------------------------------------Logged matrix
LS<-S[,c(2:6)]
  colnames(LS)=c('trait','year','N','ln mean','ln sd');LS[1,]
for (i in 1:nrow(S)){
  LS[i,4]=log(S[i,5])
  LS[i,5]=S[i,6]/S[i,5]
};LS[1:5,]
#================================================== Calculate rates 
gt=4				#Vulpes generation time following Kerk et al. 2013
idr=matrix(nrow=0,ncol=9)
  colnames(idr)=c('int','diff.sd','rate.sd.g','log.i','log|d|','log|r|',
    'sbn','wgt','tr.no.')
for(tn in 1:8){								#trait number
  for (sr in (5*tn-3):(5*tn-1)){ 					#start row
    for (i in 1:(5*tn-sr)){							#interval steps
      id<-numeric(length=9)						#temp vector
      k=LS[sr+i,2]-LS[sr,2]					#interval in yrs
      if(k<=gt){id[1]=1}else{id[1]=k/gt}				#interval in gen
      meandiff=LS[(sr+i),4]-LS[sr,4]				#mean diff.
      psd=PoolSD(LS[sr,3],LS[sr+i,3],LS[sr,5],LS[sr+i,5])#n1,n2,sd1,sd2
      id[2]=meandiff/psd							#diff.sd
      id[3]=abs(id[2])/id[1]						#rate.sd.gen
      id[4]=log10(id[1])							#log.i
      id[5]=log10(abs(id[2]))						#log.d
      id[6]=log10(id[3])							#log.r
      if(i==1){
        if(id[1]<=1){id[7]=1}else{id[7]=2}				#sbn
      }else{id[7]=3}
      id[8]=1/id[1]							#wgt
      id[9]=tn #id[4]+id[6]						#sum
      idr=rbind(idr,id)
    }
  }
}
nrow(idr)
plot(idr[,4],idr[,6])
min(idr[,9]);max(idr[,9])

#--------------------------------------------------------------------------
idrx=idr[!rowSums(!is.finite(idr)),]  #remove rows that have -Inf
idrxS=idrx[,1:8]
idrc=idrx[idrx[,9]>=1,]  #keep only rows where interval >= gentime
idrcS=idrc[,1:8]
nrow(idr);nrow(idrx);nrow(idrxS);nrow(idrcS)
writefile<-"C://R_aaROEVchapt08//8.1.24_Szuma2003_Vulpes_out.csv"
write.csv(idrx,file=writefile,na="NA",row.names=FALSE) 


