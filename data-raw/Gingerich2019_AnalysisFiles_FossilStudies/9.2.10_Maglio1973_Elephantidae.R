##Philip D. Gingerich RATES OF EVOLUTION (Cambridge University Press 2019) 
##9.2.10_Maglio1973_Elephantidae
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
file1<-"C://R_aaROEVchapt09//9.2.10_Maglio1973_Elephantidae.csv" 
M<-read.csv(file1,header=TRUE,sep=",",quote="\"",dec=".",fill=TRUE)
nrow(M);M[1:20,1:10]
#------------------------------------------------------Logged matrix
LM<-matrix(nrow=18*6,ncol=3+7*3)
  colnames(LM)=c('spec','age','tooth','C1N','C1M','C1S','C2N','C2M','C2S',
    'C3N','C3M','C3S','C4N','C4M','C4S','C5N','C5M','C5S','C6N','C6M','C6S',
    'C7N','C7M','C7S');#LM[1,]
for (s in 1:18){					#species
  for (r in 1:6){					#row within species
    LM[6*(s-1)+r,1]=M[18*(s-1)+(3*r-2),1]
    LM[6*(s-1)+r,2]=M[18*(s-1)+(3*r-2),4]
    LM[6*(s-1)+r,3]=r
    for (c in 1:7){				#column within row
      LM[6*(s-1)+r,3*c+1]=M[18*(s-1)+(3*r),c+6]
      LM[6*(s-1)+r,3*c+2]=log(M[18*(s-1)+(3*r-2),c+6])
      LM[6*(s-1)+r,3*c+3]=M[18*(s-1)+(3*r-1),c+6]/M[18*(s-1)+(3*r-2),c+6]
    }
  }
}  #LM[1:12,1:9]
#========================================================================== 
idr=matrix(nrow=0,ncol=9)
  colnames(idr)=c('int.g','diff.sd','rate.sd.g','log.i','log|d|','log|r|',
    'sbn','wgt','li+lr')
#------------------------------------------------------ species combinations
#------------------------------- P. gomphotheroides to L. atlantica (1,2,3)
k=c(1,2,2)	#species a-b combination with sbn
idt=Maglio73(LM[(k[1]*6-5):(k[1]*6),],LM[(k[2]*6-5):(k[2]*6),],k[3])
idr=rbind(idr,idt)
k=c(2,3,2)	#species a-b combination with sbn
idt=Maglio73(LM[(k[1]*6-5):(k[1]*6),],LM[(k[2]*6-5):(k[2]*6),],k[3])
idr=rbind(idr,idt)
k=c(1,3,3)	#species a-b combination with sbn
idt=Maglio73(LM[(k[1]*6-5):(k[1]*6),],LM[(k[2]*6-5):(k[2]*6),],k[3])
idr=rbind(idr,idt)
plot(idr[,4],idr[,6]);nrow(idr)
#------------------------------- P. gomphotheroides to L. africana (1,2,4)
k=c(2,4,2)	#species a-b combination with sbn
idt=Maglio73(LM[(k[1]*6-5):(k[1]*6),],LM[(k[2]*6-5):(k[2]*6),],k[3])
idr=rbind(idr,idt)
k=c(1,4,3)	#species a-b combination with sbn
idt=Maglio73(LM[(k[1]*6-5):(k[1]*6),],LM[(k[2]*6-5):(k[2]*6),],k[3])
idr=rbind(idr,idt)
plot(idr[,4],idr[,6]);nrow(idr)
#------------------------------- P. gomphotheroides to E. namadicus (1,5,6,9)
k=c(1,5,2)	#species a-b combination with sbn
idt=Maglio73(LM[(k[1]*6-5):(k[1]*6),],LM[(k[2]*6-5):(k[2]*6),],k[3])
idr=rbind(idr,idt)
k=c(5,6,2)	#species a-b combination with sbn
idt=Maglio73(LM[(k[1]*6-5):(k[1]*6),],LM[(k[2]*6-5):(k[2]*6),],k[3])
idr=rbind(idr,idt)
k=c(6,9,2)	#species a-b combination with sbn
idt=Maglio73(LM[(k[1]*6-5):(k[1]*6),],LM[(k[2]*6-5):(k[2]*6),],k[3])
idr=rbind(idr,idt)
k=c(1,6,3)	#species a-b combination with sbn
idt=Maglio73(LM[(k[1]*6-5):(k[1]*6),],LM[(k[2]*6-5):(k[2]*6),],k[3])
idr=rbind(idr,idt)
k=c(1,9,3)	#species a-b combination with sbn
idt=Maglio73(LM[(k[1]*6-5):(k[1]*6),],LM[(k[2]*6-5):(k[2]*6),],k[3])
idr=rbind(idr,idt)
k=c(1,9,3)	#species a-b combination with sbn
idt=Maglio73(LM[(k[1]*6-5):(k[1]*6),],LM[(k[2]*6-5):(k[2]*6),],k[3])
idr=rbind(idr,idt)
plot(idr[,4],idr[,6]);nrow(idr)
#----------------------------- P. gomphotheroides to E. iolensis (1,5,6,7,8)
k=c(6,7,2)	#species a-b combination with sbn
idt=Maglio73(LM[(k[1]*6-5):(k[1]*6),],LM[(k[2]*6-5):(k[2]*6),],k[3])
idr=rbind(idr,idt)
k=c(7,8,2)	#species a-b combination with sbn
idt=Maglio73(LM[(k[1]*6-5):(k[1]*6),],LM[(k[2]*6-5):(k[2]*6),],k[3])
idr=rbind(idr,idt)
k=c(1,7,3)	#species a-b combination with sbn
idt=Maglio73(LM[(k[1]*6-5):(k[1]*6),],LM[(k[2]*6-5):(k[2]*6),],k[3])
idr=rbind(idr,idt)
k=c(1,8,3)	#species a-b combination with sbn
idt=Maglio73(LM[(k[1]*6-5):(k[1]*6),],LM[(k[2]*6-5):(k[2]*6),],k[3])
idr=rbind(idr,idt)
k=c(5,7,3)	#species a-b combination with sbn
idt=Maglio73(LM[(k[1]*6-5):(k[1]*6),],LM[(k[2]*6-5):(k[2]*6),],k[3])
idr=rbind(idr,idt)
k=c(5,8,3)	#species a-b combination with sbn
idt=Maglio73(LM[(k[1]*6-5):(k[1]*6),],LM[(k[2]*6-5):(k[2]*6),],k[3])
idr=rbind(idr,idt)
k=c(6,8,3)	#species a-b combination with sbn
idt=Maglio73(LM[(k[1]*6-5):(k[1]*6),],LM[(k[2]*6-5):(k[2]*6),],k[3])
idr=rbind(idr,idt)
plot(idr[,4],idr[,6]);nrow(idr)
#----------------------- P. gomphotheroides to E. hysudrindicus (1,5,12,13)
k=c(5,12,2)	#species a-b combination with sbn
idt=Maglio73(LM[(k[1]*6-5):(k[1]*6),],LM[(k[2]*6-5):(k[2]*6),],k[3])
idr=rbind(idr,idt)
k=c(12,13,2)#species a-b combination with sbn
idt=Maglio73(LM[(k[1]*6-5):(k[1]*6),],LM[(k[2]*6-5):(k[2]*6),],k[3])
idr=rbind(idr,idt)
k=c(1,12,3)	#species a-b combination with sbn
idt=Maglio73(LM[(k[1]*6-5):(k[1]*6),],LM[(k[2]*6-5):(k[2]*6),],k[3])
idr=rbind(idr,idt)
k=c(1,13,3)	#species a-b combination with sbn
idt=Maglio73(LM[(k[1]*6-5):(k[1]*6),],LM[(k[2]*6-5):(k[2]*6),],k[3])
idr=rbind(idr,idt)
k=c(5,13,3)	#species a-b combination with sbn
idt=Maglio73(LM[(k[1]*6-5):(k[1]*6),],LM[(k[2]*6-5):(k[2]*6),],k[3])
idr=rbind(idr,idt)
plot(idr[,4],idr[,6]);nrow(idr)
#----------------------- P. gomphotheroides to E. celebensis (1,5,10,11)
k=c(5,10,2)	#species a-b combination with sbn
idt=Maglio73(LM[(k[1]*6-5):(k[1]*6),],LM[(k[2]*6-5):(k[2]*6),],k[3])
idr=rbind(idr,idt)
k=c(10,11,2)#species a-b combination with sbn
idt=Maglio73(LM[(k[1]*6-5):(k[1]*6),],LM[(k[2]*6-5):(k[2]*6),],k[3])
idr=rbind(idr,idt)
k=c(1,10,3)	#species a-b combination with sbn
idt=Maglio73(LM[(k[1]*6-5):(k[1]*6),],LM[(k[2]*6-5):(k[2]*6),],k[3])
idr=rbind(idr,idt)
k=c(1,11,3)	#species a-b combination with sbn
idt=Maglio73(LM[(k[1]*6-5):(k[1]*6),],LM[(k[2]*6-5):(k[2]*6),],k[3])
idr=rbind(idr,idt)
k=c(5,11,3)	#species a-b combination with sbn
idt=Maglio73(LM[(k[1]*6-5):(k[1]*6),],LM[(k[2]*6-5):(k[2]*6),],k[3])
idr=rbind(idr,idt)
plot(idr[,4],idr[,6]);nrow(idr)
#------------------------- P. gomphotheroides to M. africanvus (1,14,15)
k=c(1,14,2)	#species a-b combination with sbn
idt=Maglio73(LM[(k[1]*6-5):(k[1]*6),],LM[(k[2]*6-5):(k[2]*6),],k[3])
idr=rbind(idr,idt)
k=c(14,15,2)#species a-b combination with sbn
idt=Maglio73(LM[(k[1]*6-5):(k[1]*6),],LM[(k[2]*6-5):(k[2]*6),],k[3])
idr=rbind(idr,idt)
k=c(1,15,3)	#species a-b combination with sbn
idt=Maglio73(LM[(k[1]*6-5):(k[1]*6),],LM[(k[2]*6-5):(k[2]*6),],k[3])
idr=rbind(idr,idt)
plot(idr[,4],idr[,6]);nrow(idr)
#------------------------- P. gomphotheroides to M primigenius (1,14,15)
k=c(14,16,2)	#species a-b combination with sbn
idt=Maglio73(LM[(k[1]*6-5):(k[1]*6),],LM[(k[2]*6-5):(k[2]*6),],k[3])
idr=rbind(idr,idt)
k=c(16,17,2)#species a-b combination with sbn
idt=Maglio73(LM[(k[1]*6-5):(k[1]*6),],LM[(k[2]*6-5):(k[2]*6),],k[3])
idr=rbind(idr,idt)
k=c(17,18,2)#species a-b combination with sbn
idt=Maglio73(LM[(k[1]*6-5):(k[1]*6),],LM[(k[2]*6-5):(k[2]*6),],k[3])
idr=rbind(idr,idt)
k=c(1,16,3)	#species a-b combination with sbn
idt=Maglio73(LM[(k[1]*6-5):(k[1]*6),],LM[(k[2]*6-5):(k[2]*6),],k[3])
idr=rbind(idr,idt)
k=c(1,17,3)	#species a-b combination with sbn
idt=Maglio73(LM[(k[1]*6-5):(k[1]*6),],LM[(k[2]*6-5):(k[2]*6),],k[3])
idr=rbind(idr,idt)
k=c(1,18,3)	#species a-b combination with sbn
idt=Maglio73(LM[(k[1]*6-5):(k[1]*6),],LM[(k[2]*6-5):(k[2]*6),],k[3])
idr=rbind(idr,idt)
k=c(14,17,3)	#species a-b combination with sbn
idt=Maglio73(LM[(k[1]*6-5):(k[1]*6),],LM[(k[2]*6-5):(k[2]*6),],k[3])
idr=rbind(idr,idt)
k=c(14,18,3)	#species a-b combination with sbn
idt=Maglio73(LM[(k[1]*6-5):(k[1]*6),],LM[(k[2]*6-5):(k[2]*6),],k[3])
idr=rbind(idr,idt)
k=c(16,18,3)	#species a-b combination with sbn
idt=Maglio73(LM[(k[1]*6-5):(k[1]*6),],LM[(k[2]*6-5):(k[2]*6),],k[3])
idr=rbind(idr,idt)
plot(idr[,4],idr[,6]);nrow(idr)
#--------------------------------------------------------------------------
idrx=idr[!rowSums(!is.finite(idr)),]  #remove rows that have -Inf
nrow(idr);nrow(idrx)
writefile<-"C://R_aaROEVchapt09//9.2.10_Maglio1973_Elephantidae_out.csv"
write.csv(idrx,file=writefile,na="NA",row.names=FALSE) 


