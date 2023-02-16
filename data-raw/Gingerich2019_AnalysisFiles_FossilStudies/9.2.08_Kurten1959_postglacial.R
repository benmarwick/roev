##Philip D. Gingerich RATES OF EVOLUTION (Cambridge University Press 2019)
##9.2.08_Kurten1959_postglacial
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
file1<-"../Gingerich2019_AnalysisFiles_FossilStudies/9.2.08_Kurten1959_postglacial.csv"
K<-read.csv(file1,header=TRUE,sep=",",quote="\"",dec=".",fill=TRUE)
K[1:5,];nrow(K)
#------------------------------------------------------Logged matrix
LK=K[,c(1,4:6)]
  colnames(LK)=c('species','gentime','age','LnM');LK[1,]
LK[,4]=log(K[,6])
LK
#==========================================================================
idr=matrix(nrow=0,ncol=9)
  colnames(idr)=c('int.g','diff.sd','rate.sd.g','log.i','log|d|','log|r|',
    'sbn','wgt','li+lr')
#------------------------------------------------------------ studies
idt=Kurten59(LK[1:4,]);idr=rbind(idr,idt)
idt=Kurten59(LK[5:6,]);idr=rbind(idr,idt)
idt=Kurten59(LK[7:8,]);idr=rbind(idr,idt)
idt=Kurten59(LK[9:11,]);idr=rbind(idr,idt)
idt=Kurten59(LK[12:13,]);idr=rbind(idr,idt)
idt=Kurten59(LK[14:17,]);idr=rbind(idr,idt)
idt=Kurten59(LK[18:19,]);idr=rbind(idr,idt)
idt=Kurten59(LK[20:23,]);idr=rbind(idr,idt)
idt=Kurten59(LK[24:25,]);idr=rbind(idr,idt)
idt=Kurten59(LK[26:28,]);idr=rbind(idr,idt)
idt=Kurten59(LK[29:31,]);idr=rbind(idr,idt)
idt=Kurten59(LK[32:33,]);idr=rbind(idr,idt)

#---------------------------------------------------------------------
nrow(idr)
plot(idr[,4],idr[,6])
min(idr[,9]);max(idr[,9])
#--------------------------------------------------------------------------
idrx=idr[!rowSums(!is.finite(idr)),]  #remove rows that have -Inf
nrow(idr);nrow(idrx)
writefile<-"C://R_aaROEVchapt09//9.2.08_Kurten1959_postglacial_out.csv"
write.csv(idrx,file=writefile,na="NA",row.names=FALSE)


