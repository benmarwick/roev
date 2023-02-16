##Philip D. Gingerich RATES OF EVOLUTION (Cambridge University Press 2019)
##9.2.07_Kurten1959_Ursidae
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
file1<-"../Gingerich2019_AnalysisFiles_FossilStudies//9.2.07_Kurten1959_Ursidae.csv"
K<-read.csv(file1,header=TRUE,sep=",",quote="\"",dec=".",fill=TRUE)
K;nrow(K)
#------------------------------------------------------Logged matrix
LK=K[,c(1,3:7)]
  colnames(LK)=c('species','locality','N','LnM','LnS','Age_rev');LK[1,]
LK[,4]=log(K[,5]);LK[,5]=.05
K;LK
#==========================================================================
idr=matrix(nrow=0,ncol=9)
  colnames(idr)=c('int.g','diff.sd','rate.sd.g','log.i','log|d|','log|r|',
    'sbn','wgt','li+lr')
#----------------------------------------intervals in line 1:paracone hgt
gentime=13.5 #from COSEWIC (2012)
col=4					#UM1L
nr=nrow(LK)				#number of rows
nn=.5*(nr-1)*nr;nn
id<-numeric(length=9)
nc=0					#n count
for (i in 1:(nr-1)){		#increment
  for (sr in 1:(nr-i)){		#start row
    id<-numeric(length=9)
    nc=nc+1
    iyr=1000*(LK[(sr),6]-LK[sr+i,6])	#interval in years
    id[1]=iyr/gentime				#interval in generations
    meandiff=LK[(sr+i),col]-LK[sr,col]	#mean diff.
    id[2]=meandiff/0.05				#diff.sd
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
#---------------------------------------------------------------------
nrow(idr)
plot(idr[,4],idr[,6])
min(idr[,9]);max(idr[,9])
#--------------------------------------------------------------------------
idrx=idr[!rowSums(!is.finite(idr)),]  #remove rows that have -Inf
nrow(idr);nrow(idrx)
writefile<-"C://R_aaROEVchapt09//9.2.07_Kurten1959_Ursidae_out.csv"
write.csv(idrx,file=writefile,na="NA",row.names=FALSE)


