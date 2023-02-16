##Philip D. Gingerich RATES OF EVOLUTION (Cambridge University Press 2019)
##9.2.06_Kurten1955_bears
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
file1<-"../Gingerich2019_AnalysisFiles_FossilStudies//9.2.06_Kurten1955_bears.csv"
K<-read.csv(file1,header=TRUE,sep=",",quote="\"",dec=".",fill=TRUE)
K;nrow(K)
#------------------------------------------------------Logged matrix
LK=K;colnames(LK)=c('species','source','N','LnM','LnS','V',
  'Age_rev');LK[1,]
LK[,4]=log(K[,4]);LK[,5]=K[,5]/K[,4]
K;LK
#==========================================================================
idr=matrix(nrow=0,ncol=9)
  colnames(idr)=c('int.g','diff.sd','rate.sd.g','log.i','log|d|','log|r|',
    'sbn','wgt','li+lr')
#----------------------------------------intervals in line 1:paracone hgt
col=4					#UM1L
nr=nrow(LK)				#number of rows
nn=.5*(nr-1)*nr;nn
id<-numeric(length=9)
wgtgen<-matrix(nrow=15,ncol=6)#error check
  colnames(wgtgen)=c('nc','m1','w1g','w2g','wg','gt')
nc=0					#n count
for (i in 1:(nr-1)){		#increment
  for (sr in 1:(nr-i)){		#start row
    nc=nc+1
    iyr=1000000*(LK[(sr),7]-LK[sr+i,7])	#interval in years
    w1g=10^(1.2652*log10(exp(LK[sr,4]))+3.4863)#ursid regr. Van Valkenburgh (1990)
    w2g=10^(1.2652*log10(exp(LK[sr+1,4]))+3.4863)
    wg=exp(.5*(log(w1g)+log(w2g)))		#mean wgt
    gt=10^(.266*log10(wg)-.553)		#Eq 3.2 generation time in years
    wgtgen[nc,1:6]=c(nc,exp(LK[sr,4]),w1g,w2g,wg,gt)	#error check
    id[1]=iyr/gt					#interval in generations
    meandiff=LK[(sr+i),col]-LK[sr,col]	#mean diff.
    psd=PoolSD(LK[sr,3],LK[sr+i,3],LK[sr,col+1],LK[sr+i,col+1])#n1,n2,sd1,sd2
    id[2]=meandiff/psd				#diff.sd
    id[3]=abs(id[2])/id[1]			#rate.sd.gen
    id[4]=log10(id[1])				#log.i
    id[5]=log10(abs(id[2]))			#log.d
    id[6]=log10(id[3])				#log.r
    if(i==1){id[7]=2}else{id[7]=3}		#sbn
    id[8]=1/id[1]					#wgt
    id[9]=id[4]+id[6]				#sum
    print(id)
    idr=rbind(idr,id)
  }
}
wgtgen	#error check
#---------------------------------------------------------------------
nrow(idr)
plot(idr[,4],idr[,6])
min(idr[,9]);max(idr[,9])
#--------------------------------------------------------------------------
idrx=idr[!rowSums(!is.finite(idr)),]  #remove rows that have -Inf
nrow(idr);nrow(idrx)
writefile<-"C://R_aaROEVchapt09//9.2.06_Kurten1955_bears_out.csv"
write.csv(idrx,file=writefile,na="NA",row.names=FALSE)


