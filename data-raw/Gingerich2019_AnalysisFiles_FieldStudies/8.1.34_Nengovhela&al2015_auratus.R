##Philip D. Gingerich RATES OF EVOLUTION (Cambridge University Press 2019) 
##8.1.34_Nengovhela&al2015_auratus
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
#============================================================ Load file
file1<-"C://R_aaROEVchapt08//8.1.34_Nengovhela&al2015_auratus.csv"
N<-read.csv(file1,header=TRUE,sep=",",quote="\"",dec=".",fill=TRUE);N[1:3,]
LN=N;LN[,3]=log(N[,3]);LN[1:3,]
#--------------------------------------------------Multiple regression
LNmr=lm(LN$GLS~LN$TWCLS+LN$Latitude+LN$Longitude,data=LN)	#multiple regress.
summary(LNmr)
anova(LNmr);sqrt(anova(LNmr)[4,3])
LN=cbind(LN,LNmr[2]);LN[1:3,]				#cbind of residuals to LN
nrow(LN)
plot(LN[,c(2,7)],pch=19,col=LN[,4]-1)
LNr=lm(LN[,7]~LN[,2])			#regression of residuals on year
summary(LNr)
b=coef(summary(LNr))[1,1];b
m=coef(summary(LNr))[2,1];m
#--------------------------------------------------LN statistics
LN[1:3,]
LNst<-matrix(nrow=5,ncol=4)
  colnames(LNst)=c('yr','N','mean','stdev')
table(LN[,2],LN[,2])
for (i in 1:5){
  if (i==1){r=c(1,26)}					#1907 to 1931
  if (i==2){r=c(27,45)}					#1950 to 1963
  if (i==3){r=c(46,104)}				#1970 to 1979
  if (i==4){r=c(105,125)}				#1988 to 1994
  if (i==5){r=c(126,131)}				#2013 to 2013
  LNst[i,1]=median(LN[r[1]:r[2],2])
  LNst[i,2]=length(LN[r[1]:r[2],3])
  LNst[i,3]=mean(LN[r[1]:r[2],3])
  LNst[i,4]=sd(LN[r[1]:r[2],3])
};LNst
#============================================= Rate calculation
gentime=1								#generation time
LNst[1:5,]
n=nrow(LNst);n
nn=.5*(n-1)*n;nn
idr=matrix(nrow=nn,ncol=9)
  colnames(idr)=c('int','diff.sd','rate.sd.g','log.i','log|d|','log|r|',
    'sbn','wgt','sum')		#fgen is fraction of a generation
nc=0
for (k in 1:(n-1)){		#run length
  for (i in 1:(n-k)){  		#starting position
    nc=nc+1
    if(k<gentime){
      idr[nc,1]=1
    }else{
      idr[nc,1]=(LNst[(i+k),1]-LNst[i,1])/gentime
    }
    meandiff=LNst[(i+k),3]-LNst[i,3]					#mean diff.
    poolsd=PoolSD(LNst[i+k,2],LNst[i,2],LNst[i+k,4],LNst[i,4])	#n1,n2,sd1,sd2
    idr[nc,2]=meandiff/poolsd						#diff.sd
    idr[nc,3]=abs(idr[nc,2])/idr[nc,1]				#rate.sd.gen
    idr[nc,4]=log10(idr[nc,1])					#log.i
    idr[nc,5]=log10(abs(idr[nc,2]))					#log.d
    idr[nc,6]=log10(idr[nc,3])					#log.r
    if(k==1){
      idr[nc,7]=2								#sbn
    }else{
      idr[nc,7]=3
    }
    idr[nc,8]=1/idr[nc,1]						#wgt
    idr[nc,9]=idr[nc,4]+idr[nc,6]					#sum
  }
}
#------------------------------- add idr line from overall regression
lines(c(min(LN[,2]),max(LN[,2])),c(m*min(LN[,2])+b,m*max(LN[,2])+b),
  lty=1,lwd=2,col=2)
idra=numeric(length=9)
idra[1]=(max(LN[,2])-min(LN[,2]))/gentime			#diff.g
  meandiff=(m*max(LN[,2])+b)-(m*min(LN[,2])+b)
  poolsd=sqrt(anova(LNmr)[4,3])				#from ancova
idra[2]=meandiff/poolsd						#diff.sd
idra[3]=abs(idra[2])/idra[1]					#rate.sd.g
idra[4]=log10(idra[1])						#log.i
idra[5]=log10(abs(idra[2]))					#log.d
idra[6]=log10(idra[3])						#log.r
idra[7]=3								#sbn
idra[8]=1/idra[1]							#wgt
idra[9]=idra[4]+idra[6]						#sum
idr=rbind(idr,idra)
idr[1:3,]
idrx=idr[!rowSums(!is.finite(idr)),]  #remove rows that have -Inf
nrow(idr);nrow(idrx)
writefile<-"C://R_aaROEVchapt08//8.1.34_Nengovhela&al2015_auratus_out.csv"
write.csv(idrx,file=writefile,na="NA",row.names=FALSE) 
warnings()


