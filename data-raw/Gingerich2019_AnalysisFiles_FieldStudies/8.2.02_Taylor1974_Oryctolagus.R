##Philip D. Gingerich RATES OF EVOLUTION (Cambridge University Press 2019) 
##8.2.02_Taylor1974_Oryctolagus
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
file1<-"C://R_aaROEVchapt08//8.2.02_Taylor1974_Oryctolagus.csv"
T<-read.csv(file1,header=TRUE,sep=",",quote="\"",dec=".",fill=TRUE)
T[1:4,]
#----------------------------------------------------------------------
T[1,]
LT<-matrix(nrow=nrow(T),ncol=6*3)		#Ln matrix for sclateri
  colnames(LT)=c('N1','M1','S1','N2','M2','S2','N3','M3','S3',
    'N4','M4','S4','N5','M5','S5','N6','M6','S6')
for (i in 1:6){
  for (j in 1:nrow(LT)){
    LT[j,3*i-2]=T[j,1+i*5-3]				#N
    LT[j,3*i-1]=log(T[j,1+i*5-2])			#ln mean
    LT[j,3*i]=.01*T[j,1+i*5]				#ln sd
  }
}
LT[1:3,]
#============================= Rate calc w. Werribee, Pine Plains, etc. roots
gentime=1.5
n=6
nn=.5*(n-1)*n;nn
idr=matrix(nrow=nn*nrow(LT),ncol=9)
  colnames(idr)=c('int','diff.sd','rate.sd.g','log.i','log|d|','log|r|',
    'sbn','wgt','LI+LR')
nc=0
for (i in 1:5){	
  if(i==1){dt=105}	#divergence time in years before start of study
  if(i==2){dt=95}
  if(i==3){dt=85}
  if(i==4){dt=85}
  if(i==5){dt=55}
  for (j in (i+1):6){
    for (tr in 1:nrow(LT)){	#trait number
      nc=nc+1
      idr[nc,1]=(2*dt)/gentime					#interval gen
      meandiff=LT[tr,3*i-1]-LT[tr,3*j-1]				#mean diff.
      poolsd=PoolSD(LT[tr,3*i-2],LT[tr,3*j-2],
        LT[tr,3*i],LT[tr,3*j])					#n1,n2,sd1,sd2
      idr[nc,2]=meandiff/poolsd					#diff.sd
      idr[nc,3]=abs(idr[nc,2])/idr[nc,1]				#rate.sd.gen
      idr[nc,4]=log10(idr[nc,1])					#log.i
      idr[nc,5]=log10(abs(idr[nc,2]))				#log.d
      idr[nc,6]=log10(idr[nc,3])					#log.r
      idr[nc,7]=2								#sbn
      idr[nc,8]=1/idr[nc,1]						#wgt
      idr[nc,9]=idr[nc,4]+idr[nc,6]
    }
  }
}
idr
idrx=idr[!rowSums(!is.finite(idr)),]  #remove rows that have -Inf
nrow(idr);nrow(idrx)
plot(idr[,4],idr[,6])
writefile<-"C://R_aaROEVchapt08//8.2.02_Taylor1974_Oryctolagus_out.csv"
write.csv(idrx,file=writefile,na="NA",row.names=FALSE) 


