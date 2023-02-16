##Philip D. Gingerich RATES OF EVOLUTION (Cambridge University Press 2019) 
##9.2.20_MacFadden1988_Equidae
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
file1<-"C://R_aaROEVchapt09//9.2.20_MacFadden1988_Equidae.csv" 
M<-read.csv(file1,header=TRUE,sep=",",quote="\"",dec=".",fill=TRUE)
nrow(M);M[1:5,]
#------------------------------------------------------Logged matrix
LM<-matrix(nrow=nrow(M),ncol=8)
  colnames(LM)=c('comp','Mwgt','Imy','Dlnwsd','Dln1sd','Dln2sd','Dln3sd','Dln4sd')
for (i in 1:nrow(M)){
  LM[i,1]=M[i,1]							#comparison
  LM[i,2]=exp(.5*(log(M[i,3])+log(M[i,5])))		#geom. mean wgt kg
  LM[i,3]=M[i,6]							#int_yr
  LM[i,4]=(log(M[i,3])-log(M[i,5]))/(.5*(.103+.176))	#sd from MacFadden 1986
  LM[i,5]=M[i,6]*M[i,7]/.05					#sd for linear meas.
  LM[i,6]=M[i,6]*M[i,8]/.05					#sd for linear meas.
  LM[i,7]=M[i,6]*M[i,9]/.05					#sd for linear meas.
  LM[i,8]=M[i,6]*M[i,10]/.05					#sd for linear meas.
}
LM

#=============================================================== Rates
idr=matrix(nrow=0,ncol=9)
  colnames(idr)=c('int.g','diff.sd','rate.sd.g','log.i','log|d|','log|r|',
    'sbn','wgt','li+lr')
#------------------------------------------------------ 
LM[1:5,]
nr=nrow(LM);nr						#number of rows
zc=0							#zero count
nc=0
for (cn in seq(4,8,1)){				#columns of means
  for (i in 1:nr){				#row
      id<-numeric(length=9)
      nc=nc+1
      iyr=1000000*LM[i,3]	#interval in years
      gt=10^(.266*log10(1000*LM[i,2])-.553)	#Eq 3.2 generation time in years
      id[1]=iyr/gt					#interval in generations
      meandiff=LM[i,cn]		#mean diff.
      if (!is.na(meandiff)){if (meandiff==0){zc=zc+1}}
      id[2]=LM[i,cn]					#diff.sd
      id[3]=abs(id[2])/id[1]				#rate.sd.gen
      id[4]=log10(id[1])				#log.i
      id[5]=log10(abs(id[2]))				#log.d
      id[6]=log10(id[3])				#log.r
      id[7]=2
      id[8]=1/id[1]					#wgt
      id[9]=id[4]+id[6]					#sum
    	idr=rbind(idr,id)
  }
}
#--------------------------------------------------------------------------
nrow(idr)
#plot(idr[,4],idr[,6])
#--------------------------------------------------------------------------
idrx=idr[!rowSums(!is.finite(idr)),]  #remove rows that have -Inf
nrow(idr);nrow(idrx);zc
min(idrx[,9]);max(idrx[,9])

#--------------------------------------------------------------------------
writefile<-"C://R_aaROEVchapt09//9.2.20_MacFadden1988_Equidae_out.csv"
write.csv(idrx,file=writefile,na="NA",row.names=FALSE) 



