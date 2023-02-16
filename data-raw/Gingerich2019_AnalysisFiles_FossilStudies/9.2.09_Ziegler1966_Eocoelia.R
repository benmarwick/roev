##Philip D. Gingerich RATES OF EVOLUTION (Cambridge University Press 2019)
##9.2.09_Ziegler1966_Eocoelia
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
file1<-"../Gingerich2019_AnalysisFiles_FossilStudies/9.2.09_Ziegler1966_Eocoelia.csv"
Z<-read.csv(file1,header=TRUE,sep=",",quote="\"",dec=".",fill=TRUE)
Z;nrow(Z)
#------------------------------------------------------Logged matrix
LZ=Z[,c(1,4:12)];colnames(LZ)=c('age','Nrib','LMrib','LSrib','Nthet','Mthet',
  'Sthey','Nh/w','Mh/w','Sh/2');LZ[1,]
LZ[,3]=log(Z[,5]);LZ[,4]=Z[,6]/Z[,5]
LZ[1:3,]
#==========================================================================
idr=matrix(nrow=0,ncol=9)
  colnames(idr)=c('int.g','diff.sd','rate.sd.g','log.i','log|d|','log|r|',
    'sbn','wgt','li+lr')
#----------------------------------------intervals in line 1:paracone hgt
nr=13					#number of rows
nn=.5*(nr-1)*nr;3*nn
nc=0					#zero n count
for (c in seq(3,9,3)){		#columns to analyze
  for (i in 1:(nr-1)){		#increment
    for (sr in 1:(nr-i)){		#start row
      id<-numeric(length=9)
      nc=nc+1
      iyr=1000000*(LZ[(sr+i),1]-LZ[sr,1])		#interval in years
      gt=4 #10^(.266*log10(wg)-.553)		#Eq 3.2 generation time in years
      id[1]=iyr/gt					#interval in generations
      meandiff=LZ[(sr),c]-LZ[sr+i,c]		#mean diff.
      psd=PoolSD(LZ[sr,c-1],LZ[sr+i,c-1],LZ[sr,c+1],LZ[sr+i,c+1])#n1,n2,sd1,sd2
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
#---------------------------------------------------------------------
nrow(idr)
plot(idr[,4],idr[,6])
min(idr[,9]);max(idr[,9])
#--------------------------------------------------------------------------
idrx=idr[!rowSums(!is.finite(idr)),]  #remove rows that have -Inf
nrow(idr);nrow(idrx)
writefile<-"C://R_aaROEVchapt09//9.2.09_Ziegler1966_Eocoelia_out.csv"
write.csv(idrx,file=writefile,na="NA",row.names=FALSE)


