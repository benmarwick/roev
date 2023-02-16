##Philip D. Gingerich RATES OF EVOLUTION (Cambridge University Press 2019) 
##8.1.33_Caruso&al2014_Plethodon
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
file1<-"C://R_aaROEVchapt08//8.1.33_Caruso&al2014_Plethodon.csv" 
C<-read.csv(file1,header=TRUE,sep=",",quote="\"",dec=".",fill=TRUE)
C[1:4,]
ncol(C)
#------------------------------------------------------Logged matrix
LC<-matrix(nrow=nrow(C),ncol=6)
  colnames(LC)=c('sp.','yr','N','Ln M','Ln S','gen');LC[1:3,]
for (rn in 1:nrow(LC)){							#column no.
  LC[rn,1]=C[rn,1]
  if(C[rn,3]=='1950s'){LC[rn,2]=1955}
  if(C[rn,3]=='1960s'){LC[rn,2]=1965}
  if(C[rn,3]=='1970s'){LC[rn,2]=1975}
  if(C[rn,3]=='1980s'){LC[rn,2]=1985}
  if(C[rn,3]=='1990s'){LC[rn,2]=1995}
  if(C[rn,3]=='2000s'){LC[rn,2]=2005}
  if(C[rn,3]=='2010s'){LC[rn,2]=2015}
  LC[rn,3]=C[rn,4]
  LC[rn,4]=log(C[rn,5])
  LC[rn,5]=sqrt(C[rn,6])/C[rn,5]
  LC[rn,6]=exp(.77*log(C[rn,5])-1.14)
};LC[1:4,]
#================================================== Calculate rates 
idr=matrix(nrow=0,ncol=9)
  colnames(idr)=c('int','diff.sd','rate.sd.g','log.i','log|d|','log|r|',
    'sbn','wgt','sumIR')
#cc=0
for(sn in 2:15){
  test<-matrix(nrow=0,ncol=6)
    colnames(test)=c('sp.','yr','N','Ln SVL','s','gen')
  for(i in 1:nrow(LC)){
    if(LC[i,1]==sn){
      test=rbind(test,LC[i,])
    }
  }
  for (sr in 1:(nrow(test)-1)){ 					#start row
    for (i in 1:(nrow(test)-sr)){					#interval steps
      print(c(sr,i))
      id<-numeric(length=9)
      k=test[sr+i,2]-test[sr,2]					#interval in yrs
      gt=exp(.77*(.5*(test[sr,4]+test[sr+i,4]))-1.14)
      if(k<=gt){id[1]=1}else{id[1]=k/gt}				#interval in gen
      meandiff=test[(sr+i),4]-test[sr,4]				#mean diff.
      psd=PoolSD(test[sr,3],test[sr+i,3],test[sr,5],test[sr+i,5])#n1,n2,sd1,sd2
      id[2]=meandiff/psd							#diff.sd
      id[3]=abs(id[2])/id[1]						#rate.sd.gen
      id[4]=log10(id[1])							#log.i
      id[5]=log10(abs(id[2]))						#log.d
      id[6]=log10(id[3])							#log.r
      if(i==1){
        if(id[1]<=1){id[7]=1}else{id[7]=2}				#sbn
      }else{id[7]=3}
      id[8]=1/id[1]							#wgt
      id[9]=id[4]+id[6]							#sum
      print(id)
      idr=rbind(idr,id)
    }
  }
}
nrow(idr)
plot(idr[,4],idr[,6])
summary(idr[,9])

#--------------------------------------------------------------------------
idrx=idr[!rowSums(!is.finite(idr)),]  #remove rows that have -Inf
idrxC=idrx[,1:8]
idrc=idrx[idrx[,9]>=1,]  #keep only rows where interval >= gentime
idrcC=idrc[,1:8]
nrow(idr);nrow(idrx);nrow(idrxC);nrow(idrcC)
writefile<-"C://R_aaROEVchapt08//8.1.33_Caruso&al2014_Plethodon_out.csv"
write.csv(idrxC,file=writefile,na="NA",row.names=FALSE) 


