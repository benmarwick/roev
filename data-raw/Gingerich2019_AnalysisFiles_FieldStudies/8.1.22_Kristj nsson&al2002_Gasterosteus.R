##Philip D. Gingerich RATES OF EVOLUTION (Cambridge University Press 2019) 
##8.1.22_Kristjánsson&al2002_Gasterosteus
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
file1<-"C://R_aaROEVchapt08//8.1.22_Kristjánsson&al2002_Gasterosteus.csv" 
K<-read.csv(file1,header=TRUE,sep=",",quote="\"",dec=".",fill=TRUE)
K[1:4,]
#------------------------------------------------------Standardized traits
KS<-matrix(nrow=18,ncol=11)
  colnames(KS)=c('Age1','Age2','MarN','MarM','MarS','LavN','LavM','LavS',
    'MudN','MudM','MudS')
KS[1:18,1:11]=data.matrix(K)[1:18,2:12]
KS[1:3,]
gentime=2
idr=matrix(nrow=2*18,ncol=9)
  colnames(idr)=c('int','diff.sd','rate.sd.g','log.i','log|d|','log|r|',
    'sbn','wgt','sum')		#fgen is fraction of a generation
nc=0
for (i in 1:2){
  for (j in 1:18){
    k=(i-1)*18
    nc=nc+1
    #i=2;nc=1
    idr[nc,1]=(KS[nc-k,2]-KS[nc-k,1])/gentime			#gen
    meandiff=KS[nc-k,4+i*3]-KS[nc-k,4]
    poolsd=PoolSD(KS[nc-k,3],KS[nc-k,3+i*3],
      KS[nc-k,5],KS[nc-k,5+i*3])					#n1,n2,sd1,sd2
    idr[nc,2]=meandiff/poolsd						#diff.sd
    idr[nc,3]=abs(idr[nc,2])/idr[nc,1]				#rate.sd.gen
    idr[nc,4]=log10(idr[nc,1])					#log.i
    idr[nc,5]=log10(abs(idr[nc,2]))					#log.d
    idr[nc,6]=log10(idr[nc,3])					#log.r
    if(idr[nc,1]==1){idr[nc,7]=1}else{idr[nc,7]=2}		#sbn
    idr[nc,8]=1/idr[nc,1]						#wgt
    idr[nc,9]=idr[nc,4]+idr[nc,6]					#sum
    
  }
}
idr[1:3,]
nrow(idr)
idrxS=idr[!rowSums(!is.finite(idr)),]  #remove rows that have -Inf
nrow(idrxS)
#------------------------------------------------------Standardized traits
KT<-matrix(nrow=8,ncol=11)
  colnames(KT)=c('Age1','Age2','MarN','MarM','MarS','LavN','LavM','LavS',
    'MudN','MudM','MudS')
KT[1:8,1:11]=data.matrix(K)[19:26,2:12]
KT[1:3,]
gentime=2
idr=matrix(nrow=2*8,ncol=9)
  colnames(idr)=c('int','diff.sd','rate.sd.g','log.i','log|d|','log|r|',
    'sbn','wgt','sum')		#fgen is fraction of a generation
nc=0
for (i in 1:2){
  for (j in 1:8){
    k=(i-1)*8
    nc=nc+1
    #i=2;j=1;k=(i-1)*8;nc=j+(i-1)*8
    idr[nc,1]=(KT[nc-k,2]-KT[nc-k,1])/gentime			#gen
    meandiff=log(KT[nc-k,4+i*3])-log(KT[nc-k,4])
    poolsd=PoolSD(KT[nc-k,3],KT[nc-k,3+i*3],
      KT[nc-k,5]/KT[nc-k,4],KT[nc-k,5+i*3]/KT[nc-k,4+i*3])					#n1,n2,sd1,sd2
    idr[nc,2]=meandiff/poolsd						#diff.sd
    idr[nc,3]=abs(idr[nc,2])/idr[nc,1]				#rate.sd.gen
    idr[nc,4]=log10(idr[nc,1])					#log.i
    idr[nc,5]=log10(abs(idr[nc,2]))					#log.d
    idr[nc,6]=log10(idr[nc,3])					#log.r
    if(idr[nc,1]==1){idr[nc,7]=1}else{idr[nc,7]=2}		#sbn
    idr[nc,8]=1/idr[nc,1]						#wgt
    idr[nc,9]=idr[nc,4]+idr[nc,6]					#sum
    
  }
}
idr[1:3,]
nrow(idr)
idrxT=idr[!rowSums(!is.finite(idr)),]  #remove rows that have -Inf
nrow(idrxT)
#--------------------------------------------------------------------------
idrxA=rbind(idrxS,idrxT)
writefile<-"C://R_aaROEVchapt08//8.1.22_Kristjánsson&al2002_Gasterosteus_out.csv"
write.csv(idrxA,file=writefile,na="NA",row.names=FALSE) 
