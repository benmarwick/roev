##Philip D. Gingerich RATES OF EVOLUTION (Cambridge University Press 2019) 
##8.1.08_Diamond&al1989_Myzomelasclateri
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
file1<-"C://R_aaROEVchapt08//8.1.08_Diamond&al1989_Myzomelasclateri.csv"
D<-read.csv(file1,header=TRUE,sep=",",quote="\"",dec=".",fill=TRUE)
D
#----------------------------------------------------------------------
D[1,]
LS<-matrix(nrow=2*6,ncol=3*6)		#Ln matrix for sclateri
  colnames(LS)=c('N1','M1','S1','N2','M2','S2','N3','M3','S3',
    'N4','M4','S4','N5','M5','S5','N6','M6','S6')
for (i in 1:nrow(LS)){
  for (j in 1:(ncol(LS)/3)){
    LS[i,3*j-2]=D[i,j*3+3]				#N
    LS[i,3*j-1]=log(D[i,j*3+1])			#ln mean
    LS[i,3*j]=D[i,j*3+2]/D[i,j*3+1]			#ln sd
  }
}
nrow(LS)
LS[seq(3,11,2),1:18]
#====================================== Average non-Long sites
AS<-matrix(nrow=2,ncol=3*6)	#Ln matrix for aveaged non-Long sclateri
  colnames(LS)=c('N1','M1','S1','N2','M2','S2','N3','M3','S3',
    'N4','M4','S4','N5','M5','S5','N6','M6','S6')
Nnz<-numeric(length=6)
for (tr in 1:6){  #6){	#trait number for wing, tail, bill, culmen, tarsus, wgt
  Nsum=sum(LS[seq(3,11,2),3*tr-2],na.rm=TRUE)		#male line
  Nnz[tr]=sum(LS[seq(3,11,2),3*tr-2]!=0,na.rm=TRUE)	#non-zero entries
  AS[1,3*tr-2]=Nsum
  AS[1,3*tr-1]=sum(LS[seq(3,11,2),3*tr-2]*LS[seq(3,11,2),3*tr-1],na.rm=TRUE)/
    Nsum
  AS[1,3*tr]=sqrt(sum((LS[seq(3,11,2),3*tr-2]-1)*LS[seq(3,11,2),3*tr]^2,
    na.rm=TRUE)/(Nsum-Nnz[tr]))

  Nsum=sum(LS[seq(4,12,2),3*tr-2],na.rm=TRUE)		#female line
  Nnz[tr]=sum(LS[seq(4,12,2),3*tr-2]!=0,na.rm=TRUE)	#non-zero entries
  AS[2,3*tr-2]=Nsum
  AS[2,3*tr-1]=sum(LS[seq(4,12,2),3*tr-2]*LS[seq(4,12,2),3*tr-1],na.rm=TRUE)/
    Nsum
  AS[2,3*tr]=sqrt(sum((LS[seq(4,12,2),3*tr-2]-1)*LS[seq(4,12,2),3*tr]^2,
    na.rm=TRUE)/(Nsum-Nnz[tr]))
}
AS
#====================================== Rate calc male & female wing length
gentime=3
idr=matrix(nrow=6*2,ncol=8)
  colnames(idr)=c('int','diff.sd','rate.sd.g','log.i','log|d|','log|r|',
    'sbn','wgt')		#fgen is fraction of a generation
nc=0
for (tr in 1:6){	#trait number for wing, tail, bill, culmen, tarsus, wgt
  nc=nc+1
  idr[nc,1]=100								#interval gen
  meandiff=LS[1,3*tr-1]-AS[1,3*tr-1]				#mean diff.
  poolsd=PoolSD(LS[1,3*tr-2],AS[1,3*tr-2],
    LS[1,3*tr],AS[1,3*tr])						#n1,n2,sd1,sd2
  idr[nc,2]=meandiff/poolsd						#diff.sd
  idr[nc,3]=abs(idr[nc,2])/idr[nc,1]				#rate.sd.gen
  idr[nc,4]=log10(idr[nc,1])					#log.i
  idr[nc,5]=log10(abs(idr[nc,2]))					#log.d
  idr[nc,6]=log10(idr[nc,3])					#log.r
  idr[nc,7]=2								#sbn
  idr[nc,8]=1/idr[nc,1]						#wgt
  nc=nc+1
  idr[nc,1]=100								#interval gen
  meandiff=LS[2,3*tr-1]-AS[2,3*tr-1]				#mean diff.
  poolsd=PoolSD(LS[2,3*tr-2],AS[2,3*tr-2],
    LS[2,3*tr],AS[2,3*tr])						#n1,n2,sd1,sd2
  idr[nc,2]=meandiff/poolsd						#diff.sd
  idr[nc,3]=abs(idr[nc,2])/idr[nc,1]				#rate.sd.gen
  idr[nc,4]=log10(idr[nc,1])					#log.i
  idr[nc,5]=log10(abs(idr[nc,2]))					#log.d
  idr[nc,6]=log10(idr[nc,3])					#log.r
  idr[nc,7]=2								#sbn
  idr[nc,8]=1/idr[nc,1]						#wgt
};idr
writefile<-"C://R_aaROEVchapt08//8.1.08_Diamond&al1989_Myzomelasclateri_out.csv"
write.csv(idr,file=writefile,na="NA",row.names=FALSE) 


