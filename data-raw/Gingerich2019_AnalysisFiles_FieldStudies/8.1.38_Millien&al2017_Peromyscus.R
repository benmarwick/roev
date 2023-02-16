##Philip D. Gingerich RATES OF EVOLUTION (Cambridge University Press 2019) 
##8.1.38_Millien&al2017_Peromyscus
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
#============================================================ Load files
file1<-"C://R_aaROEVchapt08//8.1.38_Millien&al2017_Peromyscus.csv"
M<-read.csv(file1,header=TRUE,sep=",",quote="\"",dec=".",fill=TRUE)
M
#----------------------------------------------------set up matrices
idr=matrix(nrow=16,ncol=7)
  colnames(idr)=c('int','diff.sd','rate.sd.g','log.i','log|d|','log|r|','wgt')
gen=90
#----------------------------------------- leucopus females
id<-matrix(nrow=4,ncol=10)
  rownames(id)=c('LFbm','LFbl','LFa1','LFa2')
  colnames(id)=c('N1','M1','S1','N2','M2','S2',
    'diff','pool.sd','diff.sd','rate.sd.g');id
M[1:3,7:10]
hs=1;he=11;rs=12;re=27
LF<-matrix(nrow=2,ncol=12)
  rownames(LF)=c('Hist.','Recent')
  colnames(LF)=c('N1','Mean1','SD1','N2','Mean2','SD2',
    'N3','Mean3','SD3','N4','Mean4','SD4')
c<-c(1,4,7,10)	#column to fill
for (i in 1:2){
  LF[1,c[i]]=sum(!is.na(M[hs:he,(i+6)]))
  LF[1,c[i]+1]=mean(log(M[hs:he,(i+6)]),na.rm=TRUE)
  LF[1,c[i]+2]=sd(log(M[hs:he,(i+6)]),na.rm=TRUE)
  LF[2,c[i]]=sum(!is.na(M[rs:re,(i+6)]))
  LF[2,c[i]+1]=mean(log(M[rs:re,(i+6)]),na.rm=TRUE)
  LF[2,c[i]+2]=sd(log(M[rs:re,(i+6)]),na.rm=TRUE)
}
for (i in 3:4){
  LF[1,c[i]]=sum(!is.na(M[hs:he,(i+6)]))
  LF[1,c[i]+1]=mean(M[hs:he,(i+6)],na.rm=TRUE)
  LF[1,c[i]+2]=sd(M[hs:he,(i+6)],na.rm=TRUE)
  LF[2,c[i]]=sum(!is.na(M[rs:re,(i+6)]))
  LF[2,c[i]+1]=mean(M[rs:re,(i+6)],na.rm=TRUE)
  LF[2,c[i]+2]=sd(M[rs:re,(i+6)],na.rm=TRUE)
};LF
#----------------------------------------------------calculate rates
for (i in 1:4){
  id[i,1]=LF[1,i*3-2];id[i,2]=LF[1,i*3-1];id[i,3]=LF[1,i*3]
  id[i,4]=LF[2,i*3-2];id[i,5]=LF[2,i*3-1];id[i,6]=LF[2,i*3]
  id[i,7]=id[i,5]-id[i,2]
  id[i,8]=sqrt((id[i,1]*id[i,3]^2+id[i,4]*id[i,6]^2)/(id[i,1]+id[i,4]))
  id[i,9]=id[i,7]/id[i,8];id[i,10]=abs(id[i,9])/gen
}
for (i in 1:4){
  idr[i,1]=gen;idr[i,2]=id[i,9];idr[i,3]=id[i,10]
  idr[i,4]=log10(idr[i,1]);idr[i,5]=log10(abs(idr[i,2]))
  idr[i,6]=log10(idr[i,3]);idr[i,7]=1/idr[i,1]
};idr
#----------------------------------------- leucopus males
id<-matrix(nrow=4,ncol=10)
  rownames(id)=c('LMbm','LMbl','LMa1','LMa2')
  colnames(id)=c('N1','M1','S1','N2','M2','S2',
    'diff','pool.sd','diff.sd','rate.sd.g');id
M[1:3,7:10]
hs=28;he=42;rs=43;re=62
LM<-matrix(nrow=2,ncol=12)
  rownames(LM)=c('Hist.','Recent')
  colnames(LM)=c('N1','Mean1','SD1','N2','Mean2','SD2',
    'N3','Mean3','SD3','N4','Mean4','SD4')
c<-c(1,4,7,10)	#column to fill
for (i in 1:2){
  LM[1,c[i]]=sum(!is.na(M[hs:he,(i+6)]))
  LM[1,c[i]+1]=mean(log(M[hs:he,(i+6)]),na.rm=TRUE)
  LM[1,c[i]+2]=sd(log(M[hs:he,(i+6)]),na.rm=TRUE)
  LM[2,c[i]]=sum(!is.na(M[rs:re,(i+6)]))
  LM[2,c[i]+1]=mean(log(M[rs:re,(i+6)]),na.rm=TRUE)
  LM[2,c[i]+2]=sd(log(M[rs:re,(i+6)]),na.rm=TRUE)
}
for (i in 3:4){
  LM[1,c[i]]=sum(!is.na(M[hs:he,(i+6)]))
  LM[1,c[i]+1]=mean(M[hs:he,(i+6)],na.rm=TRUE)
  LM[1,c[i]+2]=sd(M[hs:he,(i+6)],na.rm=TRUE)
  LM[2,c[i]]=sum(!is.na(M[rs:re,(i+6)]))
  LM[2,c[i]+1]=mean(M[rs:re,(i+6)],na.rm=TRUE)
  LM[2,c[i]+2]=sd(M[rs:re,(i+6)],na.rm=TRUE)
};LM
#----------------------------------------------------calculate rates
for (i in 1:4){
  id[i,1]=LM[1,i*3-2];id[i,2]=LM[1,i*3-1];id[i,3]=LM[1,i*3]
  id[i,4]=LM[2,i*3-2];id[i,5]=LM[2,i*3-1];id[i,6]=LM[2,i*3]
  id[i,7]=id[i,5]-id[i,2]
  id[i,8]=sqrt((id[i,1]*id[i,3]^2+id[i,4]*id[i,6]^2)/(id[i,1]+id[i,4]))
  id[i,9]=id[i,7]/id[i,8];id[i,10]=abs(id[i,9])/gen
}
for (i in 1:4){
  idr[i+4,1]=gen;idr[i+4,2]=id[i,9];idr[i+4,3]=id[i,10]
  idr[i+4,4]=log10(idr[i+4,1]);idr[i+4,5]=log10(abs(idr[i+4,2]))
  idr[i+4,6]=log10(idr[i+4,3]);idr[i+4,7]=1/idr[i+4,1]
};idr
#----------------------------------------- maniculatus females
id<-matrix(nrow=4,ncol=10)
  rownames(id)=c('MFbm','MFbl','MFa1','MFa2')
  colnames(id)=c('N1','M1','S1','N2','M2','S2',
    'diff','pool.sd','diff.sd','rate.sd.g');id
M[1:3,7:10]
hs=63;he=73;rs=74;re=85
MF<-matrix(nrow=2,ncol=12)
  rownames(MF)=c('Hist.','Recent')
  colnames(MF)=c('N1','Mean1','SD1','N2','Mean2','SD2',
    'N3','Mean3','SD3','N4','Mean4','SD4')
c<-c(1,4,7,10)	#column to fill
for (i in 1:2){
  MF[1,c[i]]=sum(!is.na(M[hs:he,(i+6)]))
  MF[1,c[i]+1]=mean(log(M[hs:he,(i+6)]),na.rm=TRUE)
  MF[1,c[i]+2]=sd(log(M[hs:he,(i+6)]),na.rm=TRUE)
  MF[2,c[i]]=sum(!is.na(M[rs:re,(i+6)]))
  MF[2,c[i]+1]=mean(log(M[rs:re,(i+6)]),na.rm=TRUE)
  MF[2,c[i]+2]=sd(log(M[rs:re,(i+6)]),na.rm=TRUE)
}
for (i in 3:4){
  MF[1,c[i]]=sum(!is.na(M[hs:he,(i+6)]))
  MF[1,c[i]+1]=mean(M[hs:he,(i+6)],na.rm=TRUE)
  MF[1,c[i]+2]=sd(M[hs:he,(i+6)],na.rm=TRUE)
  MF[2,c[i]]=sum(!is.na(M[rs:re,(i+6)]))
  MF[2,c[i]+1]=mean(M[rs:re,(i+6)],na.rm=TRUE)
  MF[2,c[i]+2]=sd(M[rs:re,(i+6)],na.rm=TRUE)
};MF
#----------------------------------------------------calculate rates
for (i in 1:4){
  id[i,1]=MF[1,i*3-2];id[i,2]=MF[1,i*3-1];id[i,3]=MF[1,i*3]
  id[i,4]=MF[2,i*3-2];id[i,5]=MF[2,i*3-1];id[i,6]=MF[2,i*3]
  id[i,7]=id[i,5]-id[i,2]
  id[i,8]=sqrt((id[i,1]*id[i,3]^2+id[i,4]*id[i,6]^2)/(id[i,1]+id[i,4]))
  id[i,9]=id[i,7]/id[i,8];id[i,10]=abs(id[i,9])/gen
}
for (i in 1:4){
  idr[i+8,1]=gen;idr[i+8,2]=id[i,9];idr[i+8,3]=id[i,10]
  idr[i+8,4]=log10(idr[i+8,1]);idr[i+8,5]=log10(abs(idr[i+8,2]))
  idr[i+8,6]=log10(idr[i+8,3]);idr[i+8,7]=1/idr[i+8,1]
};idr
#----------------------------------------- maniculatus males
id<-matrix(nrow=4,ncol=10)
  rownames(id)=c('MMbm','MMbl','MMa1','MMa2')
  colnames(id)=c('N1','M1','S1','N2','M2','S2',
    'diff','pool.sd','diff.sd','rate.sd.g');id
M[1:3,7:10]
hs=86;he=101;rs=102;re=106
MM<-matrix(nrow=2,ncol=12)
  rownames(MM)=c('Hist.','Recent')
  colnames(MM)=c('N1','Mean1','SD1','N2','Mean2','SD2',
    'N3','Mean3','SD3','N4','Mean4','SD4')
c<-c(1,4,7,10)	#column to fill
for (i in 1:2){
  MM[1,c[i]]=sum(!is.na(M[hs:he,(i+6)]))
  MM[1,c[i]+1]=mean(log(M[hs:he,(i+6)]),na.rm=TRUE)
  MM[1,c[i]+2]=sd(log(M[hs:he,(i+6)]),na.rm=TRUE)
  MM[2,c[i]]=sum(!is.na(M[rs:re,(i+6)]))
  MM[2,c[i]+1]=mean(log(M[rs:re,(i+6)]),na.rm=TRUE)
  MM[2,c[i]+2]=sd(log(M[rs:re,(i+6)]),na.rm=TRUE)
}
for (i in 3:4){
  MM[1,c[i]]=sum(!is.na(M[hs:he,(i+6)]))
  MM[1,c[i]+1]=mean(M[hs:he,(i+6)],na.rm=TRUE)
  MM[1,c[i]+2]=sd(M[hs:he,(i+6)],na.rm=TRUE)
  MM[2,c[i]]=sum(!is.na(M[rs:re,(i+6)]))
  MM[2,c[i]+1]=mean(M[rs:re,(i+6)],na.rm=TRUE)
  MM[2,c[i]+2]=sd(M[rs:re,(i+6)],na.rm=TRUE)
};MM
#----------------------------------------------------calculate rates
for (i in 1:4){
  id[i,1]=MM[1,i*3-2];id[i,2]=MM[1,i*3-1];id[i,3]=MM[1,i*3]
  id[i,4]=MM[2,i*3-2];id[i,5]=MM[2,i*3-1];id[i,6]=MM[2,i*3]
  id[i,7]=id[i,5]-id[i,2]
  id[i,8]=sqrt((id[i,1]*id[i,3]^2+id[i,4]*id[i,6]^2)/(id[i,1]+id[i,4]))
  id[i,9]=id[i,7]/id[i,8];id[i,10]=abs(id[i,9])/gen
}
for (i in 1:4){
  idr[i+12,1]=gen;idr[i+12,2]=id[i,9];idr[i+12,3]=id[i,10]
  idr[i+12,4]=log10(idr[i+12,1]);idr[i+12,5]=log10(abs(idr[i+12,2]))
  idr[i+12,6]=log10(idr[i+12,3]);idr[i+12,7]=1/idr[i+12,1]
};idr

#================================================Bind in one idr matrix
idrx=idr[!rowSums(!is.finite(idr)),]  #remove rows that have -Inf
writefile<-"C://R_aaROEVchapt08//8.1.38_Millien&al2017_Peromyscus_out.csv"
write.csv(idrx,file=writefile,na="NA",row.names=FALSE) 

