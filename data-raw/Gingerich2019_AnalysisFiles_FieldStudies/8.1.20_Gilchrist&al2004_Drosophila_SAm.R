##Philip D. Gingerich RATES OF EVOLUTION (Cambridge University Press 2019) 
##8.1.20_Gilchrist&al2004_Drosophila_SAm
#==========================================================================
cat(rep("\n",50))							#clear console
#Philip D. Gingerich May 12, 2016 RATES OF EVOLUTION
print (date()) #Ctr-s to save, Ctr-a to select all, Ctr-r to run
#dim() for dimensions, summary() for summary, dir() for directory
rm(list=ls(all=TRUE))#remove/clear all prev. variables
assign('last.warning',NULL,envir=baseenv())
ptm<-proc.time()
##======================================================================##
##library("devtools");library(roxygen2)
##setwd("c:/R_aaPackages/RATES");document();setwd("..");install("RATES") 
#=====================================================================Setup
library(RATES);library(MASS)
#--------------------------------------------------------------------------
file1<-"C://R_aaROEVchapt08//8.1.20_Gilchrist&al2004_Drosophila_SAm.csv"
G<-read.csv(file1,header=TRUE,sep=",",quote="\"",dec=".",fill=TRUE)
G
#------------------------------------------------------Logged matrix
LG<-matrix(nrow=nrow(G),ncol=14)
  rownames(LG)=c('DA','EU','ME','SA','CE','BE','DA','EU','ME','SA','CE','BE')
  colnames(LG)=c('year','N','Ln_L1mf','Ln_L1sf','Ln_L2mf','Ln_L2sf',
    'Ln_Wmf','Ln_Wsf','Ln_L1mm','Ln_L1sm','Ln_L2mm','Ln_L2sm',
    'Ln_Wmm','Ln_Wsm')
LG[1,]
for (i in 1:5){
  LG[i,1]=G[i,5];LG[i,2]=G[i,6]
  LG[i,3]=G[i,7];LG[i,4]=sqrt(G[i,6])*G[i,8]
  LG[i,5]=G[i,9];LG[i,6]=sqrt(G[i,6])*G[i,10]
  LG[i,7]=G[i,11];LG[i,8]=sqrt(G[i,6])*G[i,12]
  LG[i,9]=G[i,13];LG[i,10]=sqrt(G[i,6])*G[i,14]
  LG[i,11]=G[i,15];LG[i,12]=sqrt(G[i,6])*G[i,16]
  LG[i,13]=G[i,17];LG[i,14]=sqrt(G[i,6])*G[i,18]
}
for (i in 6:10){
  LG[i,1]=G[i,5];LG[i,2]=G[i,6]

  LG[i,3]=log(G[i,7])
  Stdev=sqrt(G[i,6])*G[i,8]
  LG[i,4]=Stdev/G[i,7]

  LG[i,5]=log(G[i,9])
  Stdev=sqrt(G[i,6])*G[i,10] 
  LG[i,6]=Stdev/G[i,9]

  LG[i,7]=log(G[i,11])
  Stdev=sqrt(G[i,6])*G[i,12]
  LG[i,8]=Stdev/G[i,11]

  LG[i,9]=log(G[i,13])
  Stdev=sqrt(G[i,6])*G[i,14]
  LG[i,10]=Stdev/G[i,13]

  LG[i,11]=log(G[i,15])
  Stdev=sqrt(G[i,6])*G[i,16] 
  LG[i,12]=Stdev/G[i,15]

  LG[i,13]=log(G[i,17])
  Stdev=sqrt(G[i,6])*G[i,18]
  LG[i,14]=Stdev/G[i,17]
}
LG
#========================================================Calculate rates
gentime=.2		#five generations per year: Gilchrist&al 2001, p. 276
idr=matrix(nrow=30,ncol=9)
  colnames(idr)=c('int','diff.sd','rate.sd.g','log.i','log|d|','log|r|',
    'sbn','wgt','sum')
nc=0
for (i in 1:6){	#no of traits
  for (j in 1:5){		#no of sites
    nc=nc+1
    idr[nc,1]=(LG[j,1]-LG[j+5,1])/gentime		#interval.g
    diff=LG[j,2*i+1]-LG[j+5,2*i+1]
    poolsd=PoolSD(LG[j,2],LG[j+5,2],LG[j,2*i+2],LG[j+5,2*i+2])#n1,n2,sd1,sd2
    idr[nc,2]=diff/poolsd				#diff.sd
    idr[nc,3]=idr[nc,2]/idr[nc,1]		#rate.sd.g   
    idr[nc,4]=log10(idr[nc,1])			#log.i
    idr[nc,5]=log10(abs(idr[nc,2]))			#log.d
    idr[nc,6]=log10(idr[nc,3])			#log.r
    idr[nc,7]=2						#sbn
    idr[nc,8]=1/idr[nc,1]				#wgt
    idr[nc,9]=idr[nc,4]+idr[nc,6]
  }
}
idr
#--------------------------------------------------------------------------
writefile<-"C://R_aaROEVchapt08//8.1.20_Gilchrist&al2004_Drosophila_SAm_out.csv"
write.csv(idr,file=writefile,na="NA",row.names=FALSE) 


