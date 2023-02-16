##Philip D. Gingerich RATES OF EVOLUTION (Cambridge University Press 2019) 
##8.1.06_Seeley1986_Littorina
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
#--------------------------------------------------Rate set-up
gentime=2				#see Hooks (2013, p. 2 for 2-year gen time)
idr=matrix(nrow=0,ncol=8)
  colnames(idr)=c('int','diff.sd','rate.sd.g','log.i','log|d|','log|r|',
    'sbn','wgt')	
#=============================================================== function
idrline<-function(S,gentime,yr1,yr2,yr3,sbn){
  #----------------------------------------------------------plot point
  plot(S$ind[S$year!=yr3],S$dep[S$year!=yr3],xlim=c(1.4,2.4),ylim=c(-2,2),
    xlab='ind',ylab='dep',pch=15,col='blue')
  points(S$ind[S$year==yr2],S$dep[S$year==yr2], pch=15, col='green')
  #------------------------------------------------analysis of covariance
  mdep=S$dep[S$year!=yr3];mind=S$ind[S$year!=yr3];myear=S$year[S$year!=yr3]
  model=lm(mdep~mind+myear);anova(model)
  stdev=sqrt(anova(model)[3,3]);stdev
  center=mean(S$ind[S$year!=yr3]);center
  regr1=lm(S$dep[S$year==yr1]~S$ind[S$year==yr1]);regr1
  mean1=regr1$coefficients[[2]]*center+regr1$coefficients[[1]];mean1
  regr2=lm(S$dep[S$year==yr2]~S$ind[S$year==yr2]);regr2
  mean2=regr2$coefficients[[2]]*center+regr2$coefficients[[1]];mean2
  id<-numeric(length=8)
  id[1]=(yr2-yr1)/gentime		#interval in generations
  id[2]=(mean2-mean1)/stdev		#diff.sd
  id[3]=abs(id[2])/id[1]		#rate.sd.g
  id[4]=log10(id[1])			#log.i
  id[5]=log10(abs(id[2]))		#log.d
  id[6]=log10(id[3])			#log.r
  id[7]=sbn					#sbn
  id[8]=1/id[1]				#wgt
  return(id)
}
#================================================== Nahant spire height
file1<-"C://R_aaROEVchapt08//8.1.06_Seeley1986_Littorina_1Nahant_spirehgt.csv"
S<-read.csv(file1,header=TRUE,sep=",",quote="\"",dec=".",fill=TRUE)
  colnames(S)=c('year','ind','dep');S[1:5,]
idr=rbind(idr,idrline(S,gentime,1898,1915,1983,2))#3rd year is one left out
idr=rbind(idr,idrline(S,gentime,1915,1983,1898,2))#3rd year is one left out
idr=rbind(idr,idrline(S,gentime,1898,1983,1915,3))#3rd year is one left out
#================================================== Nahant shell thickness
file2<-"C://R_aaROEVchapt08//8.1.06_Seeley1986_Littorina_2Nahant_shellthick.csv"
S<-read.csv(file2,header=TRUE,sep=",",quote="\"",dec=".",fill=TRUE)
  colnames(S)=c('year','ind','dep');S[1:5,]
idr=rbind(idr,idrline(S,gentime,1898,1915,1983,2))#3rd year is one left out
idr=rbind(idr,idrline(S,gentime,1915,1983,1898,2))#3rd year is one left out
idr=rbind(idr,idrline(S,gentime,1898,1983,1915,3))#3rd year is one left out
#================================================== Appledore spire height
file3<-"C://R_aaROEVchapt08//8.1.06_Seeley1986_Littorina_3Appledore_spirehgt.csv"
S<-read.csv(file3,header=TRUE,sep=",",quote="\"",dec=".",fill=TRUE)
  colnames(S)=c('year','ind','dep');S[1:5,]
idr=rbind(idr,idrline(S,gentime,1871,1983,1915,2))#3rd year is one left out
#================================================ Appledore shell thickness
file4<-"C://R_aaROEVchapt08//8.1.06_Seeley1986_Littorina_4Appledore_shellthick.csv"
S<-read.csv(file4,header=TRUE,sep=",",quote="\"",dec=".",fill=TRUE)
  colnames(S)=c('year','ind','dep');S[1:5,]
idr=rbind(idr,idrline(S,gentime,1871,1983,1915,2))#3rd year is one left out
#================================================= Isle au Haut spire height
file5<-"C://R_aaROEVchapt08//8.1.06_Seeley1986_Littorina_5IsleHaut_spirehgt.csv"
S<-read.csv(file5,header=TRUE,sep=",",quote="\"",dec=".",fill=TRUE)
  colnames(S)=c('year','ind','dep');S[1:5,]
idr=rbind(idr,idrline(S,gentime,1893,1983,1915,2))#3rd year is one left out
#============================================== Isle au Haut shell thickness
file6<-"C://R_aaROEVchapt08//8.1.06_Seeley1986_Littorina_6IsleHaut_shellthick.csv"
S<-read.csv(file6,header=TRUE,sep=",",quote="\"",dec=".",fill=TRUE)
  colnames(S)=c('year','ind','dep');S[1:5,]
idr=rbind(idr,idrline(S,gentime,1893,1983,1915,2))#3rd year is one left out
idr
#---------------------------------------------------------------------
writefile<-"C://R_aaROEVchapt08//8.1.06_Seeley1986_Littorina_out.csv"
write.csv(idr,file=writefile,na="NA",row.names=FALSE) 


