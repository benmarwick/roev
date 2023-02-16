##Philip D. Gingerich RATES OF EVOLUTION (Cambridge University Press 2019) 
##9.2.33_Kimura&al2015_Murinae
#==========================================================================
cat(rep("\n",50))							#clear console
print (date()) #Ctr-s to save, Ctr-a to select all, Ctr-r to run
rm(list=ls(all=TRUE))#remove/clear all prev. variables
assign('last.warning',NULL,envir=baseenv())
ptm<-proc.time()
####################### Toggle traits to analyze ##########################
out='size';traits=c(3,5);traits				#analyze size traits
#out='shape';traits=c(9,11,13,15,17,19);traits		#analyze shape traits
##=======================================================================##
##library("devtools");library(roxygen2)
##setwd("c:/R_aaPackages/RATES");document();setwd("..");install("RATES") 
#=====================================================================Setup
library(RATES);library(MASS)
#----------------------------------------------------------------------Plot
win.graph(width=10, height=10, pointsize=12)
xr<-c(0,20);yr<-c(-2,18)		#xrange;yrange for plot axes
plot(xr,yr,					#set up plot
	type='n',				#type 'n' means no plotting
	pin=c(10,4),			#plot dimensions x,y in inches
	asp=1, 				#aspect ratio (y/x)
	col=1,				#color black
	las=1,				#axis labels always horizontal
	mgp=c(2,.3,0),			#margin for axis title/labels/tickline
	tck=-0.01,				#tick-mark length
	xaxp=c(xr[1],xr[2],xr[1]+xr[2]),		#x-axis ticks, no. ticks
	yaxp=c(yr[1],yr[2],abs(yr[1]-yr[2])),	#y-axis ticks, no. ticks
	cex.axis=.8,
	cex.lab=1.1)			#label size
dev.size('in')
#============================================================ Load file
file1<-"C://R_aaROEVchapt09//9.2.33_Kimura&al2015_Murinae.csv" 
K<-read.csv(file1,header=TRUE,sep=",",quote="\"",dec=".",fill=TRUE)
K[1:5,];nrow(K);ncol(K)
#------------------------------------------------------Logged matrices
KK=apply(as.matrix.noquote(K[1:nrow(K),]),2,as.numeric)
KK[1:3,]
#######################################################################
FillSK<-function(a,b,KK){
  #temp=KK[a:b,9:14]
  temp=KK[a:b,c(9:14,41:42)]
  x1=KK[a,5];x2=b-a+1
  x3=mean(log(temp[,1]));x4=sd(log(temp[,1]))
  x5=mean(log(temp[,2]));x6=sd(log(temp[,2]))
  x7=mean(log(temp[,1]*temp[,2]));x8=sd(log(temp[,1]*temp[,2]))
  x9=mean(temp[,3]);x10=sd(temp[,3])
  x11=mean(temp[,4]);x12=sd(temp[,4])
  x13=mean(temp[,5]);x14=sd(temp[,5])
  x15=mean(temp[,6]);x16=sd(temp[,6])
  x17=mean(temp[,7]);x18=sd(temp[,8])
  x19=mean(temp[,7]);x20=sd(temp[,8])
  return(c(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,
    x17,x18,x19,x20))
}
######################################################################
SK<-matrix(nrow=0,ncol=20)
  colnames(SK)=c('age','N','lnMl','lnSl','lnMw','lnSw','lnMa','lnSa',
    'Mvd','Svd','Mr','Sr','Ma1','Sa1','Ma2','Sa2',
    'Mp1','Sp1','Mp2','Sp2');SK
#---------------------------------------------------------------------
SK=rbind(SK,FillSK(301,304,KK))	#Antemus
SK=rbind(SK,FillSK(289,300,KK))
SK=rbind(SK,FillSK(284,288,KK))
SK=rbind(SK,FillSK(279,283,KK))
SK=rbind(SK,FillSK(275,278,KK))
SK=rbind(SK,FillSK(271,274,KK))
i=270;SK=rbind(SK,c(KK[i,5],1,log(KK[i,9]),NA,log(KK[i,10]),NA,log(KK[i,9]*KK[i,10]),NA,
  KK[i,11],NA,KK[i,12],NA,KK[i,13],NA,KK[i,14],NA,KK[i,41],NA,KK[i,42],NA) )
SK=rbind(SK,FillSK(268,269,KK))
SK=rbind(SK,FillSK(263,267,KK))	#Progonomys
SK=rbind(SK,FillSK(261,262,KK))
SK=rbind(SK,FillSK(259,260,KK))
i=258;SK=rbind(SK,c(KK[i,5],1,log(KK[i,9]),NA,log(KK[i,10]),NA,log(KK[i,9]*KK[i,10]),NA,
  KK[i,11],NA,KK[i,12],NA,KK[i,13],NA,KK[i,14],NA,KK[i,41],NA,KK[i,42],NA) )
SK=rbind(SK,FillSK(246,257,KK))
SK=rbind(SK,FillSK(237,239,KK))
SK=rbind(SK,FillSK(234,236,KK))
SK=rbind(SK,FillSK(226,231,KK))
SK=rbind(SK,FillSK(212,221,KK))
SK=rbind(SK,FillSK(190,192,KK))
SK=rbind(SK,FillSK(157,171,KK))
SK=rbind(SK,FillSK(120,122,KK))
SK=rbind(SK,FillSK(103,116,KK))
SK=rbind(SK,FillSK(84,94,KK))
SK=rbind(SK,FillSK(65,66,KK))
SK=rbind(SK,FillSK(58,63,KK))
SK=rbind(SK,FillSK(232,233,KK))
SK=rbind(SK,FillSK(222,225,KK))
SK=rbind(SK,FillSK(200,208,KK))
SK=rbind(SK,FillSK(196,199,KK))
SK=rbind(SK,FillSK(173,181,KK))
i=172;SK=rbind(SK,c(KK[i,5],1,log(KK[i,9]),NA,log(KK[i,10]),NA,log(KK[i,9]*KK[i,10]),NA,
  KK[i,11],NA,KK[i,12],NA,KK[i,13],NA,KK[i,14],NA,KK[i,41],NA,KK[i,42],NA) )
SK=rbind(SK,FillSK(123,156,KK))
SK=rbind(SK,FillSK(117,119,KK))
SK=rbind(SK,FillSK(97,101,KK))
SK=rbind(SK,FillSK(69,82,KK))
SK=rbind(SK,FillSK(31,45,KK))
SK=rbind(SK,FillSK(1,10,KK))
SK=rbind(SK,FillSK(49,57,KK))
SK=rbind(SK,FillSK(11,23,KK))
i=68;SK=rbind(SK,c(KK[i,5],1,log(KK[i,9]),NA,log(KK[i,10]),NA,log(KK[i,9]*KK[i,10]),NA,
  KK[i,11],NA,KK[i,12],NA,KK[i,13],NA,KK[i,14],NA,KK[i,41],NA,KK[i,42],NA) )
i=30;SK=rbind(SK,c(KK[i,5],1,log(KK[i,9]),NA,log(KK[i,10]),NA,log(KK[i,9]*KK[i,10]),NA,
  KK[i,11],NA,KK[i,12],NA,KK[i,13],NA,KK[i,14],NA,KK[i,41],NA,KK[i,42],NA) )
i=29;SK=rbind(SK,c(KK[i,5],1,log(KK[i,9]),NA,log(KK[i,10]),NA,log(KK[i,9]*KK[i,10]),NA,
  KK[i,11],NA,KK[i,12],NA,KK[i,13],NA,KK[i,14],NA,KK[i,41],NA,KK[i,42],NA) )
SK=rbind(SK,FillSK(24,28,KK))

rownames(SK)=c('An01','An02','An03','An04','An05','An06','An07','An08','Pr01','Pr02',
  'Pr03','Pr04','Pr05','Pr06','Pr07','Pr08','Pr09','Pr10','Pr11','Pr12',
  'Pr13','Pr14','Pr15','Pr16','Ka01','Ka02','Ka03','Ka04','Ka05','Ka06',
  'Ka07','Ka08','Ka09','Ka10','Ka11','Ka12','Mu01','Mu02','Pa01','Pa02',
  'Pa03','Pa04');SK;nrow(SK)

#======================================================= Plot panel a
xos=3;xoe=18;yos=9.5;yoe=17.5			#x,y original start and end
xfs=.1;xfe=1.6;yfs=-14.5;yfe=-6.5			#x,y foreground start and end
xfm=(xoe-xos)/(xfe-xfs)				#x foreground multiplier
yfm=(yoe-yos)/(yfe-yfs)				#y foreground multiplier
lines(c(xos,xoe),c(yos,yos),lwd=1,lty=1,col=1)
lines(c(xos,xos),c(yos,yoe+.4),lwd=1,lty=1,col=1)
#---------------------------------------------------------axis labels
for (i in seq(yos+.5,yoe,1)){				#y-axis
  lines(c(xos-.15,xos),c(i,i),lwd=1,lty=1,col=1)
}
for (i in seq(yos+.5,yoe,1)){
  text(xos+.05,i,format(round(abs((i-10)-14),digits=0),nsmall=0),pos=2,cex=1,col=1)
}
for (i in seq(xos,xoe,1)){		#x-axis
  lines(c(i,i),c(yos,yos-.15),lwd=1,lty=1,col=1)
}
for (i in seq(xos+1,xoe,2)){		#x-axis
  text(i,yos,format(round(-.2+i/10,digits=1),nsmall=1),
    pos=1,cex=1,col=1)
}
text(xos,yos+8.7,expression(paste('Murinae'~
  'spp.: molar'~M^1~'crown area')),pos=4,cex=1.3,col=1)
text(xos-1.4,yos+4,'Age (Ma)',
	srt=90,cex=1.3,col=1)
text(xos+7.5,yos-1.2,
  expression(paste('Ln (L x W) of'~M^1~'(mm)')),cex=1.3,col=1)
#---------------------------------------------------------------- Plot ranges
for (i in 4:18){
  lines(c(i,i),c(9.5,16.5),lty=2,lwd=1,col=rainbow(40)[9])
}
for (i in seq(yos+.5,yoe,.5)){
  lines(c(3,18),c(i,i),lty=2,lwd=1,col=rainbow(40)[9])
}
SK[1:5,]
xpos=-1;ypos=14.5
for (i in 1:nrow(SK)){
  lines(xos+10*c(SK[i,7]-SK[i,8],SK[i,7]+SK[i,8])+xpos,
    yos-1*c(SK[i,1],SK[i,1])+ypos,lty=1,lwd=1,col=gray(6/10))
}
for (i in 1:nrow(SK)){
  points(xos+10*SK[i,7]+xpos,yos-SK[i,1]+ypos,pch=20,cex=.8,col=1)
}
lines(xos+10*SK[1:8,7]+xpos,yos-SK[1:8,1]+ypos,lty=1,lwd=1.5,col=2)
lines(xos+10*SK[9:24,7]+xpos,yos-SK[9:24,1]+ypos,lty=1,lwd=1.5,col=3)
lines(xos+10*SK[25:36,7]+xpos,yos-SK[25:36,1]+ypos,lty=1,lwd=1.5,col=4)
lines(xos+10*SK[37:38,7]+xpos,yos-SK[37:38,1]+ypos,lty=1,lwd=1.5,col=6)
lines(xos+10*SK[39:42,7]+xpos,yos-SK[39:42,1]+ypos,lty=1,lwd=1.5,col=6)

#=============================================== Genus lineages
nrow(SK);rownames(SK)
GL1=SK[1:8,];nrow(GL1)		#An		
GL2=SK[9:24,];nrow(GL2)		#Pr
GL3=SK[25:36,];nrow(GL3)	#Ka
GL4=SK[37:38,];nrow(GL4)	#Mu
GL5=SK[39:42,];nrow(GL5)	#Pa


#=============================================== Rates: An
idr=matrix(nrow=0,ncol=9)
  colnames(idr)=c('int.g','diff.sd','rate.sd.g','log.i','log|d|','log|r|',
    'sbn','wgt','li+lr')
#------------------------------------------------ Rates
GL1[1:nrow(GL1),]
nr=nrow(GL1);nr						#number of rows
nn1=.5*(nr-1)*nr;nn1
zc1=0							#zero count
nc=0
for (cn in traits){
  for (i in 1:(nr-1)){				#increment
    for (sr in 1:(nr-i)){			#start row
      id<-numeric(length=9)
      nc=nc+1
      iyr=1000000*(GL1[sr,1]-GL1[sr+i,1])		#interval in years
      w1g=exp(2.723*GL1[sr,7]+1.6732)		#Rong.regr. from FreudSuar2013
      w2g=exp(2.723*GL1[sr+1,7]+1.6732)		#GSR body weight in grams
      wg=exp(.5*(log(w1g)+log(w2g)))		#pooled weight in grams
      gt=10^(.266*log10(wg)-.553)			#Eq 3.2 gen. time in years
      id[1]=iyr/gt					#interval in generations
      meandiff=GL1[sr+i,cn]-GL1[sr,cn]		#mean diff.
      if (!is.na(meandiff)){if (meandiff==0){zc1=zc1+1}}
      psd=PoolSD(GL1[sr+i,2],GL1[sr,2],GL1[sr+i,cn+1],GL1[sr,cn+1])#n1,n2,sd1,sd2
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
#--------------------------------------------------------------------------
nrow(idr)
#--------------------------------------------------------------------------
idrx1=idr[!rowSums(!is.finite(idr)),]  #remove rows that have -Inf
6*nn1;nrow(idr);nrow(idrx1);zc1
min(idrx1[,9]);max(idrx1[,9])

#=============================================== Rates: Pr
idr=matrix(nrow=0,ncol=9)
  colnames(idr)=c('int.g','diff.sd','rate.sd.g','log.i','log|d|','log|r|',
    'sbn','wgt','li+lr')
#------------------------------------------------ Rates
GL2[1:nrow(GL2),]
nr=nrow(GL2);nr						#number of rows
nn2=.5*(nr-1)*nr;nn2
zc2=0							#zero count
nc=0
for (cn in traits){
  for (i in 1:(nr-1)){				#increment
    for (sr in 1:(nr-i)){			#start row
      id<-numeric(length=9)
      nc=nc+1
      iyr=1000000*(GL2[sr,1]-GL2[sr+i,1])		#interval in years
      w1g=exp(2.723*GL2[sr,7]+1.6732)		#Rong.regr. from FreudSuar2013
      w2g=exp(2.723*GL2[sr+1,7]+1.6732)		#GSR body weight in grams
      wg=exp(.5*(log(w1g)+log(w2g)))		#pooled weight in grams
      gt=10^(.266*log10(wg)-.553)			#Eq 3.2 gen. time in years
      id[1]=iyr/gt					#interval in generations
      meandiff=GL2[sr+i,cn]-GL2[sr,cn]		#mean diff.
      if (!is.na(meandiff)){if (meandiff==0){zc2=zc2+1}}
      psd=PoolSD(GL2[sr+i,2],GL2[sr,2],GL2[sr+i,cn+1],GL2[sr,cn+1])#n1,n2,sd1,sd2
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
#--------------------------------------------------------------------------
nrow(idr)
#--------------------------------------------------------------------------
idrx2=idr[!rowSums(!is.finite(idr)),]  #remove rows that have -Inf
6*nn2;nrow(idr);nrow(idrx2);zc2
min(idrx2[,9]);max(idrx2[,9])

#=============================================== Rates: Ka
idr=matrix(nrow=0,ncol=9)
  colnames(idr)=c('int.g','diff.sd','rate.sd.g','log.i','log|d|','log|r|',
    'sbn','wgt','li+lr')
#------------------------------------------------ Rates
GL3[1:nrow(GL3),]
nr=nrow(GL3);nr						#number of rows
nn3=.5*(nr-1)*nr;nn3
zc3=0							#zero count
nc=0
for (cn in traits){
  for (i in 1:(nr-1)){				#increment
    for (sr in 1:(nr-i)){			#start row
      id<-numeric(length=9)
      nc=nc+1
      iyr=1000000*(GL3[sr,1]-GL3[sr+i,1])		#interval in years
      w1g=exp(2.723*GL3[sr,7]+1.6732)		#Rong.regr. from FreudSuar2013
      w2g=exp(2.723*GL3[sr+1,7]+1.6732)		#GSR body weight in grams
      wg=exp(.5*(log(w1g)+log(w2g)))		#pooled weight in grams
      gt=10^(.266*log10(wg)-.553)			#Eq 3.2 gen. time in years
      id[1]=iyr/gt					#interval in generations
      meandiff=GL3[sr+i,cn]-GL3[sr,cn]		#mean diff.
      if (!is.na(meandiff)){if (meandiff==0){zc3=zc3+1}}
      psd=PoolSD(GL3[sr+i,2],GL3[sr,2],GL3[sr+i,cn+1],GL3[sr,cn+1])#n1,n2,sd1,sd2
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
#--------------------------------------------------------------------------
nrow(idr)
#--------------------------------------------------------------------------
idrx3=idr[!rowSums(!is.finite(idr)),]  #remove rows that have -Inf
6*nn3;nrow(idr);nrow(idrx3);zc3
min(idrx3[,9]);max(idrx3[,9])

#=============================================== Rates: Mu
idr=matrix(nrow=0,ncol=9)
  colnames(idr)=c('int.g','diff.sd','rate.sd.g','log.i','log|d|','log|r|',
    'sbn','wgt','li+lr')
#------------------------------------------------ Rates
GL4[1:nrow(GL4),]
nr=nrow(GL4);nr						#number of rows
nn4=.5*(nr-1)*nr;nn4
zc4=0							#zero count
nc=0
for (cn in traits){
  for (i in 1:(nr-1)){				#increment
    for (sr in 1:(nr-i)){			#start row
      id<-numeric(length=9)
      nc=nc+1
      iyr=1000000*(GL4[sr,1]-GL4[sr+i,1])		#interval in years
      w1g=exp(2.723*GL4[sr,7]+1.6732)		#Rong.regr. from FreudSuar2013
      w2g=exp(2.723*GL4[sr+1,7]+1.6732)		#GSR body weight in grams
      wg=exp(.5*(log(w1g)+log(w2g)))		#pooled weight in grams
      gt=10^(.266*log10(wg)-.553)			#Eq 3.2 gen. time in years
      id[1]=iyr/gt					#interval in generations
      meandiff=GL4[sr+i,cn]-GL4[sr,cn]		#mean diff.
      if (!is.na(meandiff)){if (meandiff==0){zc4=zc4+1}}
      psd=PoolSD(GL4[sr+i,2],GL4[sr,2],GL4[sr+i,cn+1],GL4[sr,cn+1])#n1,n2,sd1,sd2
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
#--------------------------------------------------------------------------
nrow(idr)
#--------------------------------------------------------------------------
idrx4=idr[!rowSums(!is.finite(idr)),]  #remove rows that have -Inf
6*nn4;nrow(idr);nrow(idrx4);zc4
min(idrx4[,9]);max(idrx4[,9])

#=============================================== Rates: Pa
idr=matrix(nrow=0,ncol=9)
  colnames(idr)=c('int.g','diff.sd','rate.sd.g','log.i','log|d|','log|r|',
    'sbn','wgt','li+lr')
#------------------------------------------------ Rates
GL5[1:nrow(GL5),]
#gt=2	#generation time assumed
nr=nrow(GL5);nr						#number of rows
nn5=.5*(nr-1)*nr;nn5
zc5=0							#zero count
nc=0
for (cn in traits){
  for (i in 1:(nr-1)){				#increment
    for (sr in 1:(nr-i)){			#start row
      id<-numeric(length=9)
      nc=nc+1
      iyr=1000000*(GL5[sr,1]-GL5[sr+i,1])		#interval in years
      w1g=exp(2.723*GL5[sr,7]+1.6732)		#Rong.regr. from FreudSuar2013
      w2g=exp(2.723*GL5[sr+1,7]+1.6732)		#GSR body weight in grams
      wg=exp(.5*(log(w1g)+log(w2g)))		#pooled weight in grams
      gt=10^(.266*log10(wg)-.553)			#Eq 3.2 gen. time in years
      id[1]=iyr/gt					#interval in generations
      meandiff=GL5[sr+i,cn]-GL5[sr,cn]		#mean diff.
      if (!is.na(meandiff)){if (meandiff==0){zc5=zc5+1}}
      psd=PoolSD(GL5[sr+i,2],GL5[sr,2],GL5[sr+i,cn+1],GL5[sr,cn+1])#n1,n2,sd1,sd2
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
#--------------------------------------------------------------------------
nrow(idr)
#--------------------------------------------------------------------------
idrx5=idr[!rowSums(!is.finite(idr)),]  #remove rows that have -Inf
6*nn5;nrow(idr);nrow(idrx5);zc5
min(idrx5[,9]);max(idrx5[,9])



#--------------------------------------------------------------------Summary
idrxa=rbind(idrx1,idrx2,idrx3,idrx4,idrx5);idrx12=rbind(idrx1,idrx2)
6*nn1;6*nn2;6*nn3;6*nn4;6*nn5
nrow(idrx1);nrow(idrx2);nrow(idrx3);nrow(idrx4);nrow(idrx5);nrow(idrxa)
zc1;zc2;zc3;zc4;zc5;zca=zc1+zc2+zc3+zc4+zc5
zca;nrow(idrxa);zca+nrow(idrxa);nn1+nn2+nn3+nn4+nn5
if (out=='size'){
  writefile<-"C://R_aaROEVchapt09//9.2.33_Kimura&al2015_Murinae_size_out.csv"
}else{
  writefile<-"C://R_aaROEVchapt09//9.2.33_Kimura&al2015_Murinae_shape_out.csv"
}
write.csv(idrxa,file=writefile,na="NA",row.names=FALSE) 

#======================================================== Plot LRI panel (b)
#assign 'bootn' as boot number
bootn<-1000;text(-1.7,-2.6,paste("Boot n =",bootn),pos=4,cex=1,col=4)
#assign 'mode' as "medians","all","mixed"
mode<-"all";text(2,-2.6,paste("Mode: ",mode),pos=4,cex=1,col=4)
#assign circle size for points (1.5 or 2)
psize<-1.2	#1.5/2	
#assign 'equation' position as "normal","lower","none" at end of each call
#send (1)idrx matrix, (2)mode(diff/rate), (3)panel placement coordinate x, 
#  (4)panel placement coordinate y, (5)bootn, (6)mode, (7)psize, (8)equation
bootresultd=PalPanelBC(idrxa[,1:8],"r",1,4.5,bootn,mode,psize,"normal",0)
proc.time()-ptm

if(out=='size'){
  text(9.2,5,expression(paste("Size")),pos=2,cex=1,col=1)
}else{
  text(9.2,5,expression(paste("Shape")),pos=2,cex=1,col=1)
}

#----------------------------------------------------------
text(7,-2.6,paste("Processing time: ",round(proc.time()[3]-ptm[3],digits=2),
  "seconds"),pos=4,cex=1,col=4)
#print(bootcountd);print(bootcountr)



