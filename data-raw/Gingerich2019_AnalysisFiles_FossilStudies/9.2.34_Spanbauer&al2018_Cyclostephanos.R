##Philip D. Gingerich RATES OF EVOLUTION (Cambridge University Press 2019)
##9.2.34_Spanbauer&al2018_Cyclostephanos
#==========================================================================
cat(rep("\n",50))		#clear console
print (date()) 		#Ctr-s to save, Ctr-a to select all, Ctr-r to run
rm(list=ls(all=TRUE))	#remove/clear all prev. variables
ptm<-proc.time()
##=================================================================== Tools
##library("devtools");library(roxygen2)
##setwd("c:/R_aaPackages/RATES");document();setwd("..");install("RATES")
#=====================================================================Setup
library(RATES);library(MASS)
#----------------------------------------------------------------------Plot

xr<-c(0,20);yr<-c(-2,18)		#xrange;yrange for plot axes
plot(xr,yr,					#set up plot
     ylab = "",
     xlab = "",
     xaxt='n',
     yaxt='n',
	type='n',				#type 'n' means no plotting
	pin=c(10,4),			#plot dimensions x,y in inches
	asp=1, 				#aspect ratio (y/x)
	col=1,				#color black
	las=1,				#axis labels always horizontal
	mgp=c(2,.3,0),			#margin for axis title/labels/tickline
	tck=-0.01,				#tick-mark length
	xaxp=c(xr[1],xr[2],xr[1]+xr[2]),	  	#x-axis ticks, no. ticks
	yaxp=c(yr[1],yr[2],abs(yr[1]-yr[2])),	#y-axis ticks, no. ticks
	cex.axis=.8,
	cex.lab=1.1)			#label size
dev.size('in')
#============================================================ Load file
file1<-"../Gingerich2019_AnalysisFiles_FossilStudies/9.2.34_Spanbauer&al2018_Cyclostephanos.csv"
S<-read.csv(file1,header=TRUE,sep=",",quote="\"",dec=".",fill=TRUE)
S[1:5,];nrow(S);ncol(S)
#------------------------------------------------------Logged matrix
SS=apply(as.matrix.noquote(S[1:nrow(S),]),2,as.numeric)
SS[1:3,]
LS<-matrix(nrow=nrow(SS),ncol=4)
  colnames(LS)=c('Time','N1','M1','S1')
LS[,1]=-1*SS[,5]                               # Time: time
LS[,2]=rowSums(!is.na(SS[,6:55]))              # N1: sum of measurements per row
LS[,3]=rowMeans(log(SS[,6:55]),na.rm=TRUE)     # M1: mean of logs of measurments per row
LS[,4]=apply(log(SS[,6:55]),1,sd,na.rm=TRUE)   # S1: sd of logs of measurement per row

min(LS[,1]);max(LS[,1])
min(na.omit(LS[,3]))-max(LS[,4])
max(na.omit(LS[,3]))+max(LS[,4])
(max(na.omit(LS[,3]))+max(LS[,4]))-(min(na.omit(LS[,3]))-max(LS[,4]))

#======================================================= Plot panel a
xos=2;xoe=17;yos=9.5;yoe=16.5			#x,y original start and end
xfs=1.8;xfe=4.8;yfs=-350;yfe=0			#x,y foreground start and end
xfm=(xoe-xos)/(xfe-xfs)				#x foreground multiplier
yfm=(yoe-yos)/(yfe-yfs)				#y foreground multiplier
lines(c(xos,xoe),c(yos,yos),lwd=1,lty=1,col=1)
lines(c(xos,xos),c(yos,yoe),lwd=1,lty=1,col=1)
#---------------------------------------------------------axis labels
for (i in seq(yos,yoe,1)){				#y-axis
  lines(c(xos-.15,xos),c(i,i),lwd=1,lty=1,col=1)
}
for (i in seq(yos,yoe,1)){
  text(xos+.1,i,format(round(abs(yfs+50*(i-9.5)),digits=0),nsmall=0),pos=2,cex=1,col=1)
}
for (i in xos-5*xfs+5*seq(xfs,xfe,.2)){		#x-axis
  lines(c(i,i),c(yos,yos-.15),lwd=1,lty=1,col=1)
}
for (i in xos-5*xfs+5*seq(xfs,xfe,.2)){		#x-axis
  #text(i,yos,3.1+.1*i,pos=1,cex=1,col=1)
  text(i,yos,format(round(1.4+.2*i,digits=1),nsmall=1),
    pos=1,cex=1,col=1)
}
text(xos,yos+7.6,expression(paste(italic('Cyclostephanos andinus')~
	':  valve diameter')),pos=4,cex=1.3,col=1)
text(xos-1.6,yos+3.5,'Age (Ka)',
	srt=90,cex=1.3,col=1)
text(xos+7.5,yos-1.2,
  expression(paste('Ln diameter'~'('~italic(mu)~')')),cex=1.3,col=1)
#---------------------------------------------------------------- Plot ranges
xpos=-9;ypos=7

py1=yos+.00002*(-220000)+ypos
py2=yos+.00002*(-70000)+ypos
lines(c(xos,xos+13.6),c(py1,py1),lty=2,lwd=1,col=gray(6/10))
lines(c(xos,xos+13.6),c(py2,py2),lty=2,lwd=1,col=gray(6/10))
text(xos+13.5,py1,'1st punctuation',pos=4,cex=1,col=1)
text(xos+13.5,py2,'2nd punctuation',pos=4,cex=1,col=1)

LS[1:5,]
min(na.omit(LS[,3]));max(na.omit(LS[,3]))
min(LS[,1]);max(LS[,1])
# plot mean log dia and sd log dia
for (i in 1:nrow(LS)){
  lines(xos+5*c(LS[i,3]-LS[i,4],LS[i,3]+LS[i,4])+xpos,
    yos+.00002*c(LS[i,1],LS[i,1])+ypos,lty=1,lwd=1,col=gray(6/10))
}
for (i in 1:nrow(LS)){
  points(xos+5*LS[i,3]+xpos,yos+.00002*LS[i,1]+ypos,pch=20,cex=.8,col=1)
}
lines(xos+5*LS[,3]+xpos,yos+.00002*LS[,1]+ypos,lty=1,lwd=1.5,col=1)

#==========================================================================
SP3=LS[1:130,1:4];SP3
SP2=LS[131:226,1:4];SP2
SP1=LS[227:266,1:4];SP1
TR2=LS[130:131,1:4];TR2
TR1=LS[226:227,1:4];TR1
#================================================ Rates SP1
idr=matrix(nrow=0,ncol=9)
  colnames(idr)=c('int.g','diff.sd','rate.sd.g','log.i','log|d|','log|r|',
    'sbn','wgt','li+lr')
#------------------------------------------------ Rates
gt=4	                      # generation time assumed
nr1=nrow(SP1);nr1						# number of rows
nn1=.5*(nr1-1)*nr1;nn1
zc1=0							#zero count
nc=0
SP1[1:4,]
cn=3
  for (i in 1:(nr1-1)){				#increment
    for (sr in 1:(nr1-i)){			#start row
      id<-numeric(length=9)
      nc=nc+1
      iyr=SP1[sr,1]-SP1[sr+i,1]		#interval in years
      id[1]=iyr/gt					#interval in generations
      meandiff=SP1[sr,cn]-SP1[sr+i,cn]		#mean diff.
      if (!is.na(meandiff)){if (meandiff==0){zc1=zc1+1}}
      psd=PoolSD(SP1[sr,cn-1],SP1[sr+i,cn-1],SP1[sr,cn+1],SP1[sr+i,cn+1])#n1,n2,sd1,sd2
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
#}
#--------------------------------------------------------------------------
nrow(idr)
#--------------------------------------------------------------------------
idrx1=idr[!rowSums(!is.finite(idr)),]  #remove rows that have -Inf
nn1;nrow(idr);nrow(idrx1);zc1
min(idrx1[,9]);max(idrx1[,9])

#================================================ Rates SP2
idr=matrix(nrow=0,ncol=9)
  colnames(idr)=c('int.g','diff.sd','rate.sd.g','log.i','log|d|','log|r|',
    'sbn','wgt','li+lr')
#------------------------------------------------ Rates
gt=4	#generation time assumed
nr2=nrow(SP2);nr2						#number of rows
nn2=.5*(nr2-1)*nr2;nn2
zc2=0							#zero count
nc=0
SP2[1:4,]
cn=3
  for (i in 1:(nr2-1)){				#increment
    for (sr in 1:(nr2-i)){			#start row
      id<-numeric(length=9)
      nc=nc+1
      iyr=SP2[sr,1]-SP2[sr+i,1]		#interval in years
      id[1]=iyr/gt					#interval in generations
      meandiff=SP2[sr,cn]-SP2[sr+i,cn]		#mean diff.
      if (!is.na(meandiff)){if (meandiff==0){zc2=zc2+1}}
      psd=PoolSD(SP2[sr,cn-1],SP2[sr+i,cn-1],SP2[sr,cn+1],SP2[sr+i,cn+1])#n1,n2,sd1,sd2
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
#--------------------------------------------------------------------------
nrow(idr)
#--------------------------------------------------------------------------
idrx2=idr[!rowSums(!is.finite(idr)),]  #remove rows that have -Inf
nn2;nrow(idr);nrow(idrx2);zc2
min(idrx2[,9]);max(idrx2[,9])

#================================================ Rates SP3
idr=matrix(nrow=0,ncol=9)
  colnames(idr)=c('int.g','diff.sd','rate.sd.g','log.i','log|d|','log|r|',
    'sbn','wgt','li+lr')
#------------------------------------------------ Rates
gt=4	#generation time assumed
nr3=nrow(SP3);nr3						#number of rows
nn3=.5*(nr3-1)*nr3;nn3
zc3=0							#zero count
nc=0
SP3[1:4,]
cn=3
  for (i in 1:(nr3-1)){				#increment
    for (sr in 1:(nr3-i)){			#start row
      id<-numeric(length=9)
      nc=nc+1
      iyr=SP3[sr,1]-SP3[sr+i,1]		#interval in years
      id[1]=iyr/gt					#interval in generations
      meandiff=SP3[sr,cn]-SP3[sr+i,cn]		#mean diff.
      if (!is.na(meandiff)){if (meandiff==0){zc3=zc3+1}}
      psd=PoolSD(SP3[sr,cn-1],SP3[sr+i,cn-1],SP3[sr,cn+1],SP3[sr+i,cn+1])#n1,n2,sd1,sd2
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
#--------------------------------------------------------------------------
nrow(idr)
#--------------------------------------------------------------------------
idrx3=idr[!rowSums(!is.finite(idr)),]  #remove rows that have -Inf
nn3;nrow(idr);nrow(idrx3);zc3
min(idrx3[,9]);max(idrx3[,9])

#================================================ Rates TR1
idr=matrix(nrow=0,ncol=9)
  colnames(idr)=c('int.g','diff.sd','rate.sd.g','log.i','log|d|','log|r|',
    'sbn','wgt','li+lr')
#------------------------------------------------ Rates
gt=4	#generation time assumed
nr4=nrow(TR1);nr4						#number of rows
nn4=.5*(nr4-1)*nr4;nn4
zc4=0							#zero count
nc=0
TR1[,1:4]
cn=3
  for (i in 1:(nr4-1)){				#increment
    for (sr in 1:(nr4-i)){			#start row
      id<-numeric(length=9)
      nc=nc+1
      iyr=TR1[sr,1]-TR1[sr+i,1]		#interval in years
      id[1]=iyr/gt					#interval in generations
      meandiff=TR1[sr,cn]-TR1[sr+i,cn]		#mean diff.
      if (!is.na(meandiff)){if (meandiff==0){zc4=zc4+1}}
      psd=PoolSD(TR1[sr,cn-1],TR1[sr+i,cn-1],TR1[sr,cn+1],TR1[sr+i,cn+1])#n1,n2,sd1,sd2
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
#--------------------------------------------------------------------------
nrow(idr)
#--------------------------------------------------------------------------
idrx4=idr[!rowSums(!is.finite(idr)),]  #remove rows that have -Inf
nn4;nrow(idr);nrow(idrx4);zc4
# min(idrx4[,9]);max(idrx4[,9])

#================================================ Rates TR2
idr=matrix(nrow=0,ncol=9)
  colnames(idr)=c('int.g','diff.sd','rate.sd.g','log.i','log|d|','log|r|',
    'sbn','wgt','li+lr')
#------------------------------------------------ Rates
gt=4	#generation time assumed
nr5=nrow(TR2);nr5						#number of rows
nn5=.5*(nr5-1)*nr5;nn5
zc5=0							#zero count
nc=0
TR2[,1:4]
cn=3
  for (i in 1:(nr5-1)){				#increment
    for (sr in 1:(nr5-i)){			#start row
      id<-numeric(length=9)
      nc=nc+1
      iyr=TR2[sr,1]-TR2[sr+i,1]		#interval in years
      id[1]=iyr/gt					#interval in generations
      meandiff=TR2[sr,cn]-TR2[sr+i,cn]		#mean diff.
      if (!is.na(meandiff)){if (meandiff==0){zc5=zc5+1}}
      psd=PoolSD(TR2[sr,cn-1],TR2[sr+i,cn-1],TR2[sr,cn+1],TR2[sr+i,cn+1])#n1,n2,sd1,sd2
      id[2]=meandiff/psd				#diff.sd
      id[3]=abs(id[2])/id[1]				#rate.sd.gen
      id[4]=log10(id[1])				#log.i
      id[5]=log10(abs(id[2]))				#log.d
      id[6]=log10(id[3])				#log.r
      if(i==1){id[7]=2}else{id[7]=3}		#sbn
      id[8]=1/id[1]					#wgt
      id[9]=id[4]+id[6]					#sum
    	#print(id)
    	idr=rbind(idr,id)
    }
  }
#--------------------------------------------------------------------------
nrow(idr)
#plot(idr[,4],idr[,6])
#--------------------------------------------------------------------------
idrx5=idr[!rowSums(!is.finite(idr)),]  #remove rows that have -Inf
nn5;nrow(idr);nrow(idrx5);zc5
min(idrx5[9]);max(idrx5[9])


#--------------------------------------------------------------Summary SP
idrxa=rbind(idrx1,idrx2,idrx3,idrx4,idrx5)
nrow(idrx1);nrow(idrx2);nrow(idrx3);nrow(idrxa)
zc1;zc2;zc3;zc1+zc2+zc3
zca=zc1+zc2+zc3
zca;nrow(idrxa);zca+nrow(idrxa);nn1+nn2+nn3
# writefile<-"C://R_aaROEVchapt09//9.2.34_Spanbauer&al2018_Cyclostephanos_out.csv"
# write.csv(idrxa,file=writefile,na="NA",row.names=FALSE)
#--------------------------------------------------------------Summary TR
idrxt=rbind(idrx4,idrx5)
#======================================================== Plot LRI panel (b)
#idrx12=rbind(idrx1,idrx2)
#assign 'bootn' as boot number
bootn<-1000;text(-1.7,-2.6,paste("Boot n =",bootn),pos=4,cex=1,col=4)
#assign 'mode' as "medians","all","mixed"
mode<-"all";text(2,-2.6,paste("Mode: ",mode),pos=4,cex=1,col=4)
#assign circle size for points (1.5 or 2)
psize<-1.2	#1.5/2
#assign 'equation' position as "normal","lower","none" at end of each call
#send (1)idrx matrix, (2)mode(diff/rate), (3)panel placement coordinate x,
#  (4)panel placement coordinate y, (5)bootn, (6)mode, (7)psize, (8)equation
#  (9) vadd or 'vertical addition' to position of points on plot
bootresultd=PalPanelBC(idrxa[,1:8],"d",1,4.5,bootn,mode,psize,"normal",0)
proc.time()-ptm
#======================================================== Plot LRI panel (c)
bootresultr=PalPanelBC(idrxa[,1:8],"r",13,4.5,bootn,mode,psize,"normal",1)
proc.time()-ptm

#======================================================== Add idrxt to plots
lines(c(idrxt[1,4]+1,idrxt[2,4]+3),c(idrxt[1,5],idrxt[2,5])+4.5,
  lty=1,lwd=2,col='white')
lines(c(idrxt[2,4]+1,idrxt[2,4]+3),c(idrxt[2,5],idrxt[2,5])+4.5,
  lty=1,lwd=2,col='white')
lines(c(idrxt[1,4]+1,idrxt[2,4]+3),c(idrxt[1,5],idrxt[2,5])+4.5,
  lty=1,lwd=1,col=1)
lines(c(idrxt[2,4]+1,idrxt[2,4]+3),c(idrxt[2,5],idrxt[2,5])+4.5,
  lty=1,lwd=1,col=1)
points(idrxt[,4]+1,idrxt[,5]+4.5,pch=19,cex=1.1,col=1)
text(idrxt[2,4]+2.9,idrxt[2,5]+4.55,'Punctuation',pos=4,cex=1,col=1)
text(idrxt[2,4]+3.0,idrxt[2,5]+4.05,'differences',pos=4,cex=1,col=1)

lines(c(idrxt[1,4]+13,idrxt[2,4]+15),c(idrxt[1,6],idrxt[2,6])+7.5,
  lty=1,lwd=2,col='white')
lines(c(idrxt[2,4]+13,idrxt[2,4]+15),c(idrxt[2,6],idrxt[2,6])+7.5,
  lty=1,lwd=2,col='white')
lines(c(idrxt[1,4]+13,idrxt[2,4]+15),c(idrxt[1,6],idrxt[2,6])+7.5,
  lty=1,lwd=1,col=1)
lines(c(idrxt[2,4]+13,idrxt[2,4]+15),c(idrxt[2,6],idrxt[2,6])+7.5,
  lty=1,lwd=1,col=1)
points(idrxt[,4]+13,idrxt[,6]+7.5,pch=19,cex=1.1,col=1)
text(idrxt[2,4]+14.9,idrxt[2,6]+7.55,'Punctuation',pos=4,cex=1,col=1)
text(idrxt[2,4]+15.6,idrxt[2,6]+7.05,'rates',pos=4,cex=1,col=1)

#----------------------------------------------------------
text(7,-2.6,paste("Processing time: ",round(proc.time()[3]-ptm[3],digits=2),
  "seconds"),pos=4,cex=1,col=4)


