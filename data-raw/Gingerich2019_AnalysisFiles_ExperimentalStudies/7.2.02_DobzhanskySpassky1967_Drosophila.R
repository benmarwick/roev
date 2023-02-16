##Philip D. Gingerich RATES OF EVOLUTION (Cambridge University Press 2019) 
##7.2.02_DobzhanskySpassky1967_Drosophila
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
#======================================================= Plot panel a
xos=3.5;xoe=18.5;yos=9;yoe=18			#x,y original start and end
xfs=0;xfe=17;yfs=-1.5;yfe=1.5			#x,y foreground start and end
xfm=(xoe-xos)/(xfe-xfs)				#x foreground multiplier
yfm=(yoe-yos)/(yfe-yfs)				#y foreground multiplier
lines(c(xos,xoe),c(yos,yos),lwd=1,lty=1,col=1)
lines(c(xos,xos),c(yos,yoe),lwd=1,lty=1,col=1)

#---------------------------------------------------------axis labels
for (i in seq(yos,yoe-1,8/6)){					#y-axis
	lines(c(xos-.12,xos),c(i,i),lwd=1,lty=1,col=1)
}
for (i in seq(yos,yoe-1,8/6)){
	text(xos,i,format(round(.375*(i-13),digits=2),nsmall=2),pos=2,cex=1,col=1)
}
for (i in xos+.5+.95*xfm*c(xfs:xfe)){				#x-axis
	lines(c(i,i),c(yos,yos-.1),lwd=1,lty=1,col=1)
}
for (i in c(xfs:xfe)){
	text(xos+.5+.95*xfm*i,yos,i,pos=1,cex=1,col=1)
}
rect(xos-.5,yos+8.7,xos+15,yos+9,col='white',border=NA)
text(xos,yos+8.3,expression(paste(italic('Drosophila pseudoobscura')*': phototaxis')),
	pos=4,cex=1.3,col=1)
text(xos-1.9,yos+4,'Normalized phototactic score',
	srt=90,cex=1.3,col=1)
text(xos+7.8,yos-1.2,'Generation',
	cex=1.3,col=1)
#============================================================ Load file
file1<-"C://R_aaROEVchapt07//7.2.02_DobzhanskySpassky1967_Drosophila.csv"
ZS<-read.csv(file1,header=TRUE,sep=",",quote="\"",dec=".",fill=TRUE)
ZS
#---------------------------------------------- Rescale and plot Mhigh ZS[,7:8]
MH<-matrix(nrow=18,ncol=6)
  colnames(MH)=c('gen.','N','norm.mean','norm.sd','norm.m-','norm.m+')
MH[,1]=ZS[,1];MH[,2]=ZS[,2]
for (i in 1:18){
  MH[i,3]=NormScoreAsin(1,16,ZS[i,7],sqrt(ZS[i,8]))[1]
  MH[i,4]=NormScoreAsin(1,16,ZS[i,7],sqrt(ZS[i,8]))[2]
  MH[i,5]=MH[i,3]-MH[i,4];MH[i,6]=MH[i,3]+MH[i,4]
};MH
for(i in 1:18){
  lines(xos+.47+c(MH[i,1],MH[i,1])*.95*xfm,
    13+c(MH[i,5],MH[i,6])*(8/3),
    lty=1,lwd=1,col=gray(6/10))
}
for (i in 2:18){
  lines(xos+.5+c(MH[i-1,1],MH[i,1])*.95*xfm,
    13+c(MH[i-1,3],MH[i,3])*(8/3),
    lty=1,lwd=1,col=1)
}
#----------------------------------------------- Rescale and plot Mlow [,9:10]
ML<-matrix(nrow=18,ncol=6)
  colnames(ML)=c('gen.','N','norm.mean','norm.sd','norm.m-','norm.m+')
ML[,1]=ZS[,1];ML[,2]=ZS[,2]
for (i in 1:18){
  ML[i,3]=NormScoreAsin(1,16,ZS[i,9],sqrt(ZS[i,10]))[1]
  ML[i,4]=NormScoreAsin(1,16,ZS[i,9],sqrt(ZS[i,10]))[2]
  ML[i,5]=ML[i,3]-ML[i,4];ML[i,6]=ML[i,3]+ML[i,4]
};ML
for(i in 1:18){
  lines(xos+.47+c(ML[i,1],ML[i,1])*.95*xfm,
    13+c(ML[i,5],ML[i,6])*(8/3),
    lty=1,lwd=1,col=gray(6/10))
}
for (i in 2:18){
  lines(xos+.5+c(ML[i-1,1],ML[i,1])*.95*xfm,
    13+c(ML[i-1,3],ML[i,3])*(8/3),
    lty=1,lwd=1,col=1)
}
#---------------------------------------------- Rescale and plot Fhigh ZS[,3:4]
FH<-matrix(nrow=18,ncol=6)
  colnames(FH)=c('gen.','N','norm.mean','norm.sd','norm.m-','norm.m+')
FH[,1]=ZS[,1];FH[,2]=ZS[,2]
for (i in 1:18){
  FH[i,3]=NormScoreAsin(1,16,ZS[i,3],sqrt(ZS[i,4]))[1]
  FH[i,4]=NormScoreAsin(1,16,ZS[i,3],sqrt(ZS[i,4]))[2]
  FH[i,5]=FH[i,3]-FH[i,4];FH[i,6]=FH[i,3]+FH[i,4]
};FH
for(i in 1:18){
  lines(xos+.53+c(FH[i,1],FH[i,1])*.95*xfm,
    13+c(FH[i,5],FH[i,6])*(8/3),
    lty=1,lwd=1,col=gray(6/10))
}
for (i in 2:18){
  lines(xos+.5+c(FH[i-1,1],FH[i,1])*.95*xfm,
    13+c(FH[i-1,3],FH[i,3])*(8/3),
    lty=1,lwd=1,col=1)
}
#----------------------------------------------- Rescale and plot Flow [,5:6]
FL<-matrix(nrow=18,ncol=6)
  colnames(FL)=c('gen.','N','norm.mean','norm.sd','norm.m-','norm.m+')
FL[,1]=ZS[,1];FL[,2]=ZS[,2]
for (i in 1:18){
  FL[i,3]=NormScoreAsin(1,16,ZS[i,5],sqrt(ZS[i,6]))[1]
  FL[i,4]=NormScoreAsin(1,16,ZS[i,5],sqrt(ZS[i,6]))[2]
  FL[i,5]=FL[i,3]-FL[i,4];FL[i,6]=FL[i,3]+FL[i,4]
};FL
for(i in 1:18){
  lines(xos+.53+c(FL[i,1],FL[i,1])*.95*xfm,
    13+c(FL[i,5],FL[i,6])*(8/3),
    lty=1,lwd=1,col=gray(6/10))
}
for (i in 2:18){
  lines(xos+.5+c(FL[i-1,1],FL[i,1])*.95*xfm,
    13+c(FL[i-1,3],FL[i,3])*(8/3),
    lty=1,lwd=1,col=1)
}
points(xos+.5+MH[,1]*.95*xfm,13+MH[,3]*(8/3),pch=21,bg=1,cex=1,col=1)
points(xos+.5+ML[,1]*.95*xfm,13+ML[,3]*(8/3),pch=21,bg=1,cex=1,col=1)
points(xos+.5+FH[,1]*.95*xfm,13+FH[,3]*(8/3),pch=21,bg='white',cex=1,col=1)
points(xos+.5+FL[,1]*.95*xfm,13+FL[,3]*(8/3),pch=21,bg='white',cex=1,col=1)
#================================================ Rate calc Male high
MH[1:3,]
n=length(MH[,1]);n
nn=.5*(n-1)*n;nn
idr=matrix(nrow=nn,ncol=9)
  colnames(idr)=c('int','diff.sd','rate.sd','log.i','log.d','log.r',
    'sbn','wgt','sum')
nc=0
stdev<-numeric(length=n-1)
for (k in 1:(n-1)){    #run length
  for (i in 1:(n-k)){  #starting position
    nc=nc+1
    idr[nc,1]=k							#intercept
    meandiff=MH[(i+k),3]-MH[i,3]				#mean diff.
    poolsd=PoolSD(MH[i+k,2],MH[i,2],MH[i+k,4],MH[i,4])
    idr[nc,2]=meandiff/poolsd					#diff.sd
    idr[nc,3]=idr[nc,2]/idr[nc,1]				#rate.sd
    idr[nc,4]=log10(idr[nc,1])				#log.i
    idr[nc,5]=log10(abs(idr[nc,2]))				#log.d
    idr[nc,6]=log10(abs(idr[nc,3]))				#log.r
    if(idr[nc,1]==1){idr[nc,7]=1}else{idr[nc,7]=3}	#sbn
    idr[nc,8]=1/idr[nc,1]					#wgt
    idr[nc,9]=idr[nc,4]+idr[nc,6]				#sum
    stdev[i]=poolsd
  }
}
vartot=sum(stdev^2,na.rm=TRUE);varcnt=length(table(stdev,useNA='no'))
varmean=vartot/varcnt; avestdev=sqrt(varmean);avestdev
idr[1:3,]
idrx1=idr[!rowSums(!is.finite(idr)),]  #remove rows that have -Inf
idrx1[1:3,]
#================================================ Rate calc Female high
n=length(FH[,2]);n
nn=.5*(n-1)*n;nn
idr=matrix(nrow=nn,ncol=9)
  colnames(idr)=c('int','diff.sd','rate.sd','log.i','log.d','log.r',
    'sbn','wgt','sum')
nc=0
stdev<-numeric(length=n-1)
for (k in 1:(n-1)){    #run length
  for (i in 1:(n-k)){  #starting position
    nc=nc+1
    idr[nc,1]=k							#intercept
    meandiff=FH[(i+k),3]-FH[i,3]				#mean diff.
    poolsd=PoolSD(FH[i+k,2],FH[i,2],FH[i+k,4],FH[i,4])
    idr[nc,2]=meandiff/poolsd					#diff.sd
    idr[nc,3]=idr[nc,2]/idr[nc,1]				#rate.sd
    idr[nc,4]=log10(idr[nc,1])				#log.i
    idr[nc,5]=log10(abs(idr[nc,2]))				#log.d
    idr[nc,6]=log10(abs(idr[nc,3]))				#log.r
    if(idr[nc,1]==1){idr[nc,7]=1}else{idr[nc,7]=3}	#sbn
    idr[nc,8]=1/idr[nc,1]					#wgt
    idr[nc,9]=idr[nc,4]+idr[nc,6]				#sum
    stdev[i]=poolsd
  }
}
vartot=sum(stdev^2,na.rm=TRUE);varcnt=length(table(stdev,useNA='no'))
varmean=vartot/varcnt; avestdev=sqrt(varmean);avestdev
idr[1:3,]
idrx2=idr[!rowSums(!is.finite(idr)),]  #remove rows that have -Inf
idrxH=rbind(idrx1,idrx2)
proc.time()-ptm
#================================================ Rate calc Male low
ML[1:3,]
n=length(ML[,2]);n
nn=.5*(n-1)*n;nn
idr=matrix(nrow=nn,ncol=9)
  colnames(idr)=c('int','diff.sd','rate.sd','log.i','log.d','log.r',
    'sbn','wgt','sum')
nc=0
stdev<-numeric(length=n-1)
for (k in 1:(n-1)){    #run length
  for (i in 1:(n-k)){  #starting position
    nc=nc+1
    idr[nc,1]=k							#intercept
    meandiff=ML[(i+k),3]-ML[i,3]				#mean diff.
    poolsd=PoolSD(ML[i+k,2],ML[i,2],ML[i+k,4],ML[i,4])
    idr[nc,2]=meandiff/poolsd					#diff.sd
    idr[nc,3]=idr[nc,2]/idr[nc,1]				#rate.sd
    idr[nc,4]=log10(idr[nc,1])				#log.i
    idr[nc,5]=log10(abs(idr[nc,2]))				#log.d
    idr[nc,6]=log10(abs(idr[nc,3]))				#log.r
    if(idr[nc,1]==1){idr[nc,7]=1}else{idr[nc,7]=3}	#sbn
    idr[nc,8]=1/idr[nc,1]					#wgt
    idr[nc,9]=idr[nc,4]+idr[nc,6]				#sum
    stdev[i]=poolsd
  }
}
vartot=sum(stdev^2,na.rm=TRUE);varcnt=length(table(stdev,useNA='no'))
varmean=vartot/varcnt; avestdev=sqrt(varmean);avestdev
idr[1:3,]
idrx3=idr[!rowSums(!is.finite(idr)),]  #remove rows that have -Inf
idrx3[1:3,]
#================================================ Rate calc Female low
FL[1:3,]
n=length(FL[,2]);n
nn=.5*(n-1)*n;nn
idr=matrix(nrow=nn,ncol=9)
  colnames(idr)=c('int','diff.sd','rate.sd','log.i','log.d','log.r',
    'sbn','wgt','sum')
nc=0
stdev<-numeric(length=n-1)
for (k in 1:(n-1)){    #run length
  for (i in 1:(n-k)){  #starting position
    nc=nc+1
    idr[nc,1]=k							#intercept
    meandiff=FL[(i+k),3]-FL[i,3]				#mean diff.
    poolsd=PoolSD(FL[i+k,2],FL[i,2],FL[i+k,4],FL[i,4])
    idr[nc,2]=meandiff/poolsd					#diff.sd
    idr[nc,3]=idr[nc,2]/idr[nc,1]				#rate.sd
    idr[nc,4]=log10(idr[nc,1])				#log.i
    idr[nc,5]=log10(abs(idr[nc,2]))				#log.d
    idr[nc,6]=log10(abs(idr[nc,3]))				#log.r
    if(idr[nc,1]==1){idr[nc,7]=1}else{idr[nc,7]=3}	#sbn
    idr[nc,8]=1/idr[nc,1]					#wgt
    idr[nc,9]=idr[nc,4]+idr[nc,6]				#sum
    stdev[i]=poolsd
  }
}
vartot=sum(stdev^2,na.rm=TRUE);varcnt=length(table(stdev,useNA='no'))
varmean=vartot/varcnt; avestdev=sqrt(varmean);avestdev
idr[1:3,]
idrx4=idr[!rowSums(!is.finite(idr)),]  #remove rows that have -Inf
idrxL=rbind(idrx3,idrx4)
proc.time()-ptm
#========================================================== Write rate files
idrxA=rbind(idrxH,idrxL)
writefile<-"C://R_aaROEVchapt07//7.2.02_DobzhanskySpassky1967_Drosophila_out.csv"
write.csv(idrxA,file=writefile,na="NA",row.names=FALSE) 

#========================================================== Plot LRI panel b
#assign 'bootn' as boot number
bootn<-1000;text(-1.7,-2.6,paste("Boot n =",bootn),pos=4,cex=1,col=4)
#assign 'mode' as "medians","all","mixed"
mode<-"all";text(2,-2.6,paste("Mode: ",mode),pos=4,cex=1,col=4)
#assign circle size for points (1.5 or 2)
psize<-1.5	#2	
#assign 'equation' position as "normal","lower","none" at end of each call
#send (1)idrx matrix, (2)mode(diff/rate), (3)panel placement coordinate x, 
#  (4)panel placement coordinate y, (5)bootn, (6)mode, (7)psize, (8)equation
bootresultd=TriPanelBC(idrxH,"r",1,4.5,bootn,mode,psize,"lower")
wrlb.d=bootresultd[[1]];bootmat.d=bootresultd[[2]]
proc.time()-ptm
#========================================================== Plot LRI panel c
#send idrx matrix, n, mode(diff/rate), panel placement coordinates, mode
bootresultr=TriPanelBC(idrxL,"r",13,4.5,bootn,mode,psize,"lower")
wrlb.r=bootresultr[[1]];bootmat.r=bootresultr[[2]]
proc.time()-ptm
#========================================================== Label B and C
rect(-.5,6.7,21.5,7.5,border=NA,col='white')
rect(.1,6,9,6.7,border=NA,col='white')
rect(12.1,6,21,6.7,border=NA,col='white')
text(0,6.5,"Selection high lines",pos=4,cex=1.3,col=1)
text(12,6.5,"Selection low lines",pos=4,cex=1.3,col=1)

#----------------------------------------------------------
text(7,-2.6,paste("Processing time: ",round(proc.time()[3]-ptm[3],digits=2),
  "seconds"),pos=4,cex=1,col=4)


