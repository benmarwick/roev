##Philip D. Gingerich RATES OF EVOLUTION (Cambridge University Press 2019) 
##8.1.32_Blanckenhorn2015_Scathophaga
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
xos=3.5;xoe=18;yos=9;yoe=17			#x,y original start and end
xfs=0;xfe=16;yfs=.6;yfe=1.4			#x,y foreground start and end
xfm=(xoe-xos)/(xfe-xfs)				#x foreground multiplier
yfm=(yoe-yos)/(yfe-yfs)				#y foreground multiplier
lines(c(xos,xoe),c(yos,yos),lwd=1,lty=1,col=1)
lines(c(xos,xos),c(yos,yoe),lwd=1,lty=1,col=1)

#---------------------------------------------------------axis labels
for (i in seq(yos,yoe-.5,1)){
  lines(c(xos-.12,xos),c(i,i),lwd=1,lty=1,col=1)
}
for (i in seq(yos,yoe-.5,1)){
  text(xos,i,format(round(.8/8*(i-9)+.6,digits=2),nsmall=2),pos=2,cex=1,col=1)
}
for (i in xos+c(xfs:xfe)*1*xfm){
  lines(c(i,i),c(yos,yos-.1),lwd=1,lty=1,col=1)
}
for (i in c(1994,1996,1998,2000,2002,2004,2006,2008)){
  text(xos+(i-1993)*xfm,yos,i,pos=1,cex=1,col=1)
}
text(xos,yos+8.0,expression(paste(italic('Scathophaga stercoraria')*
	': hind tibia length')),pos=4,cex=1.3,col=1)
text(xos-1.9,yos+4,'Ln length (mm)',
	srt=90,cex=1.3,col=1)
text(xos+7.5,yos-1.2,'Year',
	cex=1.3,col=1)

#============================================================ Load file
file1<-"C://R_aaROEVchapt08//8.1.32_Blanckenhorn2015_Scathophaga.csv"
LB<-read.csv(file1,header=TRUE,sep=",",quote="\"",dec=".",fill=TRUE);LB[1:3,]
LB[,1]=LB[,1]-1993
LB
sum(LB[1:10,3]);sum(LB[11:22,3])
#------------------------------------------------Plot gill raker time series
for (i in seq(yos,yoe-.5,1)){
  lines(xos+c(0,xfe)*xfm,c(i,i),lty=2,lwd=1,col=gray(6/10))
}
for (i in 1:10){
  lines(xos+c(LB[i,1],LB[i,1])*xfm,
    10*c(LB[i,4]-LB[i,5],LB[i,4]+LB[i,5])+3,
    lty=1,lwd=1,col=gray(6/10))
}
for (i in 11:22){
  lines(xos+c(LB[i,1],LB[i,1])*xfm,
    10*c(LB[i,4]-LB[i,5],LB[i,4]+LB[i,5])+3,
    lty=1,lwd=1,col=gray(6/10))
}
for (i in 2:10){
  lines(xos+c(LB[i-1,1],LB[i,1])*xfm,10*c(LB[i-1,4],LB[i,4])+3,
    lty=1,lwd=1,col=1)
}
for (i in 12:22){
  lines(xos+c(LB[i-1,1],LB[i,1])*xfm,10*c(LB[i-1,4],LB[i,4])+3,
    lty=1,lwd=1,col=1)
}
for (i in 1:nrow(LB)){
  if(LB[i,2]==1){
    points(xos+LB[i,1]*xfm,10*LB[i,4]+3,pch=21,cex=1.2,bg='white',col=1)
  }else{
    points(xos+LB[i,1]*xfm,10*LB[i,4]+3,pch=19,cex=1.1,col=1)
  }
}
#============================================= Rate calc: females
gentime=.25								#generation time
rect(3.6,9.2,9,9.6,col='white',border=NA)
text(3.5,9.4,"Generation time: 0.25 years",pos=4,cex=.9,col=1)
idrsum=0
LB[1:4,]
n=10;n
nn=.5*(n-1)*n;nn
idr=matrix(nrow=nn,ncol=9)
  colnames(idr)=c('int','diff.sd','rate.sd.g','log.i','log|d|','log|r|',
    'sbn','wgt','fgen')		#fgen is fraction of a generation
nc=0
for (k in 1:(n-1)){		#run length
  for (i in 1:(n-k)){  		#starting position
    nc=nc+1
    if(k<=gentime){idr[nc,1]=1}else{idr[nc,1]=(LB[(i+k),1]-LB[i,1])/gentime}
    meandiff=LB[(i+k),4]-LB[i,4]					#mean diff.
    poolsd=PoolSD(LB[i+k,3],LB[i,3],LB[i+k,5],LB[i,5])	#n1,n2,sd1,sd2
    idr[nc,2]=meandiff/poolsd						#diff.sd
    idr[nc,3]=abs(idr[nc,2])/idr[nc,1]				#rate.sd.gen
    idr[nc,4]=log10(idr[nc,1])					#log.i
    idr[nc,5]=log10(abs(idr[nc,2]))					#log.d
    idr[nc,6]=log10(idr[nc,3])					#log.r
    if(k==1){idr[nc,7]=2}else{idr[nc,7]=3}			#sbn
    if(k==1){idr[nc,7]=2}else{idr[nc,7]=3}			#sbn
    idr[nc,8]=1/idr[nc,1]						#wgt
    idr[nc,9]=k/gentime							#fract of gen.
  }
}
idr[1:3,]
idrsum=idrsum+nrow(idr)
idrx=idr[!rowSums(!is.finite(idr)),]  #remove rows that have -Inf
idrxF=idrx[,1:8]
idrc=idrx[idrx[,9]>=1,]  #keep only rows where interval >= gentime
idrcF=idrc[,1:8]
nrow(idr);nrow(idrx);nrow(idrxF);nrow(idrcF)
#============================================= Rate calc: males
idrsum=0
LB[1:4,]
n=12;n
nn=.5*(n-1)*n;nn
idr=matrix(nrow=nn,ncol=9)
  colnames(idr)=c('int','diff.sd','rate.sd.g','log.i','log|d|','log|r|',
    'sbn','wgt','fgen')		#fgen is fraction of a generation
nc=0
for (k in 1:(n-1)){		#run length
  for (i in 1:(n-k)){  		#starting position
    nc=nc+1
    if(k<=gentime){
      idr[nc,1]=1
    }else{
      idr[nc,1]=(LB[10+(i+k),1]-LB[10+i,1])/gentime
    }
    meandiff=LB[(10+i+k),4]-LB[10+i,4]				#mean diff.
    poolsd=PoolSD(LB[10+i+k,3],LB[10+i,3],LB[10+i+k,5],LB[10+i,5])#n1,n2,sd1,sd2
    idr[nc,2]=meandiff/poolsd						#diff.sd
    idr[nc,3]=abs(idr[nc,2])/idr[nc,1]				#rate.sd.gen
    idr[nc,4]=log10(idr[nc,1])					#log.i
    idr[nc,5]=log10(abs(idr[nc,2]))					#log.d
    idr[nc,6]=log10(idr[nc,3])					#log.r
    if(k==1){idr[nc,7]=2}else{idr[nc,7]=3}			#sbn
    if(k==1){idr[nc,7]=2}else{idr[nc,7]=3}			#sbn
    idr[nc,8]=1/idr[nc,1]						#wgt
    idr[nc,9]=k/gentime							#fract of gen.
  }
}
idr[1:3,]
idrsum=idrsum+nrow(idr)
idrx=idr[!rowSums(!is.finite(idr)),]  #remove rows that have -Inf
idrxM=idrx[,1:8]
idrc=idrx[idrx[,9]>=1,]  #keep only rows where interval >= gentime
idrcM=idrc[,1:8]
nrow(idr);nrow(idrx);nrow(idrxM);nrow(idrcM)
idrxA=rbind(idrxF,idrxM)
idrcA=rbind(idrcF,idrcM)
writefile<-"C://R_aaROEVchapt08//8.1.32_Blanckenhorn2015_Scathophaga_out.csv"
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
bootresultd=TriPanelBC(idrcA,"d",1,4.5,bootn,mode,psize,"lower")
wrlb.d=bootresultd[[1]];bootmat.d=bootresultd[[2]]
proc.time()-ptm
#========================================================== Plot LRI panel c
#send idrx matrix, n, mode(diff/rate), panel placement coordinates, mode
bootresultr=TriPanelBC(idrcA,"r",13,4.5,bootn,mode,psize,"lower")
wrlb.r=bootresultr[[1]];bootmat.r=bootresultr[[2]]
proc.time()-ptm
#========================================================== Label B and C
rect(-.5,6.7,21.5,7.5,border=NA,col='white')
rect(.1,6,9,6.7,border=NA,col='white')
rect(12.1,6,21,6.7,border=NA,col='white')
text(0,6.5,"Males and females",pos=4,cex=1.3,col=1)
text(12,6.5,"Males and females",pos=4,cex=1.3,col=1)

#----------------------------------------------------------
text(7,-2.6,paste("Processing time: ",round(proc.time()[3]-ptm[3],digits=2),
  "seconds"),pos=4,cex=1,col=4)


