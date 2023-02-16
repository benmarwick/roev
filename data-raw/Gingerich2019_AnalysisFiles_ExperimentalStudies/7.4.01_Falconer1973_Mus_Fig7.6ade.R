##Philip D. Gingerich RATES OF EVOLUTION (Cambridge University Press 2019) 
##7.4.01_Falconer1973_Mus_Fig7.6ade
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
xos=2.5;xoe=17.5;yos=9;yoe=16.5			#x,y original start and end
xfs=0;xfe=25;yfs=2.7;yfe=3.4			#x,y foreground start and end
xfm=(xoe-xos)/(xfe-xfs)				#x foreground multiplier
yfm=(yoe-yos)/(yfe-yfs)				#y foreground multiplier
lines(c(xos,xoe),c(yos,yos),lwd=1,lty=1,col=1)
lines(c(xos,xos),c(yos,yoe-.3),lwd=1,lty=1,col=1)

#---------------------------------------------------------axis labels
for (i in seq(yos,yoe-1,.5)){
	lines(c(xos-.12,xos),c(i,i),lwd=1,lty=1,col=1)
}
for (i in seq(yos,yoe-1,1)){
	text(xos,i,format(round(.2*(i-8)+2.0,digits=2),nsmall=2),pos=2,cex=1,col=1)
}
for (i in xos+c(xfs:(xfe-1/1))*1*xfm){
	lines(c(i,i),c(yos,yos-.1),lwd=1,lty=1,col=1)
}
for (i in xos+c(xfs:((xfe/2)-1))*2*xfm){
	text(i+xfm,yos,(i-2.5)/xfm,pos=1,cex=1,col=1)
}
text(xos,yos+7.0,expression(paste(italic('Mus musculus')*': 6-week body weight')),
	pos=4,cex=1.3,col=1)
text(xos-1.9,yos+3,'Ln weight (g)',
	srt=90,cex=1.3,col=1)
text(xos+7.5,yos-1.2,'Generation',
	cex=1.3,col=1)
#============================================================ Load file
file1<-"C://R_aaROEVchapt07//7.4.01_Falconer1973_Mus_wgt.csv"
F<-read.csv(file1,header=TRUE,sep=",",quote="\"",dec=".",fill=TRUE);GF[1:3,]
F[1:3,]
#-------------------------------------------------------------------------
FL<-matrix(nrow=24,ncol=22)
  colnames(FL)=c('Gen.','Csig','LnCA','LnCB','LnCC','LnCD','LnCE','LnCF',
    'Ssig','LnSA','LnSB','LnSC','LnSD','LnSE','LnSF',
    'Lsig','LnLA','LnLB','LnLC','LnLD','LnLE','LnLF')
Cstats<-matrix(nrow=24,ncol=5)
  colnames(Cstats)=c('gen','mean','stdev','min','max')
Sstats<-matrix(nrow=24,ncol=5)
  colnames(Sstats)=c('gen','mean','stdev','min','max')
Lstats<-matrix(nrow=24,ncol=5)
  colnames(Lstats)=c('gen','mean','stdev','min','max')
for (i in 1:length(FL[,1])){
  FL[i,1]=F[i,1];FL[i,2]=F[1,2]/mean(as.numeric(F[1,3:8]))
  Cstats[i,1]=FL[i,1]; Cstats[i,2]=log(mean(as.numeric(F[i,3:8])))
  Cstats[i,3]=FL[1,2]; Cstats[i,4]=log(min(as.numeric(F[i,3:8])))
  Cstats[i,5]=log(max(as.numeric(F[i,3:8])))
  FL[i,3]=log(F[i,3]);FL[i,4]=log(F[i,4]);FL[i,5]=log(F[i,5])
  FL[i,6]=log(F[i,6]);FL[i,7]=log(F[i,7]);FL[i,8]=log(F[i,8])
  FL[i,9]=F[1,9]/mean(as.numeric(F[1,10:15]))
  Sstats[i,1]=FL[i,1]; Sstats[i,2]=log(mean(as.numeric(F[i,10:15])))
  Sstats[i,3]=FL[1,2]; Sstats[i,4]=log(min(as.numeric(F[i,10:15])))
  Sstats[i,5]=log(max(as.numeric(F[i,10:15])))
  FL[i,10]=log(F[i,10]);FL[i,11]=log(F[i,11]);FL[i,12]=log(F[i,12])
  FL[i,13]=log(F[i,13]);FL[i,14]=log(F[i,14]);FL[i,15]=log(F[i,15])
  FL[i,16]=F[1,16]/mean(as.numeric(F[1,17:22]))
  Lstats[i,1]=FL[i,1]; Lstats[i,2]=log(mean(as.numeric(F[i,17:22])))
  Lstats[i,3]=FL[1,2]; Lstats[i,4]=log(min(as.numeric(F[i,17:22])))
  Lstats[i,5]=log(max(as.numeric(F[i,17:22])))
  FL[i,17]=log(F[i,17]);FL[i,18]=log(F[i,18]);FL[i,19]=log(F[i,19])
  FL[i,20]=log(F[i,20]);FL[i,21]=log(F[i,21]);FL[i,22]=log(F[i,22])
};FL[1:3,]
print(c(Cstats[1,2],Sstats[1,2],Lstats[1,2]))
print(c(Cstats[1,4],Sstats[1,4],Lstats[1,4]))
print(c(Cstats[1,5],Sstats[1,5],Lstats[1,5]))
#----------------------------------------------- Plot all 18 lines
Cxpoly=c(Cstats[1:24,1],Cstats[24:1,1]);Cypoly=c(Cstats[1:24,4],Cstats[24:1,5])
polygon(2.5+xfm+Cxpoly*xfm, 5*Cypoly-2.5, border=NA,col=gray(9/10))
Sxpoly=c(Sstats[1:24,1],Sstats[24:1,1]);Sypoly=c(Sstats[1:24,4],Sstats[24:1,5])
polygon(2.5+xfm+Sxpoly*xfm, 5*Sypoly-2.5, border=NA,col=gray(9/10))
Lxpoly=c(Cstats[1:24,1],Lstats[24:1,1]);Lypoly=c(Lstats[1:24,4],Lstats[24:1,5])
polygon(2.5+xfm+Lxpoly*xfm, 5*Lypoly-2.5, border=NA,col=gray(9/10))
for (i in c(3:8,10:15,17:22)){
  for (j in 2:24){
    lines(2.5+xfm+c(FL[j-1,1],FL[j,1])*xfm,
      5*c(FL[j-1,i],FL[j,i])-2.5,
      lty=1,lwd=1,col=gray(7/10))
  }
}
#----------------------------------------------- Plot mean control line 
for (i in 1:24){
  lines(2.5+xfm+c(Cstats[i,1],Cstats[i,1])*xfm,
    5*c(Cstats[i,2]-Cstats[i,3],Cstats[i,2]+Cstats[i,3])-2.5,
    lty=1,lwd=1,col=gray(6/10))
}
for (i in 2:24){
  lines(2.5+xfm+c(Cstats[i-1,1],Cstats[i,1])*xfm,
    5*c(Cstats[i-1,2],Cstats[i,2])-2.5,
    lty=1,lwd=1,col=1)
}
#------------------------------------------- Plot mean low selection 
for (i in 1:24){
  lines(2.46+xfm+c(Sstats[i,1],Sstats[i,1])*xfm,
    5*c(Sstats[i,2]-Sstats[i,3],Sstats[i,2]+Sstats[i,3])-2.5,
    lty=1,lwd=1,col=gray(6/10))
}
for (i in 2:24){
  lines(2.5+xfm+c(Sstats[i-1,1],Sstats[i,1])*xfm,
    5*c(Sstats[i-1,2],Sstats[i,2])-2.5,
    lty=1,lwd=1,col=1)
}
#------------------------------------------ Plot mean high selection line 
for (i in 1:24){
  lines(2.54+xfm+c(Lstats[i,1],Lstats[i,1])*xfm,
    5*c(Lstats[i,2]-Lstats[i,3],Lstats[i,2]+Lstats[i,3])-2.5,
    lty=1,lwd=1,col=gray(6/10))
}
for (i in 2:24){
  lines(2.5+xfm+c(Lstats[i-1,1],Lstats[i,1])*xfm,
    5*c(Lstats[i-1,2],Lstats[i,2])-2.5,
    lty=1,lwd=1,col=1)
}
#------------------------------------------------Overlay points
points(2.5+xfm+Sstats[,1]*xfm,
  5*Sstats[,2]-2.5,
  pch=21,bg=1,cex=1,col=1)
points(2.5+xfm+Lstats[,1]*xfm,
  5*Lstats[,2]-2.5,
  pch=21,bg=1,cex=1,col=1)
points(2.5+xfm+Cstats[,1]*xfm,
  5*Cstats[,2]-2.5,
  pch=21,bg='white',cex=1.1,col=1)
#============================================= Rate calc: low selection lines
FL[1:3,]
n=length(FL[,1]);n
nn=6*.5*(n-1)*n;nn
idr=matrix(nrow=nn,ncol=9)
  colnames(idr)=c('int','diff.sd','rate.sd','log.i','log.d','log.r',
    'sbn','wgt','sum')
nc=0
for (k in 1:(n-1)){    #run length
  for (i in 1:(n-k)){  #starting position
    for (rep in 10:15){
      nc=nc+1
      idr[nc,1]=k							#intercept
      meandiff=FL[(i+k),rep]-FL[i,rep]			#mean diff.
      poolsd=FL[1,9] #PoolSD(FL[i+k,3],FL[i,3],FL[i+k,5],FL[i,5])
      idr[nc,2]=meandiff/poolsd				#diff.sd
      idr[nc,3]=idr[nc,2]/idr[nc,1]				#rate.sd
      idr[nc,4]=log10(idr[nc,1])				#log.i
      idr[nc,5]=log10(abs(idr[nc,2]))			#log.d
      idr[nc,6]=log10(abs(idr[nc,3]))			#log.r
      if(idr[nc,1]==1){idr[nc,7]=1}else{idr[nc,7]=3}	#sbn
      idr[nc,8]=1/idr[nc,1]					#wgt
      idr[nc,9]=idr[nc,4]+idr[nc,6]				#sum
    }
  }
}
idr[1:3,]
idrxS=idr[!rowSums(!is.finite(idr)),]  #remove rows that have -Inf
#============================================= Rate calc: high selection lines
FL[1:3,]
n=length(FL[,1]);n
nn=6*.5*(n-1)*n;nn
idr=matrix(nrow=nn,ncol=9)
  colnames(idr)=c('int','diff.sd','rate.sd','log.i','log.d','log.r',
    'sbn','wgt','sum')
nc=0
for (k in 1:(n-1)){    #run length
  for (i in 1:(n-k)){  #starting position
    for (rep in 17:22){
      nc=nc+1
      idr[nc,1]=k							#intercept
      meandiff=FL[(i+k),rep]-FL[i,rep]			#mean diff.
      poolsd=FL[1,2] #PoolSD(FL[i+k,3],FL[i,3],FL[i+k,5],FL[i,5])
      idr[nc,2]=meandiff/poolsd				#diff.sd
      idr[nc,3]=idr[nc,2]/idr[nc,1]				#rate.sd
      idr[nc,4]=log10(idr[nc,1])				#log.i
      idr[nc,5]=log10(abs(idr[nc,2]))			#log.d
      idr[nc,6]=log10(abs(idr[nc,3]))			#log.r
      if(idr[nc,1]==1){idr[nc,7]=1}else{idr[nc,7]=3}	#sbn
      idr[nc,8]=1/idr[nc,1]					#wgt
      idr[nc,9]=idr[nc,4]+idr[nc,6]				#sum
    }
  }
}
idr[1:3,]
idrxL=idr[!rowSums(!is.finite(idr)),]  #remove rows that have -Inf
idrxSL=rbind(idrxS,idrxL)
#============================================= Rate calc: control lines
FL[1:3,]
n=length(FL[,1]);n
nn=6*.5*(n-1)*n;nn
idr=matrix(nrow=nn,ncol=9)
  colnames(idr)=c('int','diff.sd','rate.sd','log.i','log.d','log.r',
    'sbn','wgt','sum')
nc=0
for (k in 1:(n-1)){    #run length
  for (i in 1:(n-k)){  #starting position
    for (rep in 3:8){
      nc=nc+1
      idr[nc,1]=k							#intercept
      meandiff=FL[(i+k),rep]-FL[i,rep]			#mean diff.
      poolsd=FL[1,16] #PoolSD(FL[i+k,3],FL[i,3],FL[i+k,5],FL[i,5])
      idr[nc,2]=meandiff/poolsd				#diff.sd
      idr[nc,3]=idr[nc,2]/idr[nc,1]				#rate.sd
      idr[nc,4]=log10(idr[nc,1])				#log.i
      idr[nc,5]=log10(abs(idr[nc,2]))			#log.d
      idr[nc,6]=log10(abs(idr[nc,3]))			#log.r
      if(idr[nc,1]==1){idr[nc,7]=1}else{idr[nc,7]=3}	#sbn
      idr[nc,8]=1/idr[nc,1]					#wgt
      idr[nc,9]=idr[nc,4]+idr[nc,6]				#sum
    }
  }
}
idr[1:3,]
idrxC=idr[!rowSums(!is.finite(idr)),]  #remove rows that have -Inf
#=========================================--- Rate calc: averaged high lines
Lstats[1:3,]
n=length(Lstats[,1]);n
nn=.5*(n-1)*n;nn
idr=matrix(nrow=nn,ncol=9)
  colnames(idr)=c('int','diff.sd','rate.sd','log.i','log.d','log.r',
    'sbn','wgt','sum')
nc=0
for (k in 1:(n-1)){    #run length
  for (i in 1:(n-k)){  #starting position
    nc=nc+1
    idr[nc,1]=k							#intercept
    meandiff=Lstats[(i+k),2]-Lstats[i,2]			#mean diff.
    poolsd=FL[1,2] #PoolSD(FL[i+k,3],FL[i,3],FL[i+k,5],FL[i,5])
    idr[nc,2]=meandiff/poolsd					#diff.sd
    idr[nc,3]=idr[nc,2]/idr[nc,1]				#rate.sd
    idr[nc,4]=log10(idr[nc,1])				#log.i
    idr[nc,5]=log10(abs(idr[nc,2]))				#log.d
    idr[nc,6]=log10(abs(idr[nc,3]))				#log.r
      if(idr[nc,1]==1){idr[nc,7]=1}else{idr[nc,7]=3}	#sbn
      idr[nc,8]=1/idr[nc,1]					#wgt
      idr[nc,9]=idr[nc,4]+idr[nc,6]				#sum
  }
}
idr[1:3,]
idrxLA=idr[!rowSums(!is.finite(idr)),]  #remove rows that have -Inf
#============================================ Rate calc: averaged low lines
Sstats[1:3,]
n=length(Sstats[,1]);n
nn=.5*(n-1)*n;nn
idr=matrix(nrow=nn,ncol=9)
  colnames(idr)=c('int','diff.sd','rate.sd','log.i','log.d','log.r',
    'sbn','wgt','sum')
nc=0
for (k in 1:(n-1)){    #run length
  for (i in 1:(n-k)){  #starting position
    nc=nc+1
    idr[nc,1]=k							#intercept
    meandiff=Sstats[(i+k),2]-Sstats[i,2]			#mean diff.
    poolsd=FL[1,2] #PoolSD(FL[i+k,3],FL[i,3],FL[i+k,5],FL[i,5])
    idr[nc,2]=meandiff/poolsd					#diff.sd
    idr[nc,3]=idr[nc,2]/idr[nc,1]				#rate.sd
    idr[nc,4]=log10(idr[nc,1])				#log.i
    idr[nc,5]=log10(abs(idr[nc,2]))				#log.d
    idr[nc,6]=log10(abs(idr[nc,3]))				#log.r
      if(idr[nc,1]==1){idr[nc,7]=1}else{idr[nc,7]=3}	#sbn
      idr[nc,8]=1/idr[nc,1]					#wgt
      idr[nc,9]=idr[nc,4]+idr[nc,6]				#sum
  }
}
idr[1:3,]
idrxSA=idr[!rowSums(!is.finite(idr)),]  #remove rows that have -Inf
idrxLSA=rbind(idrxLA,idrxSA)
#========================================================== Write rate files
idrxA=rbind(idrxL,idrxS,idrxC)
#writefile<-"C://R_aaROEVchapt07//7.4.01_Falconer1973_Mus_wgt_out.csv"
#write.csv(idrxA,file=writefile,na="NA",row.names=FALSE) 

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
bootresultd=TriPanelBC(idrxLSA,"r",1,4.5,bootn,mode,psize,"lower")
wrlb.d=bootresultd[[1]];bootmat.d=bootresultd[[2]]
proc.time()-ptm
#========================================================== Plot LRI panel c
#send idrx matrix, n, mode(diff/rate), panel placement coordinates, mode
bootresultr=TriPanelBC(idrxC,"r",13,4.5,bootn,mode,psize,"lower")
wrlb.r=bootresultr[[1]];bootmat.r=bootresultr[[2]]
proc.time()-ptm
#========================================================== Label B and C
rect(-.5,6.7,21.5,7.5,border=NA,col='white')
rect(.1,6,9,6.7,border=NA,col='white')
rect(12.1,6,21,6.7,border=NA,col='white')
text(0,6.5,"Ave. high and low selection lines",pos=4,cex=1.3,col=1)
text(12,6.5,"Six control lines",pos=4,cex=1.3,col=1)

#----------------------------------------------------------
text(7,-2.6,paste("Processing time: ",round(proc.time()[3]-ptm[3],digits=2),
  "seconds"),pos=4,cex=1,col=4)
#print(bootcountd);print(bootcountr)
#======================================================== Likelihoods
#Last Draw.Likelihoods number is a vertical scale multiplier relative to 800
#Draw.Likelihoods(bootmat.d,bootmat.r,1)


