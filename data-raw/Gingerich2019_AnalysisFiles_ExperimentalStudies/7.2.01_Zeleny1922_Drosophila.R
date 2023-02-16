##Philip D. Gingerich RATES OF EVOLUTION (Cambridge University Press 2019)
##7.2.01_Zeleny1922_Drosophila
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
##setwd("c:/R_aaPackages/ROEV");document();setwd("..");install("ROEV")
#=====================================================================Setup
library(ROEV);library(MASS)
#----------------------------------------------------------------------Plot
win.graph(width=10, height=10, pointsize=12)
#win.graph(width=10, height=6.5, pointsize=12)
#Device size is size of window on the screen.
#dev.new(width=10,height=5.5)		#set device size in inches
xr<-c(0,20);yr<-c(-2,18)		#xrange;yrange for plot axes
plot(xr,yr,					#set up plot
	type='n',				#type 'n' means no plotting
	pin=c(10,4),			#plot dimensions x,y in inches
	asp=1, 				#aspect ratio (y/x)
	#asp=16*.5/1400,			#aspect ratio asp=x/y)
	col=1,				#color black
	las=1,				#axis labels always horizontal
	#mgp=c(2.2,.3,0),			#margin for axis title/labels/tickline
	mgp=c(2,.3,0),			#margin for axis title/labels/tickline
	tck=-0.01,				#tick-mark length
#	xlab=expression('Log'[10]*' body weight (g)'),
#	xaxp=c(27,53,26),
	xaxp=c(xr[1],xr[2],xr[1]+xr[2]),
#	ylab=expression('Log'[10]*' generation time (yr)'),
	yaxp=c(yr[1],yr[2],abs(yr[1]-yr[2])),		#y-axis extreme ticks and number of ticks
	cex.axis=.8,
	cex.lab=1.1)			#label size
dev.size('in')
#======================================================= Plot panel a
xos=2.5;xoe=17.5;yos=9;yoe=18	#x,y original start and end
xfs=0;xfe=42;yfs=3.5;yfe=6		#x,y foreground start and end
xfm=(xoe-xos)/(xfe-xfs)				#x foreground multiplier
yfm=(yoe-yos)/(yfe-yfs)				#y foreground multiplier
lines(c(xos,xoe),c(yos,yos),lwd=1,lty=1,col=1)
lines(c(xos,xos),c(yos,yoe),lwd=1,lty=1,col=1)

#---------------------------------------------------------axis labels
for (i in seq(yos,yoe-1,1)){
	lines(c(xos-.12,xos),c(i,i),lwd=1,lty=1,col=1)
}
for (i in seq(yos,yoe-1,1)){
	text(xos,i,format(round(.25*(i-8)+3.25,digits=2),nsmall=2),pos=2,cex=1,col=1)
}
for (i in xos+c(xfs:(xfe/5))*5*xfm){
	lines(c(i,i),c(yos,yos-.1),lwd=1,lty=1,col=1)
}
for (i in xos+c(xfs:(xfe/5))*5*xfm){
	text(i,yos,(i-2.5)/xfm,pos=1,cex=1,col=1)
}
text(xos,yos+8.8,expression(paste(italic('Drosophila melanogaster')*': eye facet number')),
	pos=4,cex=1.3,col=1)
text(xos-1.8,yos+4,'Eye facet number (ln count)',
	srt=90,cex=1.3,col=1)
text(xos+7.5,yos-1.2,'Generation',
	cex=1.3,col=1)
#============================================================ Load file
file1<-"C://R_aaROEVchapt07//7.2.01_Zeleny1922_Drosophila_male_high.csv"
file2<-"C://R_aaROEVchapt07//7.2.01_Zeleny1922_Drosophila_male_low.csv"
file3<-"C://R_aaROEVchapt07//7.2.01_Zeleny1922_Drosophila_female_high.csv"
file4<-"C://R_aaROEVchapt07//7.2.01_Zeleny1922_Drosophila_female_low.csv"
ZMH<-read.csv(file1,header=TRUE,sep=",",quote="\"",dec=".",fill=TRUE);ZMH
ZML<-read.csv(file2,header=TRUE,sep=",",quote="\"",dec=".",fill=TRUE);ZML
ZFH<-read.csv(file3,header=TRUE,sep=",",quote="\"",dec=".",fill=TRUE);ZFH
ZFL<-read.csv(file4,header=TRUE,sep=",",quote="\"",dec=".",fill=TRUE);ZFL
#----------------------------------------------- Replace and plot ZMH[,6:7]
ZMH[,7]=.0953*ZMH[,5]+4.7131
ZMH[,8]=.5*((.0953*(ZMH[,5]+ZMH[,6])+4.7131)-
  (.0953*(ZMH[,5]-ZMH[,6])+4.7131));ZMH
for (i in 2:43){
  lines(2.475+c(ZMH[i,1],ZMH[i,1])*xfm,
    4*c(ZMH[i,7]-ZMH[i,8],ZMH[i,7]+ZMH[i,8])-5,
    lty=1,lwd=1,col=gray(6/10))
  lines(2.5+c(ZMH[i-1,1],ZMH[i,1])*xfm,4*c(ZMH[i-1,7],ZMH[i,7])-5,
    lty=1,lwd=1,col=1)
}
#----------------------------------------------- Replace and plot ZML[,6:7]
ZML[,7]=.0953*ZML[,5]+4.7131
ZML[,8]=.5*((.0953*(ZML[,5]+ZML[,6])+4.7131)-
  (.0953*(ZML[,5]-ZML[,6])+4.7131));ZML
for (i in 2:43){
  lines(2.475-.01+c(ZML[i,1],ZML[i,1])*xfm,
    4*c(ZML[i,7]-ZML[i,8],ZML[i,7]+ZML[i,8])-5,
    lty=1,lwd=1,col=gray(6/10))
  lines(2.5+c(ZML[i-1,1],ZML[i,1])*xfm,4*c(ZML[i-1,7],ZML[i,7])-5,
    lty=1,lwd=1,col=1)
}
#----------------------------------------------- Replace and plot ZFH[,6:7]
ZFH[,7]=.0953*ZFH[,5]+4.0741
ZFH[,8]=.5*((.0953*(ZFH[,5]+ZFH[,6])+4.7131)-
  (.0953*(ZFH[,5]-ZFH[,6])+4.7131));ZFH
for (i in 2:43){
  lines(2.525+c(ZFH[i,1],ZFH[i,1])*xfm,
    4*c(ZFH[i,7]-ZFH[i,8],ZFH[i,7]+ZFH[i,8])-5,
    lty=1,lwd=1,col=gray(6/10))
  lines(2.5+c(ZFH[i-1,1],ZFH[i,1])*xfm,4*c(ZFH[i-1,7],ZFH[i,7])-5,
    lty=1,lwd=1,col=1)
}
#----------------------------------------------- Replace and plot ZFL[,6:7]
ZFL[,7]=.0953*ZFL[,5]+4.0741
ZFL[,8]=.5*((.0953*(ZFL[,5]+ZFL[,6])+4.7131)-
  (.0953*(ZFL[,5]-ZFL[,6])+4.7131));ZFL
for (i in 2:43){
  lines(2.525+c(ZFL[i,1],ZFL[i,1])*xfm,
    4*c(ZFL[i,7]-ZFL[i,8],ZFL[i,7]+ZFL[i,8])-5,
    lty=1,lwd=1,col=gray(6/10))
  lines(2.5+c(ZFL[i-1,1],ZFL[i,1])*xfm,4*c(ZFL[i-1,7],ZFL[i,7])-5,
    lty=1,lwd=1,col=1)
}
points(2.5+ZMH[,1]*xfm,4*ZMH[,7]-5,pch=21,bg=1,cex=1,col=1)
points(2.5+ZML[,1]*xfm,4*ZML[,7]-5,pch=21,bg=1,cex=1,col=1)
points(2.5+ZFH[,1]*xfm,4*ZFH[,7]-5,pch=21,bg='white',cex=1,col=1)
points(2.5+ZFL[,1]*xfm,4*ZFL[,7]-5,pch=21,bg='white',cex=1,col=1)
#================================================ Rate calc ZMH
ZMH[1:3,]
n=length(ZMH[,1]);n
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
    meandiff=ZMH[(i+k),7]-ZMH[i,7]				#mean diff.
    poolsd=PoolSD(ZMH[i+k,2],ZMH[i,2],ZMH[i+k,8],ZMH[i,8])
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
#================================================ Rate calc ZFH
n=length(ZFH[,2]);n
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
    meandiff=ZFH[(i+k),7]-ZFH[i,7]				#mean diff.
    poolsd=PoolSD(ZFH[i+k,2],ZFH[i,2],ZFH[i+k,8],ZFH[i,8])
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
#================================================ Rate calc ZML
ZML[1:3,]
n=length(ZML[,2]);n
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
    meandiff=ZML[(i+k),7]-ZML[i,7]				#mean diff.
    poolsd=PoolSD(ZML[i+k,2],ZML[i,2],ZML[i+k,8],ZML[i,8])
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
#================================================ Rate calc ZFL
n=length(ZFL[,2]);n
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
    meandiff=ZFL[(i+k),7]-ZFL[i,7]				#mean diff.
    poolsd=PoolSD(ZFL[i+k,2],ZFL[i,2],ZFL[i+k,8],ZFL[i,8])
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
writefile<-"C://R_aaROEVchapt07//7.2.01_Zeleny1922_Drosophila_out.csv"
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


