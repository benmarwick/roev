##Philip D. Gingerich RATES OF EVOLUTION (Cambridge University Press 2019) 
##7.3.02_GoweFairfull1985_Gallus_eggs
#==========================================================================
cat(rep("\n",50))							#clear console
print (date()) 
rm(list=ls(all=TRUE))	#remove/clear all prev. variables
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
xos=2.5;xoe=17.5;yos=9;yoe=18	#x,y original start and end
xfs=0;xfe=30;yfs=3.6;yfe=5.0			#x,y foreground start and end
xfm=(xoe-xos)/(xfe-xfs)				#x foreground multiplier
yfm=(yoe-yos)/(yfe-yfs)				#y foreground multiplier
lines(c(xos,xoe),c(yos,yos),lwd=1,lty=1,col=1)
lines(c(xos,xos),c(yos,yoe),lwd=1,lty=1,col=1)

#---------------------------------------------------------axis labels
for (i in seq(yos,yoe-1,1)){
	lines(c(xos-.12,xos),c(i,i),lwd=1,lty=1,col=1)
}
for (i in seq(yos,yoe-1,1)){
	text(xos,i,format(round(.2*(i-8)+3.2,digits=2),nsmall=2),pos=2,cex=1,col=1)
}
for (i in xos+c(xfs:(xfe/1))*1*xfm){
	lines(c(i,i),c(yos,yos-.1),lwd=1,lty=1,col=1)
}
for (i in xos+c(xfs:(xfe/2))*2*xfm){
	text(i,yos,(i-2.5)/xfm,pos=1,cex=1,col=1)
}
text(xos,yos+8.8,expression(paste(italic('Gallus gallus')*': egg production')),
	pos=4,cex=1.3,col=1)
text(xos-1.9,yos+4,'Ln egg production',
	srt=90,cex=1.3,col=1)
text(xos+7.5,yos-1.2,'Generation',
	cex=1.3,col=1)
#============================================================ Load file
file1<-"C://R_aaROEVchapt07//7.3.02_GoweFairfull1985_Gallus_eggs.csv"
GF<-read.csv(file1,header=TRUE,sep=",",quote="\"",dec=".",fill=TRUE)
GF[1:3,]
#------------------------------------------------------------------------
GFL<-matrix(nrow=31,ncol=9)
  colnames(GFL)=c('Year','Gen.','N','Ln-3','SD-3','Ln-4','SD-4','Ln-C','SD-C')
for (i in 1:length(GF[,1])){
  GFL[i,1]=GF[i,1];GFL[i,2]=GF[i,2];GFL[i,3]=GF[i,3]
  GFL[i,4]=log(GF[i,4]);GFL[i,5]=sqrt(GF[i,5])/GF[i,4]
  GFL[i,6]=log(GF[i,6]);GFL[i,7]=sqrt(GF[i,7])/GF[i,6]
  GFL[i,8]=log(GF[i,8]);GFL[i,9]=sqrt(.5*(GFL[i,5]^2+GFL[i,7]^2))
}
##control variance averaged from lines 3 and 4
GFL[1:3,]
lines(c(2.5,17.5),5*c(GFL[1,8],GFL[1,8])-8,lty=2,lwd=1.4,col=gray(6/10))
#----------------------------------------------- Plot control line [,2&8-9]
for (i in 2:31){
  lines(2.5+c(GFL[i,2],GFL[i,2])*xfm,
    5*c(GFL[i,8]-GFL[i,9],GFL[i,8]+GFL[i,9])-8,
    lty=1,lwd=1,col=gray(6/10))
  lines(2.5+c(GFL[i-1,2],GFL[i,2])*xfm,
    5*c(GFL[i-1,8],GFL[i,8])-8,
    lty=1,lwd=1,col=1)
}
#------------------------------------------- Plot selection line 3 [,2&4-5]
for (i in 2:31){
  lines(2.45+c(GFL[i,2],GFL[i,2])*xfm,
    5*c(GFL[i,4]-GFL[i,5],GFL[i,4]+GFL[i,5])-8,
    lty=1,lwd=1,col=gray(6/10))
  lines(2.5+c(GFL[i-1,2],GFL[i,2])*xfm,
    5*c(GFL[i-1,4],GFL[i,4])-8,
    lty=1,lwd=1,col=1)
}
#------------------------------------------- Plot selection line 4 [,2&6-7]
for (i in 2:31){
  lines(2.55+c(GFL[i,2],GFL[i,2])*xfm,
    5*c(GFL[i,6]-GFL[i,7],GFL[i,6]+GFL[i,7])-8,
    lty=1,lwd=1,col=gray(6/10))
  lines(2.5+c(GFL[i-1,2],GFL[i,2])*xfm,
    5*c(GFL[i-1,6],GFL[i,6])-8,
    lty=1,lwd=1,col=1)
}
#------------------------------------------------Overlay points
points(2.5+GFL[,2]*xfm,
  5*GFL[,4]-8,
  pch=21,bg=1,cex=1,col=1)
points(2.5+GFL[,2]*xfm,
  5*GFL[,6]-8,
  pch=21,bg=1,cex=1,col=1)
points(2.5+GFL[,2]*xfm,
  5*GFL[,8]-8,
  pch=21,bg='white',cex=1.1,col=1)
GFL[20:21,3:9]<-NA;GFL	#remove years affected by Marek's disease outbreak
#============================================= Rate calc Strain 3 selection
GFL[1:3,]
n=length(GFL[,2]);n
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
    meandiff=GFL[(i+k),4]-GFL[i,4]				#mean diff.
    poolsd=PoolSD(GFL[i+k,3],GFL[i,3],GFL[i+k,5],GFL[i,5])
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
#============================================= Rate calc Strain 4 selection
GFL[1:3,]
n=length(GFL[,2]);n
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
    meandiff=GFL[(i+k),6]-GFL[i,6]				#mean diff.
    poolsd=PoolSD(GFL[i+k,3],GFL[i,3],GFL[i+k,7],GFL[i,7])
    idr[nc,2]=meandiff/poolsd					#diff.sd
    idr[nc,3]=idr[nc,2]/idr[nc,1]				#rate.sd
    idr[nc,4]=log10(idr[nc,1])				#log.i
    idr[nc,5]=log10(abs(idr[nc,2]))				#log.d
    idr[nc,6]=log10(abs(idr[nc,3]))				#log.r
    if(idr[nc,1]==1){idr[nc,7]=1}else{idr[nc,7]=3}	#sbn
    idr[nc,8]=1/idr[nc,1]					#wgt
    idr[nc,9]=idr[nc,4]+idr[nc,6]				#sum
    stdev[i]=poolsd					#wgt
  }
}
vartot=sum(stdev^2,na.rm=TRUE);varcnt=length(table(stdev,useNA='no'))
varmean=vartot/varcnt; avestdev=sqrt(varmean);avestdev
idr[1:3,]
idrx2=idr[!rowSums(!is.finite(idr)),]  #remove rows that have -Inf
idrx12=rbind(idrx1,idrx2)

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
bootresultd=TriPanelBC(idrx12,"r",1,4.5,bootn,mode,psize,"lower")
wrlb.d=bootresultd[[1]];bootmat.d=bootresultd[[2]]
proc.time()-ptm

#================================================ Rate calc GFLcontrol
GFL[1:3,]
n=length(GFL[,2]);n
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
    meandiff=GFL[(i+k),8]-GFL[i,8]				#mean diff.
    poolsd=PoolSD(GFL[i+k,3],GFL[i,3],GFL[i+k,9],GFL[i,9])
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
idrxC=idr[!rowSums(!is.finite(idr)),]  #remove rows that have -Inf
#========================================================== Write rate files
idrxA=rbind(idrx12,idrxC)
writefile<-"C://R_aaROEVchapt07//7.3.02_GoweFairfull1985_Gallus_eggs_out.csv"
write.csv(idrxA,file=writefile,na="NA",row.names=FALSE) 

#========================================================== Plot LRI panel c
#send idrx matrix, n, mode(diff/rate), panel placement coordinates, mode
bootresultr=TriPanelBC(idrxC,"r",13,4.5,bootn,mode,psize,"lower")
wrlb.r=bootresultr[[1]];bootmat.r=bootresultr[[2]]
proc.time()-ptm
#========================================================== Label B and C
rect(-.5,6.7,21.5,7.5,border=NA,col='white')
rect(.1,6,9,6.7,border=NA,col='white')
rect(12.1,6,21,6.7,border=NA,col='white')
text(0,6.5,"Selected lines",pos=4,cex=1.3,col=1)
text(12,6.5,"Control line",pos=4,cex=1.3,col=1)

#----------------------------------------------------------
text(7,-2.6,paste("Processing time: ",round(proc.time()[3]-ptm[3],digits=2),
  "seconds"),pos=4,cex=1,col=4)


