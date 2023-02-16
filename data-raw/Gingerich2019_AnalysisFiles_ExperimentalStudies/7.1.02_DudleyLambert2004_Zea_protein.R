##Philip D. Gingerich RATES OF EVOLUTION (Cambridge University Press 2019)
##7.1.02_DudleyLambert2004_Zea_protein
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
xos=2.5;xoe=17.5;yos=10;yoe=18		#x,y original start and end
xfs=0;xfe=100;yfs=.2;yfe=.4			#x,y foreground start and end
xfm=(xoe-xos)/(xfe-xfs)				#x foreground multiplier
yfm=(yoe-yos)/(yfe-yfs)				#y foreground multiplier
lines(c(xos,xoe),c(yos,yos),lwd=1,lty=1,col=1)
lines(c(xos,xos),c(yos,yoe),lwd=1,lty=1,col=1)

#---------------------------------------------------------axis labels
for (i in seq(yos,yoe-1,1)){
	lines(c(xos-.12,xos),c(i,i),lwd=1,lty=1,col=1)
}
for (i in seq(yos,yoe-1,2)){
	text(xos,i,format(round(2.5*i-5,digits=0),nsmall=0),pos=2,cex=1,col=1)
}
for (i in xos+c(xfs:(xfe/5))*5*xfm){
	lines(c(i,i),c(yos,yos-.1),lwd=1,lty=1,col=1)
}
for (i in xos+c(xfs:(xfe/10))*10*xfm){
	text(i,yos,(i-2.5)/xfm,pos=1,cex=1,col=1)
}
text(xos,yos+7.8,expression(paste(italic('Zea mays')*': Illinois High Protein')),
	pos=4,cex=1.3,col=1)
text(xos-1.3,yos+4,'Protein (normalized %)',
	srt=90,cex=1.3,col=1)
text(xos+7.5,yos-1.2,'Generation',
	cex=1.3,col=1)
#---------------------------------------------------------Load file
file<-"../Gingerich2019_AnalysisFiles_ExperimentalStudies/7.1.02_DudleyLambert2004_Zea_protein.csv"
TT<-read.csv(file,header=TRUE,sep=",",quote="\"",dec=".",fill=TRUE)
TT[1:3,]
#------------------------------------------------ Selection differentials
TD<-matrix(nrow=101,ncol=10)
  colnames(TD)=c('Year','Gen','N','pcMean','pcSeed','pcSD',
    'Mean','Seed','Diff','SD')
TD[,1]=TT[,1]; TD[,2]=TT[,2]; TD[,3]=TT[,3]; TD[,4]=TT[,5];
TD[,5]=TT[,5]+TT[,6] ; TD[,6]=TT[,7]; TD[,7]=Pc2Ar(TD[,4]);
TD[,8]=Pc2Ar(TD[,5]) ; TD[,9]=TD[,8]-TD[,7]; TD[,10]=Pcsd2Arsd(TT[,5],TT[,7])
TD[1:3,]
mean(TD[,9]);mean(TD[,10]);mean(TD[,9])/mean(TD[,10])
median(TD[,9]);median(TD[,10]);median(TD[,9])/median(TD[,10])

#======================================================== Read file
TS<-matrix(nrow=101,ncol=9)
  colnames(TS)=c('Year','Gen','N','pcMean','pcSD','Mean','SD','M-SD','M+SD')
TS[,1]=TT[,1]; TS[,2]=TT[,2]; TS[,3]=TT[,3]; TS[,4]=TT[,5]; TS[,5]=TT[,7]
TS[,6]=Pc2Ar(TT[,5]); TS[,7]=Pcsd2Arsd(TT[,5],TT[,7])
TS[,8]=TS[,6]-TS[,7]; TS[,9]=TS[,6]+TS[,7]
TS[1:3,];TS[99:101,]
TS[101,6]; TS[1,6]; TS[101,6]-TS[1,6]
mean(TS[,7])
#===================================================== Plot values
for (i in 2:length(TS[,2])){
  lines(xos+c(TS[(i-1),2],TS[i-1,2])*xfm,
  2+c(TS[(i-1),8],TS[i-1,9])*yfm,
  lty=1,lwd=1,col=gray(6/10))
}
for (i in 2:length(TS[,2])){
  lines(xos+c(TS[(i-1),2],TS[i,2])*xfm,2+c(TS[(i-1),6],TS[i,6])*yfm,
		lty=1,lwd=1,col=1)
}
for(i in 1:length(TS[,2])){
	points(xos+TS[i,2]*xfm,2+(TS[i,6]*yfm),pch=20,cex=1,col=1)
}
#==================================================== Rate calculations
n=length(TS[,2]);n
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
    meandiff=TS[(i+k),6]-TS[i,6]				#mean diff.
    poolsd=PoolSD(TS[i+k,3],TS[i,3],TS[i+k,7],TS[i,7])
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
nc;idr[1:3,];idr[(nc-2):nc,]
idrx=idr[!rowSums(!is.finite(idr)),]  #remove rows that have -Inf
#========================================================== Write rate files
writefile<-"C://R_aaROEVchapt07//7.1.02_DudleyLambert2004_Zea_protein_out.csv"
write.csv(idrx,file=writefile,na="NA",row.names=FALSE)

proc.time()-ptm
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
bootresultd=TriPanelBC(idrx,"d",1,4.5,bootn,mode,psize,"lower")
wrlb.d=bootresultd[[1]];bootmat.d=bootresultd[[2]]
proc.time()-ptm
#========================================================== Plot LRI panel c
#send idrx matrix, n, mode(diff/rate), panel placement coordinates, mode
bootresultr=TriPanelBC(idrx,"r",13,4.5,bootn,mode,psize,"normal")
wrlb.r=bootresultr[[1]];bootmat.r=bootresultr[[2]]
proc.time()-ptm
#----------------------------------------------------------
text(7,-2.6,paste("Processing time: ",round(proc.time()[3]-ptm[3],digits=2),
  "seconds"),pos=4,cex=1,col=4)


