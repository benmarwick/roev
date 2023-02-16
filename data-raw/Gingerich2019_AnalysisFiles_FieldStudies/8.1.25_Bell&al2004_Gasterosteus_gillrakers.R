##Philip D. Gingerich RATES OF EVOLUTION (Cambridge University Press 2019) 
##8.1.25_Bell&al2004_Gasterosteus_gillrakers
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
	xaxp=c(xr[1],xr[2],xr[1]+xr[2]),
	yaxp=c(yr[1],yr[2],abs(yr[1]-yr[2])),		#y-axis extreme ticks and number of ticks
	cex.axis=.8,
	cex.lab=1.1)			#label size
dev.size('in')
#======================================================= Plot panel a
xos=3.5;xoe=18;yos=9;yoe=17			#x,y original start and end
xfs=0;xfe=28;yfs=2.9;yfe=3.3			#x,y foreground start and end
xfm=(xoe-xos)/(xfe-xfs)				#x foreground multiplier
yfm=(yoe-yos)/(yfe-yfs)				#y foreground multiplier
lines(c(xos,xoe),c(yos,yos),lwd=1,lty=1,col=1)
lines(c(xos,xos),c(yos,yoe),lwd=1,lty=1,col=1)

#---------------------------------------------------------axis labels
for (i in seq(yos,yoe-.05,1)){
  lines(c(xos-.12,xos),c(i,i),lwd=1,lty=1,col=1)
}
for (i in seq(yos,yoe-.05,1)){
  text(xos,i,format(round(.05*(i-7)+2.8,digits=2),nsmall=2),pos=2,cex=1,col=1)
}
#-----------------------------------------------
for (i in xos+c(xfs:xfe)*1*xfm){
  lines(c(i,i),c(yos,yos-.1),lwd=1,lty=1,col=1)
}
for (i in c(1990,1995,2000,2005,2010,2015)){
  text(xos+(i-1990+1)*xfm,yos,i,pos=1,cex=1,col=1)
}
#-----------------------------------------------
text(xos,yos+8.0,expression(paste(italic('Gasterosteus aculeatus')*
	': gill raker count')),pos=4,cex=1.3,col=1)
text(xos-1.9,yos+4,'Ln count',
	srt=90,cex=1.3,col=1)
text(xos+7.5,yos-1.2,'Year',
	cex=1.3,col=1)

#============================================================ Load file
file1<-"C://R_aaROEVchapt08//8.1.25_Bell&al2004_Gasterosteus_gillrakers.csv"
B<-read.csv(file1,header=TRUE,sep=",",quote="\"",dec=".",fill=TRUE);B[1:3,]
LB<-matrix(nrow=nrow(B),ncol=4)
  colnames(LB)=c('Gen.','N','LnB.x','LnB.sd')
for (i in 1:nrow(B)){
  LB[i,1]=i
  LB[i,2]=B[i,2]
  LB[i,3]=log(B[i,3])
  LB[i,4]=(sqrt(LB[i,2])*B[i,4])/B[i,3]
}
LB[c(1:2,nrow(LB)-1,nrow(LB)),]
#------------------------------------------------Plot gill raker time series
for (i in seq(yos,yoe-.5,1)){
  lines(xos+c(0,nrow(LB))*xfm,c(i,i),lty=2,lwd=1,col=gray(6/10))
}
for (i in 1:nrow(LB)){
  lines(xos-.04+c(LB[i,1],LB[i,1])*xfm,
    20*c(LB[i,3]-LB[i,4],LB[i,3]+LB[i,4])-49,
    lty=1,lwd=1,col=gray(6/10))
}
for (i in 2:nrow(LB)){
  lines(xos+c(LB[i-1,1],LB[i,1])*xfm,20*c(LB[i-1,3],LB[i,3])-49,
	lty=1,lwd=1,col=1)
}
for (i in 1:nrow(LB)){
  points(xos+LB[i,1]*xfm,20*LB[i,3]-49,pch=19,col=1)
}
#============================================= Rate calc: gill rakers
gentime=2								#generation time
rect(3.6,9.2,9,9.6,col='white',border=NA)
text(3.5,9.4,"Generation time: 2 years",pos=4,cex=.9,col=1)
idrsum=0
LB[1:3,]
n=nrow(LB);n
nn=.5*(n-1)*n;nn
idr=matrix(nrow=nn,ncol=9)
  colnames(idr)=c('int','diff.sd','rate.sd.g','log.i','log|d|','log|r|',
    'sbn','wgt','fgen')		#fgen is fraction of a generation
nc=0
for (k in 1:(n-1)){		#run length
  for (i in 1:(n-k)){  		#starting position
    nc=nc+1
    if(k<=gentime){idr[nc,1]=1} else {idr[nc,1]=k/gentime}	#interval
    meandiff=LB[(i+k),3]-LB[i,3]					#mean diff.
    poolsd=PoolSD(LB[i+k,2],LB[i,2],LB[i+k,4],LB[i,4])	#n1,n2,sd1,sd2
    idr[nc,2]=meandiff/poolsd						#diff.sd
    idr[nc,3]=abs(idr[nc,2])/idr[nc,1]				#rate.sd.gen
    idr[nc,4]=log10(idr[nc,1])					#log.i
    idr[nc,5]=log10(abs(idr[nc,2]))					#log.d
    idr[nc,6]=log10(idr[nc,3])					#log.r
    if(k==1){idr[nc,7]=1}else{idr[nc,7]=3}			#sbn
    if(k==1){idr[nc,7]=1}else{idr[nc,7]=3}			#sbn
    idr[nc,8]=1/idr[nc,1]						#wgt
    idr[nc,9]=k/gentime							#fract of gen.
  }
}
idr[1:3,]
idrsum=idrsum+nrow(idr)
idrx=idr[!rowSums(!is.finite(idr)),]  #remove rows that have -Inf
idrxL=idrx[,1:8]
idrc=idrx[idrx[,9]>=1,]  #keep only rows where interval >= gentime
idrcL=idrc[,1:8]
nrow(idr);nrow(idrx);nrow(idrxL);nrow(idrcL)
writefile<-"C://R_aaROEVchapt08//8.1.25_Bell&al2004_Gasterosteus_gillrakers_out.csv"
write.csv(idrxL,file=writefile,na="NA",row.names=FALSE) 

#========================================================== Plot LRI panel b
#assign 'bootn' as boot number
bootn<-1000;text(-1.7,-2.6,paste("Boot n =",bootn),pos=4,cex=1,col=4)
#assign 'mode' as "medians","all","mixed"
mode<-"all";text(2,-2.6,paste("Mode: ",mode),pos=4,cex=1,col=4)
#assign circle size for points (1.5 or 2)
psize<-1.5	
#assign 'equation' position as "normal","lower","none" at end of each call
#send (1)idrx matrix, (2)mode(diff/rate), (3)panel placement coordinate x, 
#  (4)panel placement coordinate y, (5)bootn, (6)mode, (7)psize, (8)equation
bootresultd=TriPanelBC(idrcL,"d",1,4.5,bootn,mode,psize,"lower")
wrlb.d=bootresultd[[1]];bootmat.d=bootresultd[[2]]
proc.time()-ptm
#========================================================== Plot LRI panel c
#send idrx matrix, n, mode(diff/rate), panel placement coordinates, mode
bootresultr=TriPanelBC(idrcL,"r",13,4.5,bootn,mode,psize,"lower")
wrlb.r=bootresultr[[1]];bootmat.r=bootresultr[[2]]
proc.time()-ptm
#========================================================== Label B and C
rect(-.5,6.7,21.5,7.5,border=NA,col='white')
rect(.1,6,9,6.7,border=NA,col='white')
rect(12.1,6,21,6.7,border=NA,col='white')
text(0,6.5,"Gill raker count",pos=4,cex=1.3,col=1)
text(12,6.5,"Gill raker count",pos=4,cex=1.3,col=1)

#----------------------------------------------------------
text(7,-2.6,paste("Processing time: ",round(proc.time()[3]-ptm[3],digits=2),
  "seconds"),pos=4,cex=1,col=4)


