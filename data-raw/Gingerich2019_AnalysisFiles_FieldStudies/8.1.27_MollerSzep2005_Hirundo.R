##Philip D. Gingerich RATES OF EVOLUTION (Cambridge University Press 2019) 
##8.1.27_MollerSzep2005_Hirundo
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
xos=3.5;xoe=17;yos=9;yoe=17			#x,y original start and end
xfs=-1;xfe=19;yfs=2.8;yfe=4.0			#x,y foreground start and end
xfm=(xoe-xos)/(xfe-xfs)				#x foreground multiplier
yfm=(yoe-yos)/(yfe-yfs)				#y foreground multiplier
lines(c(xos,xoe),c(yos,yos),lwd=1,lty=1,col=1)
lines(c(xos,xos),c(yos,yoe),lwd=1,lty=1,col=1)

#---------------------------------------------------------axis labels
for (i in seq(yos,yoe-.5,1)){
  lines(c(xos-.12,xos),c(i,i),lwd=1,lty=1,col=1)
}
for (i in seq(yos,yoe-1,2)){
  text(xos,i,format(round(.05*(i-7)+4.3,digits=2),nsmall=2),pos=2,cex=1,col=1)
}
for (i in xos+c(xfs:xfe)*1*xfm){
  lines(c(i,i),c(yos,yos-.1),lwd=1,lty=1,col=1)
}
for (i in c(1984,1987,1990,1993,1996,1999,2002)){
  text(xos+(i-1984+1)*xfm,yos,i,pos=1,cex=1,col=1)
}
text(xos,yos+8.0,expression(paste(italic('Hirundo rustica')*
	': tail length')),pos=4,cex=1.3,col=1)
text(xos-1.9,yos+4,'Ln length (mm)',
	srt=90,cex=1.3,col=1)
text(xos+7,yos-1.2,'Year',
	cex=1.3,col=1)
#============================================================ Load file
file1<-"C://R_aaROEVchapt08//8.1.27_MollerSzep2005_Hirundo.csv"
M<-read.csv(file1,header=TRUE,sep=",",quote="\"",dec=".",fill=TRUE);M
#----------------------------------------------------------------------
ML<-matrix(nrow=nrow(M),ncol=7)
  colnames(ML)=c('Year','Nf','Ln.xf','Ln.sdf','Nm','Ln.xm','Ln.sdm')
for (i in 1:nrow(ML)){
  ML[i,1]=i						#interval
  ML[i,2]=M[i,5]					#female number
  ML[i,3]=log(M[i,6])				#female ln mean
  ML[i,4]=M[i,7]/M[i,6]				#female ln stdev
  ML[i,5]=M[i,2]					#male number
  ML[i,6]=log(M[i,3])				#male ln mean
  ML[i,7]=M[i,4]/M[i,3]				#male ln stdev
}
ML[c(1:3,(nrow(ML)-2):nrow(ML)),]
#------------------------------------------------Plot ln male feather lengths
for (i in 1:nrow(ML)){
  lines(xos+c(ML[i,1]-.03,ML[i,1]-.03)*xfm,
    20*c(ML[i,6]-ML[i,7],ML[i,6]+ML[i,7])-79,
    lty=1,lwd=1,col=gray(6/10))
}
for (i in 2:nrow(ML)){
  lines(xos+c(ML[i-1,1],ML[i,1])*xfm,20*c(ML[i-1,6],ML[i,6])-79,
	lty=1,lwd=1,col=1)
}
for (i in 1:nrow(ML)){
  points(xos+ML[i,1]*xfm,20*ML[i,6]-79,pch=19,cex=1.1,col=1)
}
#-------------------------------------------Plot ln female feather lengths
for (i in 1:nrow(ML)){
  lines(xos+c(ML[i,1]+.03,ML[i,1]+.03)*xfm,
    20*c(ML[i,3]-ML[i,4],ML[i,3]+ML[i,4])-79,
    lty=1,lwd=1,col=gray(6/10))
}
for (i in 2:nrow(ML)){
  lines(xos+c(ML[i-1,1],ML[i,1])*xfm,20*c(ML[i-1,3],ML[i,3])-79,
	lty=1,lwd=1,col=1)
}
for (i in 1:nrow(ML)){
  points(xos+ML[i,1]*xfm,20*ML[i,3]-79,pch=21,bg='white',cex=1.2,col=1)
}
#============================================ Rate calc female feather length
gentime=2					#generation time fide Moller-Szek 2005
rect(3.6,9.2,9,9.6,col='white',border=NA)
text(3.5,9.4,"Generation time: 2 years",pos=4,cex=.9,col=1)
idrsum=0
ML[1:3,]
n=nrow(ML);n
nn=.5*(n-1)*n;nn
idr=matrix(nrow=nn,ncol=9)
  colnames(idr)=c('int','diff.sd','rate.sd.g','log.i','log|d|','log|r|',
    'sbn','wgt','fgen')		#fgen is fraction of a generation
nc=0
for (k in 1:(n-1)){	#run length
  for (i in 1:(n-k)){  						#starting position
    nc=nc+1
    if(k<=gentime){idr[nc,1]=1} else {idr[nc,1]=k/gentime}	#interval
    meandiff=ML[(i+k),3]-ML[i,3]					#mean diff.
    poolsd=PoolSD(ML[i+k,2],ML[i,2],ML[i+k,4],ML[i,4])	#n1,n2,sd1,sd2
    idr[nc,2]=meandiff/poolsd						#diff.sd
    idr[nc,3]=abs(idr[nc,2])/idr[nc,1]				#rate.sd.gen
    idr[nc,4]=log10(idr[nc,1])					#log.i
    idr[nc,5]=log10(abs(idr[nc,2]))					#log.d
    idr[nc,6]=log10(idr[nc,3])					#log.r
    if(k==1){idr[nc,7]=1}else{idr[nc,7]=3}			#sbn
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
writefile<-"C://R_aaRateBook//RB_Chapter08//2005_MollerSzek_femtail_out.csv"
write.csv(idrxF,file=writefile,na="NA",row.names=FALSE) 
#============================================ Rate calc male feather length
n=nrow(ML);n
nn=.5*(n-1)*n;nn
idr=matrix(nrow=nn,ncol=9)
  colnames(idr)=c('int','diff.sd','rate.sd.g','log.i','log|d|','log|r|',
    'sbn','wgt','fgen')		#fgen is fraction of a generation
nc=0
for (k in 1:(n-1)){	#run length
  for (i in 1:(n-k)){  						#starting position
    nc=nc+1
    if(k<=gentime){idr[nc,1]=1} else {idr[nc,1]=k/gentime}	#interval
    meandiff=ML[(i+k),6]-ML[i,6]					#mean diff.
    poolsd=PoolSD(ML[i+k,5],ML[i,5],ML[i+k,7],ML[i,7])	#n1,n2,sd1,sd2
    idr[nc,2]=meandiff/poolsd						#diff.sd
    idr[nc,3]=abs(idr[nc,2])/idr[nc,1]				#rate.sd.gen
    idr[nc,4]=log10(idr[nc,1])					#log.i
    idr[nc,5]=log10(abs(idr[nc,2]))					#log.d
    idr[nc,6]=log10(idr[nc,3])					#log.r
    if(k==1){idr[nc,7]=1}else{idr[nc,7]=3}			#sbn
    idr[nc,8]=1/idr[nc,1]						#wgt
    idr[nc,9]=k/gentime							#fract of gen.
    #stdev[i]=poolsd
  }
}
idr[1:3,]
idrsum=idrsum+nrow(idr)
idrx=idr[!rowSums(!is.finite(idr)),]  #remove rows that have -Inf
idrxM=idrx[,1:8]
idrc=idrx[idrx[,9]>=1,]  #keep only rows where interval >= gentime
idrcM=idrc[,1:8]
nrow(idr);nrow(idrx);nrow(idrxM);nrow(idrcM)
writefile<-"C://R_aaRateBook//RB_Chapter08//2005_MollerSzek_maletail_out.csv"
write.csv(idrxM,file=writefile,na="NA",row.names=FALSE) 
#========================================================== Write 'all' file
idrxA=rbind(idrxF,idrxM)
writefile<-"C://R_aaROEVchapt08//8.1.27_MollerSzep2005_Hirundo_out.csv"
write.csv(idrxA,file=writefile,na="NA",row.names=FALSE) 

idrsum;nrow(idrxA)
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
bootresultd=TriPanelBC(idrcF,"r",1,4.5,bootn,mode,psize,"lower")
wrlb.d=bootresultd[[1]];bootmat.d=bootresultd[[2]]
proc.time()-ptm
#========================================================== Plot LRI panel c
#send idrx matrix, n, mode(diff/rate), panel placement coordinates, mode
bootresultr=TriPanelBC(idrcM,"r",13,4.5,bootn,mode,psize,"lower")
wrlb.r=bootresultr[[1]];bootmat.r=bootresultr[[2]]
proc.time()-ptm
#========================================================== Label B and C
rect(-.5,6.7,21.5,7.5,border=NA,col='white')
rect(.1,6,9,6.7,border=NA,col='white')
rect(12.1,6,21,6.7,border=NA,col='white')
text(0,6.5,"Female tail length",pos=4,cex=1.3,col=1)
text(12,6.5,"Male tail length",pos=4,cex=1.3,col=1)

#----------------------------------------------------------
text(7,-2.6,paste("Processing time: ",round(proc.time()[3]-ptm[3],digits=2),
  "seconds"),pos=4,cex=1,col=4)


