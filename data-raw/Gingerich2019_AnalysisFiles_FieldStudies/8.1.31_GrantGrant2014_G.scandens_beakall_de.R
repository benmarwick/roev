##Philip D. Gingerich RATES OF EVOLUTION (Cambridge University Press 2019) 
##8.1.31_GrantGrant2014_G.scandens_beakall_de
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
xfs=0;xfe=40;yfs=2.25;yfe=2.75		#x,y foreground start and end
xfm=(xoe-xos)/(xfe-xfs)				#x foreground multiplier
yfm=(yoe-yos)/(yfe-yfs)				#y foreground multiplier
lines(c(xos,xoe),c(yos,yos),lwd=1,lty=1,col=1)
lines(c(xos,xos),c(yos,yoe),lwd=1,lty=1,col=1)

#----------------------------------------------------- y axis labels
for (i in seq(yos,yoe-.5,.8)){
  lines(c(xos-.12,xos),c(i,i),lwd=1,lty=1,col=1)
}
for (i in seq(yos,yoe-5.5,.8)){
  text(xos,i+6*.8,format(round(.5/8*(i-9)+2.55,digits=2),nsmall=2),pos=2,cex=1,col=1)
}
for (i in seq(yos,yoe-4.5,.8)){
  text(xos,i,format(round(.5/8*(i-9)+2.05,digits=2),nsmall=2),pos=2,cex=1,col=1)
}
#---------------------------------------------------- x axis labels
for (i in xos+c(xfs:xfe)*1*xfm){
  lines(c(i,i),c(yos,yos-.1),lwd=1,lty=1,col=1)
}
for (i in c(1975,1980,1985,1990,1995,2000,2005,2010)){
  text(xos+(i-1975+3)*xfm,yos,i,pos=1,cex=1,col=1)
}
text(xos,yos+8.0,expression(paste(italic('Geospiza scandens')*
	': beak length, depth, width')),pos=4,cex=1.3,col=1)
text(xos-1.9,yos+4,'Ln measurement (mm)',
	srt=90,cex=1.3,col=1)
text(xos+7.5,yos-1.2,'Year',
	cex=1.3,col=1)
rxll=xos+4.5*xfm;ryll=9.05;rxur=xos+5.5*xfm;ryur=16.6	#1977
rect(rxll,ryll,rxur,ryur,col=gray(9/10),border=NA)
rxll=xos+12.5*xfm;ryll=9.05;rxur=xos+13.5*xfm;ryur=16.6	#1985
rect(rxll,ryll,rxur,ryur,col=gray(9/10),border=NA)
rxll=xos+15.5*xfm;ryll=9.05;rxur=xos+17.5*xfm;ryur=16.6	#1988-89
rect(rxll,ryll,rxur,ryur,col=gray(9/10),border=NA)
rxll=xos+23.5*xfm;ryll=9.05;rxur=xos+24.5*xfm;ryur=16.6	#1996
rect(rxll,ryll,rxur,ryur,col=gray(9/10),border=NA)
rxll=xos+31.5*xfm;ryll=9.05;rxur=xos+32.5*xfm;ryur=16.6	#2004
rect(rxll,ryll,rxur,ryur,col=gray(9/10),border=NA)
#============================================================ Load file
file1<-"C://R_aaROEVchapt08//8.1.31_GrantGrant2014_G.scandens_beakall.csv"
G<-read.csv(file1,header=TRUE,sep=",",quote="\"",dec=".",fill=TRUE);G[1:3,]
GL<-matrix(nrow=40,ncol=8)
  colnames(GL)=c('Step','N','LnB.lx','LnB.lsd','LnB.dx','LnB.dsd',
	'LnB.wx','LnB.wsd')
for (i in 1:length(GL[,1])){
  GL[i,1]=G[i,1]-1972;GL[i,2]=G[i,3]
  GL[i,3]=log(G[i,4]);GL[i,4]=(.5*G[i,7]*sqrt(G[i,3]))/G[i,4]
  GL[i,5]=log(G[i,5]);GL[i,6]=(.5*G[i,8]*sqrt(G[i,3]))/G[i,5]
  GL[i,7]=log(G[i,6]);GL[i,8]=(.5*G[i,9]*sqrt(G[i,3]))/G[i,6]
}
GL[c(1:2,39:40),]
#------------------------------------------------Plot beak length
for (i in 1:40){
  lines(xos+c(GL[i,1],GL[i,1])*xfm,
    16*c(GL[i,3]-GL[i,4],GL[i,3]+GL[i,4])-27,
    lty=1,lwd=1,col=gray(6/10))
}
for (i in 2:40){
  lines(xos+c(GL[i-1,1],GL[i,1])*xfm,16*c(GL[i-1,3],GL[i,3])-27,
	lty=1,lwd=1,col=1)
}
for (i in 1:40){
  points(xos+GL[i,1]*xfm,16*GL[i,3]-27,pch=19,col=1)
}
#------------------------------------------------Plot beak depth
for (i in 1:40){
  lines(xos+.03+c(GL[i,1],GL[i,1])*xfm,
    16*c(GL[i,5]-GL[i,6],GL[i,5]+GL[i,6])-23.8,
    lty=1,lwd=1,col=gray(6/10))
}
for (i in 2:40){
  lines(xos+c(GL[i-1,1],GL[i,1])*xfm,16*c(GL[i-1,5],GL[i,5])-23.8,
	lty=1,lwd=1,col=1)
}
for (i in 1:40){
  points(xos+GL[i,1]*xfm,16*GL[i,5]-23.8,pch=19,col=1)
}
#------------------------------------------------Plot beak width
for (i in 1:40){
  lines(xos-.03+c(GL[i,1],GL[i,1])*xfm,
    16*c(GL[i,7]-GL[i,8],GL[i,7]+GL[i,8])-23.8,
    lty=1,lwd=1,col=gray(6/10))
}
for (i in 2:40){
  lines(xos+c(GL[i-1,1],GL[i,1])*xfm,16*c(GL[i-1,7],GL[i,7])-23.8,
	lty=1,lwd=1,col=1)
}
for (i in 1:40){
  points(xos+GL[i,1]*xfm,16*GL[i,7]-23.8,pch=19,col=1)
}
#------------------------------------------------Superimpose points
for (i in 1:40){
  points(xos+GL[i,1]*xfm,16*GL[i,3]-27,pch=19,col=1)
}
for (i in 1:40){
  points(xos+GL[i,1]*xfm,16*GL[i,5]-23.8,pch=19,col=1)
}
for (i in 1:40){
  points(xos+GL[i,1]*xfm,16*GL[i,7]-23.8,pch=19,col=1)
}
#============================================= Rate calc: beak length
gentime=5								#generation time
rect(3.6,9.2,9,9.6,col='white',border=NA)
text(3.5,9.4,"Generation time: 5 years",pos=4,cex=.9,col=1)
idrsum=0
GL[1:3,]
n=nrow(GL);n
nn=.5*(n-1)*n;nn
idr=matrix(nrow=nn,ncol=9)
  colnames(idr)=c('int','diff.sd','rate.sd.g','log.i','log|d|','log|r|',
    'sbn','wgt','fgen')		#fgen is fraction of a generation
nc=0
for (k in 1:(n-1)){	#run length
  for (i in 1:(n-k)){  			#starting position
    nc=nc+1
    if(k<=gentime){idr[nc,1]=1} else {idr[nc,1]=k/gentime}	#interval
    meandiff=GL[(i+k),3]-GL[i,3]					#mean diff.
    poolsd=PoolSD(GL[i+k,2],GL[i,2],GL[i+k,4],GL[i,4])	#n1,n2,sd1,sd2
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
#============================================= Rate calc: beak depth
GL[1:3,]
n=nrow(GL);n
nn=.5*(n-1)*n;nn
idr=matrix(nrow=nn,ncol=9)
  colnames(idr)=c('int','diff.sd','rate.sd.g','log.i','log|d|','log|r|',
    'sbn','wgt','fgen')		#fgen is fraction of a generation
nc=0
for (k in 1:(n-1)){    #run length
  for (i in 1:(n-k)){  #starting position
    nc=nc+1
    if(k<=gentime){idr[nc,1]=1} else {idr[nc,1]=k/gentime}	#interval
    meandiff=GL[(i+k),5]-GL[i,5]					#mean diff.
    poolsd=PoolSD(GL[i+k,2],GL[i,2],GL[i+k,6],GL[i,6])	#n1,n2,sd1,sd2
    idr[nc,2]=meandiff/poolsd						#diff.sd
    idr[nc,3]=abs(idr[nc,2])/idr[nc,1]				#rate.sd
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
idrxD=idrx[,1:8]
idrc=idrx[idrx[,9]>=1,]  #keep only rows where interval >= gentime
idrcD=idrc[,1:8]
nrow(idr);nrow(idrx);nrow(idrxD);nrow(idrcD)
#============================================= Rate calc: beak width
GL[1:3,]
n=nrow(GL);n
nn=.5*(n-1)*n;nn
idr=matrix(nrow=nn,ncol=9)
  colnames(idr)=c('int','diff.sd','rate.sd.g','log.i','log|d|','log|r|',
    'sbn','wgt','fgen')		#fgen is fraction of a generation
nc=0
for (k in 1:(n-1)){    #run length
  for (i in 1:(n-k)){  #starting position
    nc=nc+1
    if(k<=gentime){idr[nc,1]=1} else {idr[nc,1]=k/gentime}	#interval
    meandiff=GL[(i+k),7]-GL[i,7]					#mean diff.
    poolsd=PoolSD(GL[i+k,2],GL[i,2],GL[i+k,8],GL[i,8])	#n1,n2,sd1,sd2
    idr[nc,2]=meandiff/poolsd						#diff.sd
    idr[nc,3]=abs(idr[nc,2])/idr[nc,1]				#rate.sd
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
idrxW=idrx[,1:8]
idrc=idrx[idrx[,9]>=1,]  #keep only rows where interval >= gentime
idrcW=idrc[,1:8]
nrow(idr);nrow(idrx);nrow(idrxW);nrow(idrcW)
#========================================================== Write 'all' file
idrxA=rbind(idrxL,idrxD,idrxW)
#writefile<-"C://ROEVchapt08//8.1.31_GrantGrant2014_G.scandens_beakall_out.csv"
#write.csv(idrxA,file=writefile,na="NA",row.names=FALSE) 

idrsum;nrow(idrxA)
idrcA=rbind(idrcL,idrcD,idrcW)
#========================================================== Plot LRI panel d
#assign 'bootn' as boot number
bootn<-1000;text(-1.7,-2.6,paste("Boot n =",bootn),pos=4,cex=1,col=4)
#assign 'mode' as "medians","all","mixed"
mode<-"all";text(2,-2.6,paste("Mode: ",mode),pos=4,cex=1,col=4)
#assign circle size for points (1.5 or 2)
psize<-1.5	#2	
#assign 'equation' position as "normal","lower","none" at end of each call
#send (1)idrx matrix, (2)mode(diff/rate), (3)panel placement coordinate x, 
#  (4)panel placement coordinate y, (5)bootn, (6)mode, (7)psize, (8)equation
bootresultd=TriPanelBC(idrcW,"r",1,4.5,bootn,mode,psize,"lower")
wrlb.d=bootresultd[[1]];bootmat.d=bootresultd[[2]]
proc.time()-ptm
#========================================================== Plot LRI panel e
#send idrx matrix, n, mode(diff/rate), panel placement coordinates, mode
bootresultr=TriPanelBC(idrcA,"r",13,4.5,bootn,mode,psize,"lower")
wrlb.r=bootresultr[[1]];bootmat.r=bootresultr[[2]]
proc.time()-ptm
#========================================================== Label D and E
rect(-.5,6.7,21.5,7.5,border=NA,col='white')
rect(.1,6,9,6.7,border=NA,col='white')
rect(12.1,6,21,6.7,border=NA,col='white')
text(0,6.5,"Beak width",pos=4,cex=1.3,col=1)
text(12,6.5,"All rates",pos=4,cex=1.3,col=1)

#----------------------------------------------------------
text(7,-2.6,paste("Processing time: ",round(proc.time()[3]-ptm[3],digits=2),
  "seconds"),pos=4,cex=1,col=4)


