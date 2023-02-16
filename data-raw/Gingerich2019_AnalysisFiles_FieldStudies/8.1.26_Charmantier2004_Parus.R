##Philip D. Gingerich RATES OF EVOLUTION (Cambridge University Press 2019) 
##8.1.39_Geiger&al2018_Musmusculus
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
xfs=0;xfe=14;yfs=2.0;yfe=2.5			#x,y foreground start and end
xfm=(xoe-xos)/(xfe-xfs)				#x foreground multiplier
yfm=(yoe-yos)/(yfe-yfs)				#y foreground multiplier
lines(c(xos,xoe),c(yos,yos),lwd=1,lty=1,col=1)
lines(c(xos,xos),c(yos,yoe),lwd=1,lty=1,col=1)

#---------------------------------------------------------axis labels
for (i in seq(yos,yoe-.5,1)){
  lines(c(xos-.12,xos),c(i,i),lwd=1,lty=1,col=1)
}
for (i in seq(yos,yoe-.5,1)){
#	text(xos,i,format(round(.05*(i-7)+2,digits=2),nsmall=2),pos=2,cex=1,col=1)
  text(xos,i,format(round(.2*(i-9)+1.9,digits=2),nsmall=2),pos=2,cex=1,col=1)
}
for (i in xos+c(xfs:xfe)*1*xfm){
  lines(c(i,i),c(yos,yos-.1),lwd=1,lty=1,col=1)
}
for (i in c(1990,1992,1994,1996,1998,2000,2002)){
  text(xos+(i-1989+1)*xfm,yos,i,pos=1,cex=1,col=1)
}
text(xos,yos+8.0,expression(paste(italic('Parus caeruleus')*
	': tarsus length')),pos=4,cex=1.3,col=1)
text(xos-1.9,yos+4,'Ln length (mm)',
	srt=90,cex=1.3,col=1)
text(xos+7.5,yos-1.2,'Year',
	cex=1.3,col=1)
#============================================================ Load file
file1<-"C://R_aaROEVchapt08//8.1.26_Charmantier2004_Parus.csv"
C<-read.csv(file1,header=TRUE,sep=",",quote="\"",dec=".",fill=TRUE)
C[1:3,]
#----------------------------------------------------------------------
LPt<-matrix(nrow=length(na.omit(C[,2])),ncol=3)
  colnames(LPt)=c('Yr.','LnL_m','LnL_s')
for (i in 1:nrow(LPt)){
  LPt[i,1]=C[i+1,1]
  LPt[i,2]=C[i+1,2]
  LPt[i,3]=C[i+1,3]
};LPt
LRt<-matrix(nrow=length(na.omit(C[,4])),ncol=3)
  colnames(LRt)=c('Yr.','LnL_m','LnL_s')
for (i in 1:nrow(LRt)){
  LRt[i,1]=C[i+2,1]
  LRt[i,2]=C[i+2,4]
  LRt[i,3]=C[i+2,5]
};LRt
LMt<-matrix(nrow=length(na.omit(C[,6])),ncol=3)
  colnames(LMt)=c('Yr.','LnL_m','LnL_s')
for (i in 1:nrow(LMt)){
  LMt[i,1]=C[i+5,1]
  LMt[i,2]=C[i+5,6]
  LMt[i,3]=C[i+5,7]
};LMt
LPw<-matrix(nrow=length(na.omit(C[,2])),ncol=3)
  colnames(LPw)=c('Yr.','LnL_m','LnL_s')
for (i in 1:nrow(LPw)){
  LPw[i,1]=C[i+1,1]
  LPw[i,2]=C[i+1,8]
  LPw[i,3]=C[i+1,9]
};LPw
LRw<-matrix(nrow=length(na.omit(C[,4])),ncol=3)
  colnames(LRw)=c('Yr.','LnL_m','LnL_s')
for (i in 1:nrow(LRw)){
  LRw[i,1]=C[i+2,1]
  LRw[i,2]=C[i+2,10]
  LRw[i,3]=C[i+2,11]
};LRw
LMw<-matrix(nrow=length(na.omit(C[,6])),ncol=3)
  colnames(LMw)=c('Yr.','LnL_m','LnL_s')
for (i in 1:nrow(LMw)){
  LMw[i,1]=C[i+5,1]
  LMw[i,2]=C[i+5,12]
  LMw[i,3]=C[i+5,13]
};LMw

#------------------------------------------------Plot LPtarsus length
for (i in 1:nrow(LPt)){
  lines(xos+2-.04+(c(LPt[i,1],LPt[i,1])-1990)*xfm,
    40*c(LPt[i,2]-LPt[i,3],LPt[i,2]+LPt[i,3])-98,
    lty=1,lwd=1,col=gray(6/10))
}
for (i in 2:nrow(LPt)){
  lines(xos+2+(c(LPt[i-1,1],LPt[i,1])-1990)*xfm,40*c(LPt[i-1,2],LPt[i,2])-98,
	lty=1,lwd=1,col=1)
}
for (i in 1:nrow(LPt)){
  points(xos+2+(LPt[i,1]-1990)*xfm,40*LPt[i,2]-98,pch=22,bg='gray',col=1)
}
#------------------------------------------------Plot LRtarsus length
for (i in 1:nrow(LRt)){
  lines(xos+2+.04+(c(LRt[i,1],LRt[i,1])-1990)*xfm,
    40*c(LRt[i,2]-LRt[i,3],LRt[i,2]+LRt[i,3])-98,
    lty=1,lwd=1,col=gray(6/10))
}
for (i in 2:nrow(LRt)){
  lines(xos+2+(c(LRt[i-1,1],LRt[i,1])-1990)*xfm,40*c(LRt[i-1,2],LRt[i,2])-98,
	lty=1,lwd=1,col=1)
}
for (i in 1:nrow(LRt)){
  points(xos+2+(LRt[i,1]-1990)*xfm,40*LRt[i,2]-98,pch=24,bg='white',col=1)
}
#------------------------------------------------Plot LMtarsus length
for (i in 1:nrow(LMt)){
  lines(xos+2+(c(LMt[i,1],LMt[i,1])-1990)*xfm,
    40*c(LMt[i,2]-LMt[i,3],LMt[i,2]+LMt[i,3])-98,
    lty=1,lwd=1,col=gray(6/10))
}
for (i in 2:nrow(LMt)){
  lines(xos+2+(c(LMt[i-1,1],LMt[i,1])-1990)*xfm,40*c(LMt[i-1,2],LMt[i,2])-98,
	lty=1,lwd=1,col=1)
}
for (i in 1:nrow(LMt)){
  points(xos+2+(LMt[i,1]-1990)*xfm,40*LMt[i,2]-98,pch=23,bg=1,col=1)
}
#============================================= Rate calc: tarsus length
gentime=c(3.00,2.00,2.00)		#generation time in years (rounded)
idr=matrix(nrow=0,ncol=9)
  colnames(idr)=c('int','diff.sd','rate.sd.g','log.i','log|d|','log|r|',
    'sbn','wgt','fgen')		#fgen is fraction of a generation
for (run in 1:3){
  if(run==1){LT=LPt}
  if(run==2){LT=LRt}
  if(run==3){LT=LMt}
  n=nrow(LT);n
  nn=.5*(n-1)*n;nn
  id<-numeric(length=9)
  for (k in 1:(n-1)){		#run length
    for (i in 1:(n-k)){  		#starting position
      if(k<=gentime[run]){id[1]=1} else {id[1]=k/gentime[run]}	#interval
      meandiff=LT[(i+k),2]-LT[i,2]					#mean diff.
      poolsd=LT[i,3]							#n1,n2,sd1,sd2
      id[2]=meandiff/poolsd						#diff.sd
      id[3]=abs(id[2])/id[1]						#rate.sd.gen
      id[4]=log10(id[1])						#log.i
      id[5]=log10(abs(id[2]))						#log.d
      id[6]=log10(id[3])						#log.r
      if(k==1){id[7]=1}else{id[7]=3}				#sbn
      if(k==1){id[7]=1}else{id[7]=3}				#sbn
      id[8]=1/id[1]							#wgt
      id[9]=k/gentime[run]						#fract of gen.
      idr=rbind(idr,id)
    }  
  }
}
idr[1:3,]
idrx=idr[!rowSums(!is.finite(idr)),]  #remove rows that have -Inf
idrxT=idrx[,1:8]
idrc=idrx[idrx[,9]>=1,]  #keep only rows where interval >= gentime
idrcT=idrc[,1:8]
nrow(idr);nrow(idrx);nrow(idrxT);nrow(idrcT)
#writefile<-"C://R_aaROEVchapt08//8.1.26_Charmantier2004_Parus_out.csv"
#write.csv(idrxT,file=writefile,na="NA",row.names=FALSE) 
#============================================= Rate calc: weight
gentime=c(3.00,2.00,2.00)		#generation time in years (rounded)
idr=matrix(nrow=0,ncol=9)
  colnames(idr)=c('int','diff.sd','rate.sd.g','log.i','log|d|','log|r|',
    'sbn','wgt','fgen')		#fgen is fraction of a generation
for (run in 1:3){
  if(run==1){LT=LPw}
  if(run==2){LT=LRw}
  if(run==3){LT=LMw}
  n=nrow(LT);n
  nn=.5*(n-1)*n;nn
  id<-numeric(length=9)
  for (k in 1:(n-1)){		#run length
    for (i in 1:(n-k)){  		#starting position
      if(k<=gentime[run]){id[1]=1} else {id[1]=k/gentime[run]}	#interval
      meandiff=LT[(i+k),2]-LT[i,2]					#mean diff.
      poolsd=LT[i,3]							#n1,n2,sd1,sd2
      id[2]=meandiff/poolsd						#diff.sd
      id[3]=abs(id[2])/id[1]						#rate.sd.gen
      id[4]=log10(id[1])						#log.i
      id[5]=log10(abs(id[2]))						#log.d
      id[6]=log10(id[3])						#log.r
      if(k==1){id[7]=1}else{id[7]=3}				#sbn
      if(k==1){id[7]=1}else{id[7]=3}				#sbn
      id[8]=1/id[1]							#wgt
      id[9]=k/gentime[run]						#fract of gen.
      idr=rbind(idr,id)
    }  
  }
}
idr[1:3,]
idrx=idr[!rowSums(!is.finite(idr)),]  #remove rows that have -Inf
idrxW=idrx[,1:8]
idrc=idrx[idrx[,9]>=1,]  #keep only rows where interval >= gentime
idrcW=idrc[,1:8]
nrow(idr);nrow(idrx);nrow(idrxW);nrow(idrcW)
#writefile<-"C://R_aaROEVchapt08//8.1.26_Charmantier2004_Parus_out.csv"
#write.csv(idrxW,file=writefile,na="NA",row.names=FALSE) 
#========================================================== Write 'all' file
idrxA=rbind(idrxT,idrxW)
writefile<-"C://R_aaROEVchapt08//8.1.26_Charmantier2004_Parus_out.csv"
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
bootresultd=TriPanelBC(idrcT,"r",1,4.5,bootn,mode,psize,"lower")
wrlb.d=bootresultd[[1]];bootmat.d=bootresultd[[2]]
proc.time()-ptm
#========================================================== Plot LRI panel c
#send idrx matrix, n, mode(diff/rate), panel placement coordinates, mode
bootresultr=TriPanelBC(idrcW,"r",13,4.5,bootn,mode,psize,"lower")
wrlb.r=bootresultr[[1]];bootmat.r=bootresultr[[2]]
proc.time()-ptm
#========================================================== Label B and C
rect(-.5,6.7,21.5,7.5,border=NA,col='white')
rect(.1,6,9,6.7,border=NA,col='white')
rect(12.1,6,21,6.7,border=NA,col='white')
text(0,6.5,"Fledgling tarsus length",pos=4,cex=1.3,col=1)
text(12,6.5,"Fledgling body weight",pos=4,cex=1.3,col=1)

#----------------------------------------------------------
text(7,-2.6,paste("Processing time: ",round(proc.time()[3]-ptm[3],digits=2),
  "seconds"),pos=4,cex=1,col=4)


