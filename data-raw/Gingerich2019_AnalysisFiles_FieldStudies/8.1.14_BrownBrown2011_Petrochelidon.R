##Philip D. Gingerich RATES OF EVOLUTION (Cambridge University Press 2019)
##8.1.14_BrownBrown2011_Petrochelidon
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
xos=3.5;xoe=17;yos=9;yoe=17			#x,y original start and end
xfs=0;xfe=11;yfs=5;yfe=5.4			#x,y foreground start and end
xfm=(xoe-xos)/(xfe-xfs)				#x foreground multiplier
yfm=(yoe-yos)/(yfe-yfs)				#y foreground multiplier
lines(c(xos,xoe),c(yos,yos),lwd=1,lty=1,col=1)
lines(c(xos,xos),c(yos,yoe),lwd=1,lty=1,col=1)
#---------------------------------------------------------axis labels
for (i in seq(yos,yoe-.5,1)){
  lines(c(xos-.12,xos),c(i,i),lwd=1,lty=1,col=1)
}
for (i in seq(yos+6,yoe-.5,1)){
  text(xos,i,format(round(.05*(i-7)+4.3-.05,digits=2),nsmall=2),pos=2,cex=1,col=1)
}
for (i in seq(yos+1,yoe-4,1)){
  text(xos,i,format(round(.05*(i-7)+2.3-.1,digits=2),nsmall=2),pos=2,cex=1,col=1)
}
#-----------------------------------------------
for (i in xos+c(xfs:xfe)*1*xfm){
  lines(c(i,i),c(yos,yos-.1),lwd=1,lty=1,col=1)
}
for (i in seq(1996,2006,2)){
  text(xos+(i-1996+1)*xfm,yos,i,pos=1,cex=1,col=1)
}
text(xos,yos+8.0,expression(paste(italic('Petrochelidon pyrrhonota')*
	': wing and tarsus length')),pos=4,cex=1.3,col=1)
text(xos-1.9,yos+4,'Ln length (mm)',
	srt=90,cex=1.3,col=1)
text(xos+7,yos-1.2,'Year',
	cex=1.3,col=1)
#============================================================ Load file
file1<-"../Gingerich2019_AnalysisFiles_FieldStudies/8.1.14_BrownBrown2011_Petrochelidon.csv"
B<-read.csv(file1,header=TRUE,sep=",",quote="\"",dec=".",fill=TRUE)
B
#======================================================================
idrxAll=matrix(nrow=0,ncol=8)
idrsum=0
for (Fig in 1:5){
  if (Fig==1){
    fig='1A';S=1;E=11
    idrc_tarslen=matrix(nrow=0,ncol=8)
  }
  if (Fig==2){
    fig='1B';S=12;E=22
    idrc_billlen=matrix(nrow=0,ncol=8)
  }
  if (Fig==3){
    fig='1C';S=23;E=33
    idrc_billwid=matrix(nrow=0,ncol=8)
  }
  if (Fig==4){
    fig='2A';S=34;E=44
    idrc_winglen=matrix(nrow=0,ncol=8)
  }
  if (Fig==5){
    fig='2B';S=45;E=55
    idrc_taillen=matrix(nrow=0,ncol=8)
  }
  #---------------------------------------------------Female and male matrix
  B[1:3,]
  LFM<-matrix(nrow=E-S+1,ncol=7)
    colnames(LFM)=c('Yr','Nf','LMf','LSf','Nm','LMm','LSm')
  R=E-S+1;SL=S-1;c(S,E,SL)					#Run length; start line
  for (i in 1:R){
    LFM[i,1]=B[SL+i,3]						#Year
    LFM[i,2]=B[SL+i,4]						#N female
    LFM[i,3]=log(B[SL+i,5])					#Mean female
    LFM[i,4]=(sqrt(B[SL+i,4])*B[SL+i,6])/B[SL+i,5]	#Stdev/mean female
    LFM[i,5]=B[SL+i,7]						#N male
    LFM[i,6]=log(B[SL+i,8])					#Mean male
    LFM[i,7]=(sqrt(B[SL+i,7])*B[SL+i,9])/B[SL+i,8]	#Stdev/mean male
  };LFM
  #=======================================================Plot measurements
  scale=20
  if (Fig==4){
    for (i in 1:(E-S+1)){			#Male wing lengths
      lines(xos+.03+c(i,i)*xfm,
        scale*c(LFM[i,6]-LFM[i,7],LFM[i,6]+LFM[i,7])-78,
        lty=1,lwd=1,col=gray(6/10))
    }
    for (i in 2:(E-S+1)){
      lines(xos+c((i-1),i)*xfm,scale*c(LFM[i-1,6],LFM[i,6])-78,
    	lty=1,lwd=1,col=1)
    }
    for (i in 1:(E-S+1)){
      points(xos+i*xfm,scale*LFM[i,6]-78,pch=19,cex=1.1,col=1)
    }
    for (i in 1:(E-S+1)){			#Female wing lengths
      lines(xos-.03+c(i,i)*xfm,
        scale*c(LFM[i,3]-LFM[i,4],LFM[i,3]+LFM[i,4])-78,
        lty=1,lwd=1,col=gray(6/10))
    }
    for (i in 2:(E-S+1)){
      lines(xos+c((i-1),i)*xfm,scale*c(LFM[i-1,3],LFM[i,3])-78,
  	  lty=1,lwd=1,col=1)
    }
    for (i in 1:(E-S+1)){
      points(xos+i*xfm,scale*LFM[i,3]-78,pch=21,cex=1.2,bg='white',col=1)
    }
  }
  #-------------------------------------------------------------------
  if (Fig==1){
   for (i in 1:(E-S+1)){			#Male tarsus length
      lines(xos+.03+c(i,i)*xfm,
        scale*c(LFM[i,6]-LFM[i,7],LFM[i,6]+LFM[i,7])-37,
        lty=1,lwd=1,col=gray(6/10))
    }
    for (i in 2:(E-S+1)){
      lines(xos+c((i-1),i)*xfm,scale*c(LFM[i-1,6],LFM[i,6])-37,
    	lty=1,lwd=1,col=1)
    }
    for (i in 1:(E-S+1)){
      points(xos+i*xfm,scale*LFM[i,6]-37,pch=19,cex=1.1,col=1)
    }
    for (i in 1:(E-S+1)){			#Female tarsus length
      lines(xos-.03+c(i,i)*xfm,
        scale*c(LFM[i,3]-LFM[i,4],LFM[i,3]+LFM[i,4])-37,
        lty=1,lwd=1,col=gray(6/10))
    }
    for (i in 2:(E-S+1)){
      lines(xos+c((i-1),i)*xfm,scale*c(LFM[i-1,3],LFM[i,3])-37,
  	  lty=1,lwd=1,col=1)
    }
    for (i in 1:(E-S+1)){
      points(xos+i*xfm,scale*LFM[i,3]-37,pch=21,cex=1.2,bg='white',col=1)
    }
  }
  #============================================= Rate calc: female weights
  gentime=3					#estimate from BrownBrown 1998b table 10
  rect(3.6,9.2,9,9.6,col='white',border=NA)
  text(3.5,9.4,"Generation time: 3 years",pos=4,cex=.9,col=1)
  LFM[1:3,]
  n=nrow(LFM);n
  nn=.5*(n-1)*n;nn
  idr=matrix(nrow=nn,ncol=9)
    colnames(idr)=c('int','diff.sd','rate.sd.g','log.i','log|d|','log|r|',
      'sbn','wgt','fgen')		#fgen is fraction of a generation
  nc=0
  for (k in 1:(n-1)){	#run length
    for (i in 1:(n-k)){  						#starting position
      nc=nc+1
    	if(k<=gentime){idr[nc,1]=1} else {idr[nc,1]=k/gentime}#interval
      meandiff=LFM[(i+k),3]-LFM[i,3]				#mean diff.
      poolsd=PoolSD(LFM[i+k,2],LFM[i,2],LFM[i+k,4],LFM[i,4])#n1,n2,sd1,sd2
      idr[nc,2]=meandiff/poolsd					#diff.sd
      idr[nc,3]=abs(idr[nc,2])/idr[nc,1]				#rate.sd.gen
      idr[nc,4]=log10(idr[nc,1])					#log.i
      idr[nc,5]=log10(abs(idr[nc,2]))				#log.d
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
  #============================================= Rate calc: male weights
  LFM[1:3,]
  n=nrow(LFM);n
  nn=.5*(n-1)*n;nn
  idr=matrix(nrow=nn,ncol=9)
    colnames(idr)=c('int','diff.sd','rate.sd.g','log.i','log|d|','log|r|',
      'sbn','wgt','fgen')		#fgen is fraction of a generation
  nc=0
  for (k in 1:(n-1)){    #run length
    for (i in 1:(n-k)){  #starting position
      nc=nc+1
      if(k<=gentime){idr[nc,1]=1} else {idr[nc,1]=k/gentime}#interval
      meandiff=LFM[(i+k),6]-LFM[i,6]				#mean diff.
      poolsd=PoolSD(LFM[i+k,5],LFM[i,5],LFM[i+k,7],LFM[i,7])#n1,n2,sd1,sd2
      idr[nc,2]=meandiff/poolsd					#diff.sd
      idr[nc,3]=abs(idr[nc,2])/idr[nc,1]				#rate.sd
      idr[nc,4]=log10(idr[nc,1])					#log.i
      idr[nc,5]=log10(abs(idr[nc,2]))				#log.d
      idr[nc,6]=log10(idr[nc,3])					#log.r
      if(k==1){idr[nc,7]=1}else{idr[nc,7]=3}			#sbn
      idr[nc,8]=1/idr[nc,1]						#wgt
      idr[nc,9]=k/gentime						#fract of gen.
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
  #======================================================== Write 'all' file
  idrxAll=rbind(idrxAll,idrxF,idrxM)
  if (Fig==1){idrc_tarslen=rbind(idrcF,idrcM)}
  if (Fig==2){idrc_billlen=rbind(idrcF,idrcM)}
  if (Fig==3){idrc_billwid=rbind(idrcF,idrcM)}
  if (Fig==4){idrc_winglen=rbind(idrcF,idrcM)}
  if (Fig==5){idrc_taillen=rbind(idrcF,idrcM)}
}
# writefile<-"C://R_aaROEVchapt08//8.1.14_BrownBrown2011_Petrochelidon_out.csv"
# write.csv(idrxAll,file=writefile,na="NA",row.names=FALSE)

idrsum;nrow(idrxAll)
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
bootresultd=TriPanelBC(idrc_winglen,"r",1,4.5,bootn,mode,psize,"lower")
wrlb.d=bootresultd[[1]];bootmat.d=bootresultd[[2]]
proc.time()-ptm
#========================================================== Plot LRI panel c
#send idrx matrix, n, mode(diff/rate), panel placement coordinates, mode
bootresultr=TriPanelBC(idrc_tarslen,"r",13,4.5,bootn,mode,psize,"lower")
wrlb.r=bootresultr[[1]];bootmat.r=bootresultr[[2]]
proc.time()-ptm
#========================================================== Label B and C
rect(-.5,6.7,21.5,7.5,border=NA,col='white')
rect(.1,6,9,6.7,border=NA,col='white')
rect(12.1,6,21,6.7,border=NA,col='white')
text(0,6.5,"Wing length",pos=4,cex=1.3,col=1)
text(12,6.5,"Tarsus length",pos=4,cex=1.3,col=1)

#----------------------------------------------------------
text(7,-2.6,paste("Processing time: ",round(proc.time()[3]-ptm[3],digits=2),
  "seconds"),pos=4,cex=1,col=4)
warnings()
