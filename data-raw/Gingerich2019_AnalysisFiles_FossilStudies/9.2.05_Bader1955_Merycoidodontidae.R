##Philip D. Gingerich RATES OF EVOLUTION (Cambridge University Press 2019)
##9.2.05_Bader1955_Merycoidodontidae
#==========================================================================

ptm<-proc.time()
##======================================================================##

#=====================================================================Setup
library(RATES);library(MASS)
#----------------------------------------------------------------------Plot

xr<-c(0,20);yr<-c(-2,18)		#xrange;yrange for plot axes
plot(xr,yr,					#set up plot
  xaxt ='n',
  yaxt ='n',
  xlab = "",
  ylab = "",
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

#======================================================= Plot panel a
xos=3;xoe=18;yos=9.5;yoe=17.5			#x,y original start and end
xfs=3.4;xfe=4.9;yfs=-22;yfe=-15		#x,y foreground start and end
xfm=(xoe-xos)/(xfe-xfs)				#x foreground multiplier
yfm=(yoe-yos)/(yfe-yfs)				#y foreground multiplier
lines(c(xos,xoe),c(yos,yos),lwd=1,lty=1,col=1)
lines(c(xos,xos),c(yos,yoe),lwd=1,lty=1,col=1)

#---------------------------------------------------------axis labels
for (i in seq(yos,yoe-.5,1)){				#y-axis
  lines(c(xos-.15,xos),c(i,i),lwd=1,lty=1,col=1)
}
for (i in seq(yos,yoe-.5,1)){
  text(xos,i,31.5-i,pos=2,cex=1,col=1)
}
for (i in xos+10*(seq(xfs,xfe,.1)-xfs)){		#x-axis
  lines(c(i,i),c(yos,yos-.15),lwd=1,lty=1,col=1)
}
for (i in xos+10*(seq(xfs,xfe,.1)-xfs)){
  text(i,yos,format(round(3.1+.1*i,digits=1),nsmall=1),
    pos=1,cex=1,col=1)
}
text(xos,yos+8.0,expression(paste('Merycoidodontidae'*
	': length of upper molar row')),pos=4,cex=1.3,col=1)
text(xos-1.5,yos+4,'Geological age (Ma)',
	srt=90,cex=1.3,col=1)
text(xos+7.5,yos-1.2,'Ln measurement (mm)',
	cex=1.3,col=1)
#============================================================ Load file
file1<-"../Gingerich2019_AnalysisFiles_FossilStudies/9.2.05_Bader1955_Merycoidodontidae.csv"
B<-read.csv(file1,header=TRUE,sep=",",quote="\"",dec=".",fill=TRUE)
B[1:5,]
#-------------------------------------------------Add histograms
B[,c(1,22)]
taxon<-vector(length=8)			#name vector
  taxon[1]='B. siouense';taxon[2]='B. wilsoni';taxon[3]='M. proprius'
  taxon[4]='M. matthewi';taxon[5]='M. relictus';taxon[6]='M. minimus'
  taxon[7]='M. elegans';taxon[8]='M. crabilli'
UMR<-matrix(nrow=10,ncol=4)		#upper molar row matrix
  colnames(UMR)=c('age','N','lnM','LnS')
for (i in 1:10){
  UMR[i,1]=B[3*i-2,2]
  UMR[i,2]=B[3*i-2,22]
  UMR[i,3]=log(B[3*i-1,22])
  UMR[i,4]=.01*B[3*i,22]
};UMR
UMR=UMR[c(1,3:7,9:10),]	#remove species buwaldi and arenarum
for (i in 5:7){
  lines(xos+10*(c(UMR[i,3],UMR[i+1,3])-xfs),
    yos+22-c(UMR[i,1],UMR[i+1,1]),
    lty=1,lwd=2,col=1)
}
for (i in 1:3){
  lines(xos+10*(c(UMR[i,3],UMR[i+1,3])-xfs),
    yos+22-c(UMR[i,1],UMR[i+1,1]),
    lty=1,lwd=2,col=1)
}
for (i in 1:8){
  xbar=xos+10*(UMR[i,3]-xfs)
  stdev=10*UMR[i,4]
  base=yos+22-UMR[i,1]
  hgt=1
  DrawNorm(xbar,stdev,base,hgt)		#mean,sd,base,hgt
  text(xbar+2.8*stdev,base+.2,taxon[i],pos=4,cex=1,col=1)
}
rect(4,11.7,7,12.3,col='white',border=NA)
text(5.5,12,'Merychyus',font=3,cex=1.2,col=1)
text(14,12,'Merycochoerus-Brachycrus',font=3,cex=1.2,col=1)
#================================================= Calculate rates
LB<-matrix(nrow=10,ncol=1+3*23)
  rownames(LB)=c('siouense','buwaldi','wilsoni','proprius','matthewi',
    'relictus','elegans','arenarum','minimus','crabilli')
  colnames(LB)=c('age','C1N','C1M','C1S','C2N','C2M','C2S','C3N','C3M','C3S',
    'C4N','C4M','C4S','C5N','C5M','C5S','C6N','C6M','C6S',
    'C7N','C7M','C7S','C8N','C8M','C8S','C9N','C9M','C9S',
    'C10N','C10M','C10S','C11N','C11M','C11S','C12N','C12M','C12S',
    'C13N','C13M','C13S','C14N','C14M','C14S','C15N','C15M','C15S',
    'C16N','C16M','C16S','C17N','C17M','C17S','C18N','C18M','C18S',
    'C19N','C19M','C19S','C20N','C20M','C20S','C21N','C21M','C21S',
    'C22N','C22M','C22S','C23N','C23M','C23S')
B[1:3,1:10]
for (r in 1:10){
  LB[r,1]=B[3*r-2,2]
  for (c in 1:23){
    LB[r,1+(3*c-2)]=B[3*r-2,7+2*c-1]	#sample size
    LB[r,1+(3*c-1)]=log(B[3*r-1,7+2*c-1])	#ln mean
    LB[r,1+(3*c)]=.01*B[3*r,7+2*c-1]	#stdev

  }
}
LB[,1:6]
LB=LB[c(1,3:7,9:10),]				#delete rows: buwaldi, arenarum
LB[,1:6]
LB=LB[nrow(LB):1,]				#reverse the row order
LB[,1:6]
#------------------------------------------------------------------
idr=matrix(nrow=0,ncol=9)
  colnames(idr)=c('int.g','diff.sd','rate.sd.g','log.i','log|d|','log|r|',
    'sbn','wgt','li+lr')
#---------------------------------------------- Merychyus rates
gentime=5
nr=nrow(LB[1:4,])						#number of rows
nn=.5*(nr-1)*nr;nn
nc=0								#n count
for (col in 1:23){					#column number
  for (i in 1:(nr-1)){					#increment
    for (sr in 1:(nr-i)){				#start row
      nc=nc+1
      id<-numeric(length=9)
      iyr=1000000*(LB[sr,1]-LB[sr+i,1])		#interval in years
      id[1]=iyr/gentime					#interval in generations
      meandiff=LB[sr+i,1+3*col-1]-LB[sr,1+3*col-1]	#mean diff.
      psd=PoolSD(LB[sr,1+3*col-2],LB[sr+i,1+3*col-2],	#n1,n2,sd1,sd2
        LB[sr,1+3*col],LB[sr+i,1+3*col])
      id[2]=meandiff/psd				#diff.sd
      id[3]=abs(id[2])/id[1]				#rate.sd.gen
      id[4]=log10(id[1])				#log.i
      id[5]=log10(abs(id[2]))				#log.d
      id[6]=log10(id[3])				#log.r
      if(i==1){id[7]=2}else{id[7]=3}		#sbn
      id[8]=1/id[1]					#wgt
      id[9]=id[4]+id[6]					#sum
      print(id)
      idr=rbind(idr,id)
    }
  }
}
#---------------------------------------------- Merycochoerus rates
gentime=8
nr=nrow(LB[5:8,])						#number of rows
nn=.5*(nr-1)*nr;nn
nc=0								#n count
for (col in 1:23){					#column number
  for (i in 1:(nr-1)){					#increment
    for (sr in 5:(8-i)){				#start row to run 5 to 8
      nc=nc+1
      id<-numeric(length=9)
      iyr=1000000*(LB[sr,1]-LB[sr+i,1])		#interval in years
      id[1]=iyr/gentime					#interval in generations
      meandiff=LB[sr+i,1+3*col-1]-LB[sr,1+3*col-1]	#mean diff.
      psd=PoolSD(LB[sr,1+3*col-2],LB[sr+i,1+3*col-2],	#n1,n2,sd1,sd2
        LB[sr,1+3*col],LB[sr+i,1+3*col])
      id[2]=meandiff/psd				#diff.sd
      id[3]=abs(id[2])/id[1]				#rate.sd.gen
      id[4]=log10(id[1])				#log.i
      id[5]=log10(abs(id[2]))				#log.d
      id[6]=log10(id[3])				#log.r
      if(i==1){id[7]=2}else{id[7]=3}		#sbn
      id[8]=1/id[1]					#wgt
      id[9]=id[4]+id[6]					#sum
      print(id)
      idr=rbind(idr,id)
    }
  }
}
#---------------------------------------------------------------------
nrow(idr)
min(idr[,9]);max(idr[,9])
#--------------------------------------------------------------------------
idrx=idr[!rowSums(!is.finite(idr)),]  #remove rows that have -Inf
nrow(idr);nrow(idrx)

#======================================================== Plot LRI panel (b)
#assign 'bootn' as boot number
bootn<-1000;text(-1.7,-2.6,paste("Boot n =",bootn),pos=4,cex=1,col=4)
#assign 'mode' as "medians","all","mixed"
mode<-"all";text(2,-2.6,paste("Mode: ",mode),pos=4,cex=1,col=4)
#assign circle size for points (1.5 or 2)
psize<-1.2	#1.5/2
#assign 'equation' position as "normal","lower","none" at end of each call
#send (1)idrx matrix, (2)mode(diff/rate), (3)panel placement coordinate x,
#  (4)panel placement coordinate y, (5)bootn, (6)mode, (7)psize, (8)equation
#  (9) vadd or 'vertical addition' to position of points on plot
bootresultd=PalPanelBC(idrx[,1:8],"d",1,4.5,bootn,mode,psize,"normal",0)
proc.time()-ptm
#======================================================== Plot LRI panel (c)
#send idrx matrix, n, mode(diff/rate), panel placement coordinates, mode
bootresultr=PalPanelBC(idrx[,1:8],"r",13,4.5,bootn,mode,psize,"normal",0)
proc.time()-ptm

#----------------------------------------------------------
text(7,-2.6,paste("Processing time: ",round(proc.time()[3]-ptm[3],digits=2),
  "seconds"),pos=4,cex=1,col=4)




