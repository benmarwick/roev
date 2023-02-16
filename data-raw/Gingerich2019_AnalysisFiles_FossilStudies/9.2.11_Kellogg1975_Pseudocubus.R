##Philip D. Gingerich RATES OF EVOLUTION (Cambridge University Press 2019) 
##9.2.11_Kellogg1975_Pseudocubus
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
#============================================================ Load file
file1<-"C://R_aaROEVchapt09//9.2.11_Kellogg1975_Pseudocubus.csv"
K<-read.csv(file1,header=TRUE,sep=",",quote="\"",dec=".",fill=TRUE)
nrow(K)
#------------------------------------------------------Logged matrix
K[1:3,]
LK<-matrix(nrow=nrow(K),ncol=4);colnames(LK)=c('age','N','LnM','LnS')
for (i in 1:nrow(K)){
  LK[i,1]=K[i,2]
  LK[i,2]=K[i,3]
  LK[i,3]=log(K[i,4])
  LK[i,4]=K[i,6]/K[i,4]
}
LK[1:3,]
min(LK[,3]-3*LK[,4]);max(LK[,3]+3*LK[,4])
min(LK[,1]);max(LK[,1])

#======================================================= Plot panel a
xos=3;xoe=18;yos=9.5;yoe=16.5 #17.5			#x,y original start and end
xfs=4;xfe=5.4;yfs=-5;yfe=-2			#x,y foreground start and end
xfm=(xoe-xos)/(xfe-xfs)				#x foreground multiplier
yfm=(yoe-yos)/(yfe-yfs)				#y foreground multiplier
lines(c(xos,xoe),c(yos,yos),lwd=1,lty=1,col=1)
lines(c(xos,xos),c(yos,yoe),lwd=1,lty=1,col=1)
#---------------------------------------------------------axis labels
for (i in seq(yos,yoe-.5,1)){				#y-axis
  lines(c(xos-.15,xos),c(i,i),lwd=1,lty=1,col=1)
}
for (i in seq(yos,yoe-.5,1)){
  text(xos,i,format(round(abs(.5*i-9.75),digits=1),nsmall=1),pos=2,cex=1,col=1)
}
for (i in xos+10*(seq(xfs,xfe,.1)-xfs)){		#x-axis
  lines(c(i,i),c(yos,yos-.15),lwd=1,lty=1,col=1)
}
for (i in xos+10*(seq(xfs,xfe,.1)-xfs)){
  #text(i,yos,3.1+.1*i,pos=1,cex=1,col=1)
  text(i,yos,format(round(3.7+.1*i,digits=1),nsmall=1),
    pos=1,cex=1,col=1)
}
text(xos,yos+7.0,expression(paste(italic('Pseudocubus vema')*
	': width of thorax')),pos=4,cex=1.3,col=1)
text(xos-1.5,yos+3.5,'Geological age (Ma)',
	srt=90,cex=1.3,col=1)
text(xos+7.5,yos-1.2,'Ln measurement (microns)',
	cex=1.3,col=1)
#---------------------------------------------------------------- Plot ranges
LK[1:3,];LK[,1];LK[,3]
for (i in 1:34){
  lines(xos+10*c(LK[i,3]-LK[i,4],LK[i,3]+LK[i,4])-40,
    c(yos+10+2*LK[i,1],yos+10+2*LK[i,1]),lty=1,lwd=1,col=gray(6/10))
}
lines(xos+10*LK[,3]-40,yos+10+2*LK[,1],lty=1,lwd=2,col=1)

#================================================================ Rates 
idr=matrix(nrow=0,ncol=9)
  colnames(idr)=c('int.g','diff.sd','rate.sd.g','log.i','log|d|','log|r|',
    'sbn','wgt','li+lr')
#---------------------------------------------------------------- rates
LK[1:3,]
gt=1	#generation time 1 year fide Kellogg (1975)
nr=nrow(LK)					#number of rows
nn=.5*(nr-1)*nr;nn
nc=0
  for (i in 1:(nr-1)){			#increment
    for (sr in 1:(nr-i)){		#start row
      id<-numeric(length=9)
      nc=nc+1
      iyr=1000000*(LK[sr+i,1]-LK[sr,1])	#interval in years
      #w1g=1000*LZ[sr,5];w2g=1000*LZ[(sr+i),5];wg=exp(.5*(log(w1g)+log(w2g)))
      #10^(.266*log10(wg)-.553)		#Eq 3.2 generation time in years
      id[1]=iyr/gt					#interval in generations
      meandiff=LK[(sr+i),3]-LK[sr,3]		#mean diff.
      psd=PoolSD(LK[sr,2],LK[sr+i,2],LK[sr,4],LK[sr+i,4])#n1,n2,sd1,sd2
      id[2]=meandiff/psd				#diff.sd
      id[3]=abs(id[2])/id[1]				#rate.sd.gen
      id[4]=log10(id[1])				#log.i
      id[5]=log10(abs(id[2]))				#log.d
      id[6]=log10(id[3])				#log.r
      if(i==1){id[7]=2}else{id[7]=3}		#sbn
      id[8]=1/id[1]					#wgt
      id[9]=id[4]+id[6]					#sum
    	#print(id)
    	idr=rbind(idr,id)
    }
  }
nn;nrow(idr)
#--------------------------------------------------------------------------
idrx=idr[!rowSums(!is.finite(idr)),]  #remove rows that have -Inf
nn;nrow(idr);nrow(idrx)
min(idrx[,9]);max(idrx[,9])
writefile<-"C://R_aaROEVchapt09//9.2.11_Kellogg1975_Pseudocubus_out.csv"
write.csv(idrx,file=writefile,na="NA",row.names=FALSE)
 
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
bootresultr=PalPanelBC(idrx[,1:8],"r",13,4.5,bootn,mode,psize,"normal",2)
proc.time()-ptm
#----------------------------------------------------------
text(7,-2.6,paste("Processing time: ",round(proc.time()[3]-ptm[3],digits=2),
  "seconds"),pos=4,cex=1,col=4)



