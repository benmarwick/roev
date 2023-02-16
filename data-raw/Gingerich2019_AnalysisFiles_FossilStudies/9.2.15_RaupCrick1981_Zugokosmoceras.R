##Philip D. Gingerich RATES OF EVOLUTION (Cambridge University Press 2019)
##9.2.15_RaupCrick1981_Zugokosmoceras
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
xos=2;xoe=19;yos=9.5;yoe=17.5			#x,y original start and end
xfs=3.4;xfe=5.1;yfs=-22;yfe=-15		#x,y foreground start and end
xfm=(xoe-xos)/(xfe-xfs)				#x foreground multiplier
yfm=(yoe-yos)/(yfe-yfs)				#y foreground multiplier
lines(c(xos,xoe),c(yos,yos),lwd=1,lty=1,col=1)
lines(c(xos,xos),c(yos,yoe),lwd=1,lty=1,col=1)

#---------------------------------------------------------axis labels
for (i in seq(yos+.5,yoe,1)){				#y-axis
  lines(c(xos-.15,xos),c(i,i),lwd=1,lty=1,col=1)
}
for (i in seq(yos+.5,yoe,1)){
  text(xos+.1,i,format(round(165.2-.2*(i-10),digits=1),nsmall=1),pos=2,cex=1,col=1)
}
for (i in xos+10*(seq(xfs,xfe,.1)-xfs)){		#x-axis
  lines(c(i,i),c(yos,yos-.15),lwd=1,lty=1,col=1)
}
for (i in xos+10*(seq(xfs,xfe,.1)-xfs)){
  text(i,yos,format(round(3.2+.1*i,digits=1),nsmall=1),
    pos=1,cex=1,col=1)
}
text(xos,yos+8.0,expression(paste(italic('Kosmoceras')*
	': shell diameter')),pos=4,cex=1.3,col=1)
text(xos-2,yos+4,'Geological age (Ma)',
	srt=90,cex=1.3,col=1)
text(xos+9,yos-1.2,'Ln measurement (mm)',
	cex=1.3,col=1)
#============================================================ Load file
file1<-"../Gingerich2019_AnalysisFiles_FossilStudies/9.2.15_Brinkmann1929_Zugokosmoceras.csv"
file2<-"../Gingerich2019_AnalysisFiles_FossilStudies/9.2.15_Brinkmann1929_Anakosmoceras.csv"
B1<-read.csv(file1,header=TRUE,sep=",",quote="\"",dec=".",fill=TRUE)
B2<-read.csv(file2,header=TRUE,sep=",",quote="\"",dec=".",fill=TRUE)
#------------------------------------------------- Ln matrix for B1
B1[1:5,]
LB1<-matrix(nrow=nrow(B1),ncol=4);colnames(LB1)=c('Ma','N','LnM','LnS')
for (i in 1:nrow(B1)){
  LB1[i,1]=B1[i,8]*1.4/1400		#Ma as Ma/cm
  LB1[i,2]=B1[i,3]				#N
  LB1[i,3]=log(B1[i,4])			#LnM
  LB1[i,4]=B1[i,6]/100			#LnS
}
LB1[1:5,];LB1[(nrow(B1)-2):nrow(B1),];nrow(LB1)
#------------------------------------------------- Ln matrix for B2
B2[1:5,]
LB2<-matrix(nrow=nrow(B2),ncol=4);colnames(LB2)=c('Ma','N','LnM','LnS')
for (i in 1:nrow(B2)){
  LB2[i,1]=B2[i,8]*1.4/1400		#Ma as Ma/cm
  LB2[i,2]=B2[i,3]				#N
  LB2[i,3]=log(B2[i,4])			#LnM
  LB2[i,4]=B2[i,6]/100			#LnS
}
LB2[1:5,];LB2[(nrow(B2)-2):nrow(B2),];nrow(LB2)
#-------------------------------------------------------- Plot ranges for B1
xpos=34
for (i in 1:nrow(LB1)){
  lines(xos+10*c(LB1[i,3]-LB1[i,4],LB1[i,3]+LB1[i,4])-xpos,
  yos+.5+5*c(LB1[i,1],LB1[i,1]),
  lty=1,lwd=1,col=gray(6/10))
}
for (i in 1:nrow(LB1)){
  points(xos+10*LB1[i,3]-xpos,yos+.5+5*LB1[i,1],pch=20,cex=.8,col=1)
}
lines(xos+10*LB1[,3]-xpos,yos+.5+5*LB1[,1],lty=1,lwd=1.5,col=1)
#-------------------------------------------------------- Plot ranges for B2
for (i in 1:nrow(LB2)){
  lines(xos+10*c(LB2[i,3]-LB2[i,4],LB2[i,3]+LB2[i,4])-xpos,
  yos+.5+5*c(LB2[i,1],LB2[i,1]),
  lty=1,lwd=1,col=gray(6/10))
}
for (i in 1:nrow(LB2)){
  points(xos+10*LB2[i,3]-xpos,yos+.5+5*LB2[i,1],pch=20,cex=.8,col=1)
}
lines(xos+10*LB2[,3]-xpos,yos+.5+5*LB2[,1],lty=1,lwd=1.5,col=1)
#-------------------------------------------------------- Label species
text(5,16.6,expression(paste(italic("'Anakosmoceras'"))),pos=1,cex=1,col=1)
text(5,16.0,expression(paste("(microconch)")),pos=1,cex=1,col=1)
text(13,16.6,expression(paste(italic("'Zugokosmoceras'"))),pos=1,cex=1,col=1)
text(13,16.0,expression(paste("(macroconch)")),pos=1,cex=1,col=1)

#================================================== Calculate rates for B1
idr=matrix(nrow=0,ncol=9)
  colnames(idr)=c('int.g','diff.sd','rate.sd.g','log.i','log|d|','log|r|',
    'sbn','wgt','li+lr')
#------------------------------------------------------------------------
LB1[1:3,]
gt=15 							#generation time in years
nr=nrow(LB1)							#number of rows
nn=.5*(nr-1)*nr;nn					#expected number of rates
nc=0								#n count
for (col in 3:3){						#column number
  for (i in 1:(nr-1)){					#increment
    for (sr in 1:(nr-i)){				#start row
      nc=nc+1
      id<-numeric(length=9)
      iyr=1000000*(LB1[sr+i,1]-LB1[sr,1])		#interval in years
      id[1]=iyr/gt					#interval in generations
      meandiff=LB1[sr+i,col]-LB1[sr,col]		#mean diff.
      psd=PoolSD(LB1[sr,col-1],LB1[sr+i,col-1],	#n1,n2,sd1,sd2
        LB1[sr,col+1],LB1[sr+i,col+1])
      id[2]=meandiff/psd				#diff.sd
      id[3]=abs(id[2])/id[1]				#rate.sd.gen
      id[4]=log10(id[1])				#log.i
      id[5]=log10(abs(id[2]))				#log.d
      id[6]=log10(id[3])				#log.r
      if(i==1){id[7]=2}else{id[7]=3}		#sbn
      id[8]=1/id[1]					#wgt
      id[9]=id[4]+id[6]					#sum
      idr=rbind(idr,id)
    }
  }
}
#---------------------------------------------------------------------
nrow(idr)
#--------------------------------------------------------------------------
idrx1=idr[!rowSums(!is.finite(idr)),]  #remove rows that have -Inf
nrow(idr);nrow(idrx1)
min(idrx1[,9]);max(idrx1[,9])
#================================================== Calculate rates for B2
idr=matrix(nrow=0,ncol=9)
  colnames(idr)=c('int.g','diff.sd','rate.sd.g','log.i','log|d|','log|r|',
    'sbn','wgt','li+lr')
#------------------------------------------------------------------------
LB2[1:3,]
gt=15 							#generation time in years
nr=nrow(LB2)							#number of rows
nn=.5*(nr-1)*nr;nn					#expected number of rates
nc=0								#n count
for (col in 3:3){						#column number
  for (i in 1:(nr-1)){					#increment
    for (sr in 1:(nr-i)){				#start row
      nc=nc+1
      id<-numeric(length=9)
      iyr=1000000*(LB2[sr+i,1]-LB2[sr,1])		#interval in years
      id[1]=iyr/gt					#interval in generations
      meandiff=LB2[sr+i,col]-LB2[sr,col]		#mean diff.
      psd=PoolSD(LB2[sr,col-1],LB2[sr+i,col-1],	#n1,n2,sd1,sd2
        LB2[sr,col+1],LB2[sr+i,col+1])
      id[2]=meandiff/psd				#diff.sd
      id[3]=abs(id[2])/id[1]				#rate.sd.gen
      id[4]=log10(id[1])				#log.i
      id[5]=log10(abs(id[2]))				#log.d
      id[6]=log10(id[3])				#log.r
      if(i==1){id[7]=2}else{id[7]=3}		#sbn
      id[8]=1/id[1]					#wgt
      id[9]=id[4]+id[6]					#sum
      idr=rbind(idr,id)
    }
  }
}
#---------------------------------------------------------------------
nrow(idr)
#--------------------------------------------------------------------------
idrx2=idr[!rowSums(!is.finite(idr)),]  #remove rows that have -Inf
nrow(idr);nrow(idrx2)
min(idrx2[,9]);max(idrx2[,9])

#=======================================================================
idrx=rbind(idrx2,idrx1)

writefile<-"C://R_aaROEVchapt09//9.2.15_RaupCrick1981_Zugokosmoceras_out.csv"
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
#  (9)vert. scale/point adjustment
bootresultd=PalPanelBC(idrx2[,1:8],"r",1,4.5,bootn,mode,psize,"normal",1)
proc.time()-ptm
#======================================================== Plot LRI panel (c)
#send idrx matrix, n, mode(diff/rate), panel placement coordinates, mode
bootresultr=PalPanelBC(idrx1[,1:8],"r",13,4.5,bootn,mode,psize,"normal",1)
proc.time()-ptm

text(7,5.6,expression(paste(italic("'Anakosmoceras'"))),pos=1,cex=1,col=1)
text(19,5.6,expression(paste(italic("'Zugokosmoceras'"))),pos=1,cex=1,col=1)

#----------------------------------------------------------
text(7,-2.6,paste("Processing time: ",round(proc.time()[3]-ptm[3],digits=2),
  "seconds"),pos=4,cex=1,col=4)




