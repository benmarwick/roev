##Philip D. Gingerich RATES OF EVOLUTION (Cambridge University Press 2019)
##9.2.12_Gingerich1976_Hyopsodus
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
file1<-"../Gingerich2019_AnalysisFiles_FossilStudies//9.2.12_Gingerich1976_Hyopsodus.csv"
G<-read.csv(file1,header=TRUE,sep=",",quote="\"",dec=".",fill=TRUE)
nrow(G);G[1:20,1:10]
#------------------------------------------------------Logged matrix
LG<-matrix(nrow=nrow(G),ncol=4)
  colnames(LG)=c('loc','spec','LnLxW','age')
LG[,1]=G[,2]
LG[,2]=G[,3]
LG[,3]=G[,6]
LG[,4]=G[,10];LG[1:5,]
LGS<-matrix(nrow=0,ncol=5);colnames(LGS)=c('Samp','Age','N','LnM','LnS')
lgsc=0	#lgs count
ss<-matrix(nrow=0,ncol=4);colnames(ss)=c('loc','spec','LnLxW','age')
ss2=LG[1,2];ss4=LG[1,4];ss2;ss4	#species(2) and age(4)
for (i in 1:nrow(LG)){
  if (LG[i,2]==ss2&&LG[i,4]==ss4){
    ss=rbind(ss,LG[i,])
  }else{
    LGS=rbind(LGS,c(NA,NA,NA,NA,NA))
    lgsc=lgsc+1
    LGS[lgsc,1]=lgsc
    LGS[lgsc,2]=ss[nrow(ss),4]
    LGS[lgsc,3]=nrow(ss)
    LGS[lgsc,4]=mean(ss[,3])
    LGS[lgsc,5]=sd(ss[,3])
    ss<-matrix(nrow=0,ncol=4);colnames(ss)=c('loc','spec','LnLxW','age')
    ss=rbind(ss,LG[i,])
    ss2=LG[i,2];ss4=LG[i,4]
  }
}
LGS=rbind(LGS,c(NA,NA,NA,NA,NA))
lgsc=lgsc+1
LGS[lgsc,1]=lgsc
LGS[lgsc,2]=ss[nrow(ss),4]
LGS[lgsc,3]=nrow(ss)
LGS[lgsc,4]=mean(ss[,3])
LGS[lgsc,5]=sd(ss[,3])
LGS
#====================================================================
scat=FALSE #TRUE
if (scat==TRUE){
plot(LGS[,4],LGS[,2])
for (i in 1:nrow(LGS)){
  lines(c(LGS[i,4]-LGS[i,5],LGS[i,4]+LGS[i,5]),c(LGS[i,2],LGS[i,2]),
    lty=1,lwd=2,col=2)
}
for (i in 2:22){
  lines(c(LGS[i-1,4],LGS[i,4]),c(LGS[i-1,2],LGS[i,2]),
    lty=1,lwd=2,col=1)
}
for (i in 26:38){
  lines(c(LGS[i-1,4],LGS[i,4]),c(LGS[i-1,2],LGS[i,2]),
    lty=1,lwd=2,col=1)
}
for (i in 40:42){
  lines(c(LGS[i-1,4],LGS[i,4]),c(LGS[i-1,2],LGS[i,2]),
    lty=1,lwd=2,col=1)
}
lines(LGS[c(7,23),4],LGS[c(7,23),2],lty=2,lwd=1,col=1)
lines(LGS[c(9,24,25),4],LGS[c(9,24,25),2],lty=2,lwd=1,col=1)
lines(LGS[c(29,39),4],LGS[c(29,39),2],lty=2,lwd=1,col=1)
}
#======================================================= Plot panel a
xos=3;xoe=18;yos=9.5;yoe=17.5			#x,y original start and end
xfs=1.8;xfe=3.3;yfs=-55;yfe=-53			#x,y foreground start and end
xfm=(xoe-xos)/(xfe-xfs)				#x foreground multiplier
yfm=(yoe-yos)/(yfe-yfs)				#y foreground multiplier
lines(c(xos,xoe),c(yos,yos),lwd=1,lty=1,col=1)
lines(c(xos,xos),c(yos,yoe),lwd=1,lty=1,col=1)
#---------------------------------------------------------axis labels
for (i in seq(yos,yoe-.5,1)){				#y-axis
  lines(c(xos-.15,xos),c(i,i),lwd=1,lty=1,col=1)
}
for (i in seq(yos,yoe-.5,1)){
  text(xos+.1,i,format(round(abs(56.8-.2*(i-.5)),digits=1),nsmall=1),pos=2,cex=1,col=1)
}

for (i in xos+10*(seq(xfs,xfe,.1)-xfs)){		#x-axis
  lines(c(i,i),c(yos,yos-.15),lwd=1,lty=1,col=1)
}
for (i in xos+10*(seq(xfs,xfe,.1)-xfs)){
  text(i,yos,format(round(1.5+.1*i,digits=1),nsmall=1),
    pos=1,cex=1,col=1)
}
text(xos,yos+8,expression(paste(italic('Hyopsodus')~
	'spp.:  tooth size')),pos=4,cex=1.3,col=1)
text(xos-1.8,yos+4,'Geological age (Ma)',
	srt=90,cex=1.3,col=1)
text(xos+7.5,yos-1.2,
  expression(paste('Ln'~M[1]~'crown area ('*mm^2*')')),cex=1.3,col=1)
#------------------------------------------------------- plus/minus sd lines
xpos=-15;ypos=284
for (i in 1:nrow(LGS)){
  lines(xpos+10*c(LGS[i,4]-LGS[i,5],LGS[i,4]+LGS[i,5]),
    .5+ypos+5*c(LGS[i,2],LGS[i,2]),lty=1,lwd=1,col=gray(6/10))
}
for (i in 1:nrow(LGS)){
  points(xpos+10*LGS[i,4],.5+ypos+5*LGS[i,2],pch=20,cex=.8,col=1)
}

#------------------------------------------------------- mean lines
for (i in 2:22){
  lines(c(xpos+10*LGS[i-1,4],xpos+10*LGS[i,4]),
    .5+c(ypos+5*LGS[i-1,2],ypos+5*LGS[i,2]),lty=1,lwd=1.5,col=1)
}
for (i in c(26:38)){
  lines(c(xpos+10*LGS[i-1,4],xpos+10*LGS[i,4]),
    .5+c(ypos+5*LGS[i-1,2],ypos+5*LGS[i,2]),lty=1,lwd=1.5,col=1)
}
for (i in c(40:42)){
  lines(c(xpos+10*LGS[i-1,4],xpos+10*LGS[i,4]),
    .5+c(ypos+5*LGS[i-1,2],ypos+5*LGS[i,2]),lty=1,lwd=1.5,col=1)
}
lines(xpos+10*LGS[c(3,23),4],
  .5+ypos+5*LGS[c(3,23),2],lty=2,lwd=1,col=1)
lines(xpos+10*LGS[c(9,24,25),4],
  .5+ypos+5*LGS[c(9,24,25),2],lty=2,lwd=1,col=1)
lines(xpos+10*LGS[c(29,39),4],
  .5+ypos+5*LGS[c(29,39),2],lty=2,lwd=1,col=1)

text(7.9,9.8,expression(paste(italic('H. loomisi'))),pos=4,cex=1,col=1)
text(9.3,12.15,expression(paste(italic('H. latidens'))),pos=4,cex=1,col=1)
rect(3.1,12.8,6,13.3,col='white',border=NA)
text(3,13,expression(paste(italic('H. pauxillus'))),pos=4,cex=1,col=1)
rect(5.8,14.4,13,14.9,col='white',border=NA)
text(5.8,14.65,expression(paste(italic('H. minor'))),pos=4,cex=1,col=1)
text(10,14.65,expression(paste(italic('H. miticulus'))),pos=4,cex=1,col=1)
text(4.6,16.8,expression(paste(italic('H. lysitensis'))),pos=4,cex=1,col=1)
text(14.2,15.7,expression(paste(italic('H. powellianus'))),pos=4,cex=1,col=1)

#========================================================= Species lineages
SL1=LGS[1:22,];SL1
SL2=LGS[25:38,];SL2
SL3=LGS[39:42,];SL3
#=============================================================== Rates SL1
idr=matrix(nrow=0,ncol=9)
  colnames(idr)=c('int.g','diff.sd','rate.sd.g','log.i','log|d|','log|r|',
    'sbn','wgt','li+lr')
#------------------------------------------------------
SL1[1:3,]
gt=2	#generation time assumed
nr=nrow(SL1)						#number of rows
nn1=.5*(nr-1)*nr;nn1
zc1=0							#zero count
nc=0
cn=4
  for (i in 1:(nr-1)){				#increment
    for (sr in 1:(nr-i)){			#start row
      id<-numeric(length=9)
      nc=nc+1
      iyr=1000000*(SL1[sr+i,2]-SL1[sr,2])	#interval in years
      #w1g=1000*LZ[sr,5];w2g=1000*LZ[(sr+i),5];wg=exp(.5*(log(w1g)+log(w2g)))
      #10^(.266*log10(wg)-.553)		#Eq 3.2 generation time in years
      id[1]=iyr/gt					#interval in generations
      meandiff=SL1[sr+i,cn]-SL1[sr,cn]		#mean diff.
      if (!is.na(meandiff)){if (meandiff==0){zc1=zc1+1}}
      psd=PoolSD(SL1[sr+i,cn-1],SL1[sr,cn-1],SL1[sr+i,cn+1],SL1[sr,cn+1])#n1,n2,sd1,sd2
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
#--------------------------------------------------------------------------
nrow(idr)
#--------------------------------------------------------------------------
idrx1=idr[!rowSums(!is.finite(idr)),]  #remove rows that have -Inf
nn1;nrow(idr);nrow(idrx1);zc1
min(idrx1[,9]);max(idrx1[,9])
#=============================================================== Rates SL2
idr=matrix(nrow=0,ncol=9)
  colnames(idr)=c('int.g','diff.sd','rate.sd.g','log.i','log|d|','log|r|',
    'sbn','wgt','li+lr')
#------------------------------------------------------
gt=2	#generation time assumed
nr=nrow(SL2)						#number of rows
nn2=.5*(nr-1)*nr;nn2
zc2=0							#zero count
nc=0
cn=4
  for (i in 1:(nr-1)){				#increment
    for (sr in 1:(nr-i)){			#start row
      id<-numeric(length=9)
      nc=nc+1
      iyr=1000000*(SL2[sr+i,2]-SL2[sr,2])	#interval in years
      id[1]=iyr/gt					#interval in generations
      meandiff=SL2[sr+i,cn]-SL2[sr,cn]		#mean diff.
      if (!is.na(meandiff)){if (meandiff==0){zc2=zc2+1}}
      psd=PoolSD(SL2[sr+i,cn-1],SL2[sr,cn-1],SL2[sr+i,cn+1],SL2[sr,cn+1])#n1,n2,sd1,sd2
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
#--------------------------------------------------------------------------
nrow(idr)
#--------------------------------------------------------------------------
idrx2=idr[!rowSums(!is.finite(idr)),]  #remove rows that have -Inf
nn2;nrow(idr);nrow(idrx2);zc2
min(idrx2[,9]);max(idrx2[,9])
#=============================================================== Rates SL3
idr=matrix(nrow=0,ncol=9)
  colnames(idr)=c('int.g','diff.sd','rate.sd.g','log.i','log|d|','log|r|',
    'sbn','wgt','li+lr')
#------------------------------------------------------
gt=2	#generation time assumed
nr=nrow(SL3)						#number of rows
nn3=.5*(nr-1)*nr;nn3
zc3=0							#zero count
nc=0
cn=4
  for (i in 1:(nr-1)){				#increment
    for (sr in 1:(nr-i)){			#start row
      id<-numeric(length=9)
      nc=nc+1
      iyr=1000000*(SL3[sr+i,2]-SL3[sr,2])	#interval in years
      id[1]=iyr/gt					#interval in generations
      meandiff=SL3[sr+i,cn]-SL3[sr,cn]		#mean diff.
      if (!is.na(meandiff)){if (meandiff==0){zc3=zc3+1}}
      psd=PoolSD(SL3[sr+i,cn-1],SL3[sr,cn-1],SL3[sr+i,cn+1],SL3[sr,cn+1])#n1,n2,sd1,sd2
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
#--------------------------------------------------------------------------
nrow(idr)
#--------------------------------------------------------------------------
idrx3=idr[!rowSums(!is.finite(idr)),]  #remove rows that have -Inf
nn3;nrow(idr);nrow(idrx3);zc3
min(idrx3[,9]);max(idrx3[,9])
idrx=rbind(idrx1,idrx2,idrx3)
nrow(idrx);zc1+zc2+zc3;nn1+nn2+nn3
#--------------------------------------------------------------------------
writefile<-"C://R_aaROEVchapt09//9.2.12_Gingerich1976_Hyopsodus_out.csv"
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
bootresultd=PalPanelBC(idrx[,1:8],"r",1,4.5,bootn,mode,psize,"normal",1)
proc.time()-ptm
#======================================================== Plot LRI panel (c)
#send idrx matrix, n, mode(diff/rate), panel placement coordinates, mode
bootresultr=PalPanelBC(idrx2[,1:8],"r",13,4.5,bootn,mode,psize,"normal",1)
proc.time()-ptm

#----------------------------------------------------------
text(7,-2.6,paste("Processing time: ",round(proc.time()[3]-ptm[3],digits=2),
  "seconds"),pos=4,cex=1,col=4)


