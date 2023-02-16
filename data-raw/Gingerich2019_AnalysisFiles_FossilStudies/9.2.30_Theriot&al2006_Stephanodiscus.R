##Philip D. Gingerich RATES OF EVOLUTION (Cambridge University Press 2019) 
##9.2.30_Theriot&al2006_Stephanodiscus
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
	#mgp=c(2.2,.3,0),			#margin for axis title/labels/tickline
	mgp=c(2,.3,0),			#margin for axis title/labels/tickline
	tck=-0.01,				#tick-mark length
	xaxp=c(xr[1],xr[2],xr[1]+xr[2]),		#x-axis ticks, no. ticks
	yaxp=c(yr[1],yr[2],abs(yr[1]-yr[2])),	#y-axis ticks, no. ticks
	cex.axis=.8,
	cex.lab=1.1)			#label size
dev.size('in')
#============================================================ Load file
file1<-"C://R_aaROEVchapt09//9.2.30_Theriot&al2006_Stephanodiscus.csv" 
T<-read.csv(file1,header=TRUE,sep=",",quote="\"",dec=".",fill=TRUE)
T[1:5,];nrow(T);ncol(T)
#------------------------------------------------------Logged matrices
TT=apply(as.matrix.noquote(T[1:nrow(T),]),2,as.numeric)
LT<-matrix(nrow=nrow(TT),ncol=3*4+1)
  colnames(LT)=c('Age','N1','X1','S1','N2','X2','S2','N3','X3','S3',
    'N4','X4','S4')
LT[,1]=TT[,1]						#ages
LT[,c(2,5,8,11)]=TT[,c(2,2,2,2)]			#Ns
LT[,c(3,6,9)]=log(TT[,c(3,4,5)]);LT[,12]=T[,9]	#means
LT[,c(4,7,10)]=TT[,c(6,7,8)];LT[,13]=T[,10]	#stdev
TT[1:5,]
LT[1:5,]
min(LT[,1]);max(LT[,1])
min(na.omit(LT[,12]))-min(na.omit(LT[,13]))
max(na.omit(LT[,12]))+max(na.omit(LT[,13]))

#======================================================= Plot panel a
xos=2;xoe=17;yos=9.5;yoe=17.5			#x,y original start and end
xfs=0.1;xfe=.8;yfs=-14;yfe=0			#x,y foreground start and end
xfm=(xoe-xos)/(xfe-xfs)				#x foreground multiplier
yfm=(yoe-yos)/(yfe-yfs)				#y foreground multiplier
lines(c(xos,xoe),c(yos,yos),lwd=1,lty=1,col=1)
lines(c(xos,xos),c(yos,yoe-1),lwd=1,lty=1,col=1)
#---------------------------------------------------------axis labels
for (i in seq(yos,yoe-.6,1)){				#y-axis
  lines(c(xos-.15,xos),c(i,i),lwd=1,lty=1,col=1)
}
for (i in seq(yos,yoe-.6,1)){
  text(xos,i,format(round(abs(2*i-19-14),digits=0),nsmall=0),pos=2,cex=1,col=1)
}
for (i in seq(xos,xoe,1)){		#x-axis
  lines(c(i,i),c(yos,yos-.15),lwd=1,lty=1,col=1)
}
for (i in seq(xos,xoe,2)){		#x-axis
  text(i+1,yos,format(round(i/20,digits=1),nsmall=1),
    pos=1,cex=1,col=1)
}
text(xos,yos+7.5,expression(paste(italic('Stephanodiscus')*
  ': relative spine count')),pos=4,cex=1.3,col=1)
text(xos-1.2,yos+3.5,'Age (Ka)',
	srt=90,cex=1.3,col=1)
text(xos+7.5,yos-1.2,
  expression(paste('Slope of spine count on diameter')),cex=1.3,col=1)
#---------------------------------------------------------------- Plot ranges
lines(c(3,18),c(11.5,11.5),lty=2,lwd=1,col=gray(6/10))
LT[1:5,]
min(na.omit(LT[,3]));max(na.omit(LT[,3]))
min(LT[,1]);max(LT[,1])
xpos=-1;ypos=7
for (i in 1:nrow(LT)){
  lines(xos+20*c(LT[i,12]-LT[i,13],LT[i,12]+LT[i,13])+xpos,
    yos+.0005*c(LT[i,1],LT[i,1])+ypos,lty=1,lwd=1,col=gray(6/10))
}
for (i in 1:nrow(LT)){
  points(xos+20*LT[i,12]+xpos,yos+.0005*LT[i,1]+ypos,pch=20,cex=.8,col=1)
}
lines(xos+20*LT[,12]+xpos,yos+.0005*LT[,1]+ypos,lty=1,lwd=1.5,col=1)

text(xos+7,14,expression(paste(italic('S. yellowstonensis'))),
  pos=4,cex=1,col=1)
text(xos+10,11,expression(paste(italic('S. niagarae to S. yellowstonensis'))),
  pos=4,cex=1,col=1)
text(xos+13,10.5,expression(paste('transition')),pos=4,cex=1,col=1)

#====================================================== Lineages
LT1=LT[45:64,];nrow(LT1)
LT2=LT[1:44,];nrow(LT2)
#====================================================== Rates1: transition
idr=matrix(nrow=0,ncol=9)
  colnames(idr)=c('int.g','diff.sd','rate.sd.g','log.i','log|d|','log|r|',
    'sbn','wgt','li+lr')
#------------------------------------------------ Rates
LT[1:3,]
gt=4	#generation time assumed
nr1=nrow(LT1);nr1					#number of rows
nn1=.5*(nr1-1)*nr1;nn1
zc1=0							#zero count
nc=0
for (cn in c(3,6,9,12)){				#columns of means
  for (i in 1:(nr1-1)){				#increment
    for (sr in 1:(nr1-i)){			#start row
      id<-numeric(length=9)
      nc=nc+1
      iyr=LT1[sr,1]-LT1[sr+i,1]				#interval in years
      id[1]=iyr/gt					#interval in generations
      meandiff=LT1[sr,cn]-LT1[sr+i,cn]		#mean diff.
      if (!is.na(meandiff)){if (meandiff==0){zc1=zc1+1}}
      psd=PoolSD(LT1[sr,cn-1],LT1[sr+i,cn-1],LT1[sr,cn+1],LT1[sr+i,cn+1])#n1,n2,sd1,sd2
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
#--------------------------------------------------------------------------
nrow(idr)
#--------------------------------------------------------------------------
idrx1=idr[!rowSums(!is.finite(idr)),]  #remove rows that have -Inf
nn1;nrow(idr);nrow(idrx1);zc1
min(idrx1[,9]);max(idrx1[,9])

#=================================================== Rates2: yellowstonensis
idr=matrix(nrow=0,ncol=9)
  colnames(idr)=c('int.g','diff.sd','rate.sd.g','log.i','log|d|','log|r|',
    'sbn','wgt','li+lr')
#------------------------------------------------ Rates
LT2[1:3,]
gt=4	#generation time assumed
nr2=nrow(LT2);nr2					#number of rows
nn2=.5*(nr2-1)*nr2;nn2
zc2=0							#zero count
nc=0
for (cn in c(3,6,9,12)){				#columns of means
  for (i in 1:(nr2-1)){				#increment
    for (sr in 1:(nr2-i)){			#start row
      id<-numeric(length=9)
      nc=nc+1
      iyr=LT2[sr,1]-LT2[sr+i,1]				#interval in years
      id[1]=iyr/gt					#interval in generations
      meandiff=LT2[sr,cn]-LT2[sr+i,cn]		#mean diff.
      if (!is.na(meandiff)){if (meandiff==0){zc2=zc2+1}}
      psd=PoolSD(LT2[sr,cn-1],LT2[sr+i,cn-1],LT2[sr,cn+1],LT2[sr+i,cn+1])#n1,n2,sd1,sd2
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
#--------------------------------------------------------------------------
nrow(idr)
#--------------------------------------------------------------------------
idrx2=idr[!rowSums(!is.finite(idr)),]  #remove rows that have -Inf
nn2;nrow(idr);nrow(idrx2);zc2
min(idrx2[,9]);max(idrx2[,9])

#--------------------------------------------------------------------Summary
idrxa=rbind(idrx1,idrx2)
zca=zc1+zc2
zca;nrow(idrxa);zca+nrow(idrxa);nn1+nn2
writefile<-"C://R_aaROEVchapt09//9.2.30_Theriot&al2006_Stephanodiscus_out.csv"
write.csv(idrxa,file=writefile,na="NA",row.names=FALSE) 
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
bootresultd=PalPanelBC(idrx1[,1:8],"r",1,4.5,bootn,mode,psize,"normal",0)
proc.time()-ptm
#======================================================== Plot LRI panel (c)
bootresultr=PalPanelBC(idrx2[,1:8],"r",13,4.5,bootn,mode,psize,"normal",0)
proc.time()-ptm

text(9.2,4,expression(paste("Transition interval")),pos=2,cex=1,col=1)
text(21.2,4,expression(paste(italic("S. yellowstonensis"))),pos=2,cex=1,col=1)

#----------------------------------------------------------
text(7,-2.6,paste("Processing time: ",round(proc.time()[3]-ptm[3],digits=2),
  "seconds"),pos=4,cex=1,col=4)


