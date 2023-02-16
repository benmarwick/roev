##Philip D. Gingerich RATES OF EVOLUTION (Cambridge University Press 2019) 
##9.2.25_Chiba1996_Mandarina
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
file1<-"C://R_aaROEVchapt09//9.2.25_Chiba1996_Mandarina.csv" 
C<-read.csv(file1,header=TRUE,sep=",",quote="\"",dec=".",fill=TRUE)
C[1:5,];nrow(C);ncol(C)
#------------------------------------------------------Logged matrix
CC=apply(as.matrix.noquote(C[1:nrow(C),]),2,as.numeric)
CC[1:5,]
LC<-matrix(nrow=nrow(CC),ncol=1+10*3)
  colnames(LC)=c('Age','N1','X1','S1','N2','X2','S2','N3','X3','S3',
    'N4','X4','S4','N5','X5','S5','N6','X6','S6','N7','X7','S7',
    'N8','X8','S8','N9','X9','S9','N10','X10','S10')
LC[,1]=CC[,3]						#ages
LC[,seq(2,29,3)]=CC[,rep(4,10)]			#Ns
LC[,seq(3,30,3)]=log(CC[,seq(5,23,2)])		#means
LC[,seq(4,31,3)]=CC[,seq(6,24,2)]/CC[,seq(5,23,2)]	#stdev
LC[1:5,]
min(LC[,1]);max(LC[,1])
min(na.omit(LC[,3]))-min(na.omit(LC[,4]))
max(na.omit(LC[,3]))+max(na.omit(LC[,4]))

PanelA=FALSE;if(PanelA){
#======================================================= Plot panel a
xos=3;xoe=18;yos=9.5;yoe=16.5			#x,y original start and end
xfs=1.1;xfe=1.7;yfs=-56.2;yfe=53.8			#x,y foreground start and end
xfm=(xoe-xos)/(xfe-xfs)				#x foreground multiplier
yfm=(yoe-yos)/(yfe-yfs)				#y foreground multiplier
lines(c(xos,xoe),c(yos,yos),lwd=1,lty=1,col=1)
lines(c(xos,xos),c(yos,yoe-1),lwd=1,lty=1,col=1)
#---------------------------------------------------------axis labels
for (i in seq(yos+.5,yoe-1.2,.5)){				#y-axis
  lines(c(xos-.15,xos),c(i,i),lwd=1,lty=1,col=1)
}
for (i in seq(yos+.5,yoe-1.2,1)){
  text(xos+.1,i,format(round(abs(.4*i-9.9-50.1),digits=1),nsmall=1),pos=2,cex=1,col=1)
}
for (i in seq(xos,xoe,1)){		#x-axis
  lines(c(i,i),c(yos,yos-.15),lwd=1,lty=1,col=1)
}
#i=xos
for (i in seq(xos,xoe,2)){		#x-axis
  text(i+1,yos,format(round(.9+(i+1)/20,digits=1),nsmall=1),
    pos=1,cex=1,col=1)
}
text(xos,yos+6,expression(paste(italic('Cantius')~
  'spp.: molar'~M[1]~'length')),pos=4,cex=1.3,col=1)
text(xos-1.8,yos+3,'Age (Ma)',
	srt=90,cex=1.3,col=1)
text(xos+7.5,yos-1.2,
  expression(paste('Ln length (mm)')),cex=1.3,col=1)
#---------------------------------------------------------------- Plot ranges
LC[1:5,]
min(na.omit(LC[,3]));max(na.omit(LC[,3]))
min(LC[,1]);max(LC[,1])
xpos=-21;ypos=140.1
for (i in 1:nrow(LC)){
  lines(xos+20*c(LC[i,3]-LC[i,4],LC[i,3]+LC[i,4])+xpos,
    yos+.4+2.5*c(LC[i,1],LC[i,1])+ypos,lty=1,lwd=1,col=gray(6/10))
}
for (i in 1:nrow(LC)){
  points(xos+20*LC[i,3]+xpos,yos+.4+2.5*LC[i,1]+ypos,pch=20,cex=.8,col=1)
}
lines(xos+20*LC[,3]+xpos,yos+.4+2.5*LC[,1]+ypos,lty=1,lwd=1.5,col=1)
}

#============================================================= Rates: all
idr=matrix(nrow=0,ncol=9)
  colnames(idr)=c('int.g','diff.sd','rate.sd.g','log.i','log|d|','log|r|',
    'sbn','wgt','li+lr')
#------------------------------------------------ Rates
LC[1:3,]
gt=2	#generation time from Chiba 2007
nr=nrow(LC);nr					#number of rows
nn=.5*(nr-1)*nr;nn
zc=0							#zero count
nc=0
for (cn in seq(3,30,3)){			#columns of means
  for (i in 1:(nr-1)){				#increment
    for (sr in 1:(nr-i)){			#start row
      id<-numeric(length=9)
      nc=nc+1
      iyr=LC[sr,1]-LC[sr+i,1]				#interval in years
      #w1g=exp(1.495*(LC[sr,3]+LC[sr,6])+3.556)	#GSRbody weight in grams
      #w2g=exp(1.495*(LC[sr+1,3]+LC[sr+1,6])+3.556)#GSR body weight in grams
      #wg=exp(.5*(log(w1g)+log(w2g)))		#pooled weight
      #gt=10^(.266*log10(wg)-.553)			#Eq 3.2 gen. time in years
      id[1]=iyr/gt					#interval in generations
      meandiff=LC[sr+i,cn]-LC[sr,cn]		#mean diff.
      if (!is.na(meandiff)){if (meandiff==0){zc=zc+1}}
      psd=PoolSD(LC[sr+i,cn-1],LC[sr,cn-1],LC[sr+i,cn+1],LC[sr,cn+1])#n1,n2,sd1,sd2
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
}
#--------------------------------------------------------------------------
nrow(idr)
#--------------------------------------------------------------------------
idrx=idr[!rowSums(!is.finite(idr)),]  #remove rows that have -Inf
nn;nrow(idr);nrow(idrx);zc
min(idrx[,9]);max(idrx[,9])

writefile<-"C://R_aaROEVchapt09//9.2.25_Chiba1996_Mandarina_out.csv"
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
bootresultd=PalPanelBC(idrx[,1:8],"r",1,4.5,bootn,mode,psize,"normal",1)
proc.time()-ptm
#======================================================== Plot LRI panel (c)
proc.time()-ptm
#----------------------------------------------------------
text(7,-2.6,paste("Processing time: ",round(proc.time()[3]-ptm[3],digits=2),
  "seconds"),pos=4,cex=1,col=4)



