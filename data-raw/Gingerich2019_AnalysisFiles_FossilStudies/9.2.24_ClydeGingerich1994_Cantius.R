##Philip D. Gingerich RATES OF EVOLUTION (Cambridge University Press 2019)
##9.2.24_ClydeGingerich1994_Cantius
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
file1<-"..//Gingerich2019_AnalysisFiles_FossilStudies/9.2.24_ClydeGingerich1994_Cantius_Appendix01.csv"
C<-read.csv(file1,header=TRUE,sep=",",quote="\"",dec=".",fill=TRUE)
C[1:5,];nrow(C);ncol(C)
file2<-"../Gingerich2019_AnalysisFiles_FossilStudies/9.2.24_ClydeGingerich1994_Cantius_Appendix02.csv"
C2<-read.csv(file2,header=TRUE,sep=",",quote="\"",dec=".",fill=TRUE)
C2[1:5,];nrow(C2);ncol(C2)
#------------------------------------------------------Logged matrices
CC=apply(as.matrix.noquote(C[1:nrow(C),]),2,as.numeric)
CC[1:5,]
LC<-matrix(nrow=nrow(CC),ncol=1+6*3)
  colnames(LC)=c('Lev.','N1','X1','S1','N2','X2','S2','N3','X3','S3',
    'N4','X4','S4','N5','X5','S5','N6','X6','S6')
LC[,1]=-56+(CC[,1]-1520)/416				#ages
LC[,c(2,5,8,11,14,17)]=CC[,c(2,2,2,9,9,9)]	#Ns
LC[,c(3,6,9,12,15,18)]=CC[,c(3,5,7,10,12,14)]	#means
LC[,c(4,7,10,13,16,19)]=CC[,c(4,6,8,11,13,15)]	#stdev
CC[1:5,]
LC[1:5,]
min(LC[,1]);max(LC[,1])
min(na.omit(LC[,3]))-min(na.omit(LC[,4]))
max(na.omit(LC[,3]))+max(na.omit(LC[,4]))
#------------------------------------------------------
CC2=apply(as.matrix.noquote(C2[1:nrow(C2),]),2,as.numeric)
CC2[1:5,]
LC2<-matrix(nrow=nrow(CC2),ncol=1+7*3)
  colnames(LC2)=c('Lev.','N1','X1','S1','N2','X2','S2','N3','X3','S3',
    'N4','X4','S4','N5','X5','S5','N6','X6','S6','N7','X7','S7')
LC2[,1]=-56+(CC2[,1]-1520)/416				#ages
LC2[,c(2,5,8,11,14,17,20)]=CC2[,c(2,2,2,2,2,2,15)]	#Ns
LC2[,c(3,6,9,12,15,18,21)]=CC2[,c(3,5,7,9,11,13,16)]	#means
LC2[,c(4,7,10,13,16,19,22)]=CC2[,c(4,6,8,10,12,14,17)]	#stdev
CC2[1:5,]
LC2[1:5,]

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

text(5.9,10,expression(paste(italic('C. torresi'))),pos=2,cex=1,col=1)
text(9.9,10.5,expression(paste(italic('C. ralstoni'))),pos=4,cex=1,col=1)
text(12,12.4,expression(paste(italic('C. mckennai'))),pos=4,cex=1,col=1)
text(13.3,14.25,expression(paste(italic('C. trigonodus'))),pos=4,cex=1,col=1)

#============================================================= Rates: size1
idr=matrix(nrow=0,ncol=9)
  colnames(idr)=c('int.g','diff.sd','rate.sd.g','log.i','log|d|','log|r|',
    'sbn','wgt','li+lr')
#------------------------------------------------ Rates
LC[1:3,]
nr=nrow(LC);nr					#number of rows
nn1=.5*(nr-1)*nr;nn1
zc1=0							#zero count
nc=0
for (cn in c(3,6,12,15)){		#columns of means
  for (i in 1:(nr-1)){				#increment
    for (sr in 1:(nr-i)){			#start row
      id<-numeric(length=9)
      nc=nc+1
      iyr=1000000*(LC[sr+i,1]-LC[sr,1])		#interval in years
      w1g=exp(1.495*(LC[sr,3]+LC[sr,6])+3.556)	#GSRbody weight in grams
      w2g=exp(1.495*(LC[sr+1,3]+LC[sr+1,6])+3.556)#GSR body weight in grams
      wg=exp(.5*(log(w1g)+log(w2g)))		#pooled weight
      gt=10^(.266*log10(wg)-.553)			#Eq 3.2 gen. time in years
      id[1]=iyr/gt					#interval in generations
      meandiff=LC[sr+i,cn]-LC[sr,cn]		#mean diff.
      if (!is.na(meandiff)){if (meandiff==0){zc1=zc1+1}}
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
idrx1=idr[!rowSums(!is.finite(idr)),]  #remove rows that have -Inf
nn1;nrow(idr);nrow(idrx1);zc1
min(idrx1[,9]);max(idrx1[,9])

#====================================================== Rates: shape2
idr=matrix(nrow=0,ncol=9)
  colnames(idr)=c('int.g','diff.sd','rate.sd.g','log.i','log|d|','log|r|',
    'sbn','wgt','li+lr')
#------------------------------------------------ Rates
LC[1:3,]
nr=nrow(LC);nr					#number of rows
nn2=.5*(nr-1)*nr;nn2
zc2=0							#zero count
nc=0
for (cn in c(9,18)){				#columns of means
  for (i in 1:(nr-1)){				#increment
    for (sr in 1:(nr-i)){			#start row
      id<-numeric(length=9)
      nc=nc+1
      iyr=1000000*(LC[sr+i,1]-LC[sr,1])		#interval in years
      w1g=exp(1.495*(LC[sr,3]+LC[sr,6])+3.556)	#GSRbody weight in grams
      w2g=exp(1.495*(LC[sr+1,3]+LC[sr+1,6])+3.556)#GSR body weight in grams
      wg=exp(.5*(log(w1g)+log(w2g)))		#pooled weight
      gt=10^(.266*log10(wg)-.553)			#Eq 3.2 gen. time in years
      id[1]=iyr/gt					#interval in generations
      meandiff=LC[sr+i,cn]-LC[sr,cn]		#mean diff.
      if (!is.na(meandiff)){if (meandiff==0){zc2=zc2+1}}
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
idrx2=idr[!rowSums(!is.finite(idr)),]  #remove rows that have -Inf
nn2;nrow(idr);nrow(idrx2);zc2
min(idrx2[,9]);max(idrx2[,9])
#====================================================== Rates: shape3
idr=matrix(nrow=0,ncol=9)
  colnames(idr)=c('int.g','diff.sd','rate.sd.g','log.i','log|d|','log|r|',
    'sbn','wgt','li+lr')
#------------------------------------------------ Rates
LC2[1:3,]
nr=nrow(LC2);nr					#number of rows
nn3=.5*(nr-1)*nr;nn3
zc3=0							#zero count
nc=0
for (cn in c(3,6,9,12,15,18,21)){		#columns of means
  for (i in 1:(nr-1)){				#increment
    for (sr in 1:(nr-i)){			#start row
      id<-numeric(length=9)
      nc=nc+1
      iyr=1000000*(LC2[sr+i,1]-LC2[sr,1])		#interval in years
      w1g=exp(1.495*(LC2[sr,3]+LC2[sr,6])+3.556)	#GSRbody weight in grams
      w2g=exp(1.495*(LC2[sr+1,3]+LC2[sr+1,6])+3.556)#GSR body weight in grams
      wg=exp(.5*(log(w1g)+log(w2g)))		#pooled weight
      gt=10^(.266*log10(wg)-.553)			#Eq 3.2 gen. time in years
      id[1]=iyr/gt					#interval in generations
      meandiff=LC2[sr+i,cn]-LC2[sr,cn]		#mean diff.
      if (!is.na(meandiff)){if (meandiff==0){zc3=zc3+1}}
      psd=PoolSD(LC2[sr+i,cn-1],LC2[sr,cn-1],LC2[sr+i,cn+1],LC2[sr,cn+1])#n1,n2,sd1,sd2
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
idrx3=idr[!rowSums(!is.finite(idr)),]  #remove rows that have -Inf
nn3;nrow(idr);nrow(idrx3);zc3
min(idrx3[,9]);max(idrx3[,9])

#--------------------------------------------------------------------Summary
idrxa=rbind(idrx1,idrx2,idrx3);idrx23=rbind(idrx2,idrx3)
nrow(idrx1);nrow(idrx2);nrow(idrx3);nrow(idrx23);nrow(idrxa)
zc1;zc2;zc3;zc1+zc2+zc3
zca=zc1+zc2+zc3
zca;nrow(idrxa);zca+nrow(idrxa);nn1+nn2+nn3
writefile<-"C://R_aaROEVchapt09//9.2.24_ClydeGingerich1994_Cantius_out.csv"
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
bootresultd=PalPanelBC(idrx1[,1:8],"r",1,4.5,bootn,mode,psize,"normal",1)
proc.time()-ptm
#======================================================== Plot LRI panel (c)
#send idrx matrix, n, mode(diff/rate), panel placement coordinates, mode
bootresultr=PalPanelBC(idrx23[,1:8],"r",13,4.5,bootn,mode,psize,"normal",0)
proc.time()-ptm

text(9.2,4,expression(paste("Size change")),pos=2,cex=1,col=1)
text(21.2,4,expression(paste("Shape change")),pos=2,cex=1,col=1)

#----------------------------------------------------------
text(7,-2.6,paste("Processing time: ",round(proc.time()[3]-ptm[3],digits=2),
  "seconds"),pos=4,cex=1,col=4)



