##Philip D. Gingerich RATES OF EVOLUTION (Cambridge University Press 2019) 
##9.2.26_Gingerich2003_Ectocion
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
file1<-"C://R_aaROEVchapt09//9.2.26_Gingerich2003_Ectocion.csv" 
G<-read.csv(file1,header=TRUE,sep=",",quote="\"",dec=".",fill=TRUE)
G[1:5,];nrow(G);ncol(G)
#------------------------------------------------------Logged matrices
GG=apply(as.matrix.noquote(G[1:nrow(G),]),2,as.numeric)
GG[1:5,]
SG<-matrix(nrow=0,ncol=4);colnames(SG)=c('age','N','lnM','lnS')	#Stat G
st=GG[1,1];nc=0		#strat level
samp<-matrix(nrow=0,ncol=3);colnames(samp)=c('age','n','lnlxw')
for (i in 1:nrow(GG)){
  if (GG[i,1]==st){
    nc=nc+1;samp=rbind(samp,c(st,nc,log(GG[i,2])+log(GG[i,3])))
    st=GG[i,1]
  }else{
    if (nc>1){ 
      sg<-numeric(length=4)
      sg[1]=st;sg[2]=nrow(samp);sg[3]=mean(samp[,3]);sg[4]=sd(samp[,3])
      SG=rbind(SG,sg)
    }
    if (nc==1){
      SG=rbind(SG,c(samp,NA))
    }
    samp<-matrix(nrow=0,ncol=3);colnames(samp)=c('age','n','lnlxw')
    nc=1;samp=c(GG[i,1],1,log(GG[i,2])+log(GG[i,3]))
    st=GG[i,1]
  }
}
SG=rbind(SG,c(samp,NA));SG
#---------------------------------------------------------- pool small sd
nrow(SG)
redSG=SG[!is.na(SG[,4]),];nrow(redSG)
ssd=sqrt(sum(redSG[,2]*redSG[,4]^2)/sum(redSG[,2]))
SG[SG[,2]<=3,4]=ssd  #SG
#---------------------------------------------------------- convert ages
SG[,1]=.0025*SG[,1]-59.777
min(SG[,1]);max(SG[,1])

#======================================================= Plot panel a
xos=3;xoe=18;yos=9.5;yoe=16.5			#x,y original start and end
xfs=2.8;xfe=4.3;yfs=-57.6;yfe=54.6			#x,y foreground start and end
xfm=(xoe-xos)/(xfe-xfs)				#x foreground multiplier
yfm=(yoe-yos)/(yfe-yfs)				#y foreground multiplier
lines(c(xos,xoe),c(yos,yos),lwd=1,lty=1,col=1)
lines(c(xos,xos),c(yos,yoe),lwd=1,lty=1,col=1)
#---------------------------------------------------------axis labels
for (i in seq(yos,yoe,1)){				#y-axis
  lines(c(xos-.15,xos),c(i,i),lwd=1,lty=1,col=1)
}
for (i in seq(yos,yoe,1)){
  text(xos+.1,i,format(round(abs(.5*(i-9.5)-58),digits=1),nsmall=1),pos=2,cex=1,col=1)
}
for (i in seq(xos,xoe,1)){		#x-axis
  lines(c(i,i),c(yos,yos-.15),lwd=1,lty=1,col=1)
}
for (i in seq(xos+1,xoe,2)){		#x-axis
  text(i,yos,format(round(2.4+i/10,digits=1),nsmall=1),
    pos=1,cex=1,col=1)
}
text(xos,yos+7.5,expression(paste(italic('Ectocion')~
  'spp.: molar'~M[1]~'crown area')),pos=4,cex=1.3,col=1)
text(xos-1.8,yos+3.5,'Age (Ma)',
	srt=90,cex=1.3,col=1)
text(xos+7.5,yos-1.2,
  expression(paste('Ln (L x W) of'~M[1]~'(mm)')),cex=1.3,col=1)
#---------------------------------------------------------------- Plot ranges
SG[1:5,]
min(na.omit(SG[,3]));max(na.omit(SG[,3]))
min(SG[,1]);max(SG[,1])
xpos=-27;ypos=116
for (i in 1:nrow(SG)){
  lines(xos+10*c(SG[i,3]-SG[i,4],SG[i,3]+SG[i,4])+xpos,
    yos+2*c(SG[i,1],SG[i,1])+ypos,lty=1,lwd=1,col=gray(6/10))
}
for (i in 1:nrow(SG)){
  points(xos+10*SG[i,3]+xpos,yos+2*SG[i,1]+ypos,pch=20,cex=.8,col=1)
}
lines(xos+10*SG[,3]+xpos,yos+2*SG[,1]+ypos,lty=1,lwd=1.5,col=1)

text(6.7,13.52,expression(paste(italic('E. parvus'))),pos=2,cex=1,col=1)
text(13.3,15.2,expression(paste(italic('E. osbornianus'))),pos=4,cex=1,col=1)
text(14.1,12.3,expression(paste(italic('E. ralstonensis'))),pos=4,cex=1,col=1)
lines(c(13.8,17.5),c(13.53,13.53),lty=2,lwd=1,col=1)
lines(c(3,3.6),c(13.53,13.53),lty=2,lwd=1,col=1)
text(17.5,13.53,expression(paste('PETM')),pos=4,cex=1,col=1)

#=============================================== Species samples
SP1=SG[1:54,];nrow(SP1)		#Clarkforkian Ectocion ralstonensis
SP12=SG[54:55,];SP12
SP23=SG[55:56,];SP23
SP3=SG[56:78,];nrow(SP3)
nrow(SG)
#=============================================== Rates: E. ralstonensis
idr=matrix(nrow=0,ncol=9)
  colnames(idr)=c('int.g','diff.sd','rate.sd.g','log.i','log|d|','log|r|',
    'sbn','wgt','li+lr')
#------------------------------------------------ Rates
SP1[1:3,]
nr=nrow(SP1);nr						#number of rows
nn1=.5*(nr-1)*nr;nn1
zc1=0							#zero count
nc=0
cn=3
  for (i in 1:(nr-1)){				#increment
    for (sr in 1:(nr-i)){			#start row
      id<-numeric(length=9)
      nc=nc+1
      iyr=1000000*(SP1[sr+i,1]-SP1[sr,1])		#interval in years
      w1g=exp(1.5147*SP1[sr,cn]+3.6049)		#Ongul. regr. of Legendre 1989
      w2g=exp(1.5147*SP1[sr+1,cn]+3.6049)		#GSR body weight in grams
      wg=exp(.5*(log(w1g)+log(w2g)))		#pooled weight
      gt=10^(.266*log10(wg)-.553)			#Eq 3.2 gen. time in years
      id[1]=iyr/gt					#interval in generations
      meandiff=SP1[sr+i,cn]-SP1[sr,cn]		#mean diff.
      if (!is.na(meandiff)){if (meandiff==0){zc1=zc1+1}}
      psd=PoolSD(SP1[sr+i,cn-1],SP1[sr,cn-1],SP1[sr+i,cn+1],SP1[sr,cn+1])#n1,n2,sd1,sd2
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
idrx1=idr[!rowSums(!is.finite(idr)),]  #remove rows that have -Inf
nn1;nrow(idr);nrow(idrx1);zc1
min(idrx1[,9]);max(idrx1[,9])
#=============================================== Rates: E. ralst.-E. parv.
idr=matrix(nrow=0,ncol=9)
  colnames(idr)=c('int.g','diff.sd','rate.sd.g','log.i','log|d|','log|r|',
    'sbn','wgt','li+lr')
#------------------------------------------------ Rates
SP12[1:3,]
nr=nrow(SP12);nr						#number of rows
nn2=.5*(nr-1)*nr;nn2
zc2=0							#zero count
nc=0
cn=3
  for (i in 1:(nr-1)){				#increment
    for (sr in 1:(nr-i)){			#start row
      id<-numeric(length=9)
      nc=nc+1
      iyr=1000000*(SP12[sr+i,1]-SP12[sr,1])		#interval in years
      w1g=exp(1.5147*SP12[sr,cn]+3.6049)		#Ongul. regr. of Legendre 1989
      w2g=exp(1.5147*SP12[sr+1,cn]+3.6049)		#GSR body weight in grams
      wg=exp(.5*(log(w1g)+log(w2g)))		#pooled weight
      gt=10^(.266*log10(wg)-.553)			#Eq 3.2 gen. time in years
      id[1]=iyr/gt					#interval in generations
      meandiff=SP12[sr+i,cn]-SP12[sr,cn]		#mean diff.
      if (!is.na(meandiff)){if (meandiff==0){zc2=zc2+1}}
      psd=PoolSD(SP12[sr+i,cn-1],SP12[sr,cn-1],SP12[sr+i,cn+1],SP12[sr,cn+1])#n1,n2,sd1,sd2
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
#=============================================== Rates: E. parv.-E. osborn.
idr=matrix(nrow=0,ncol=9)
  colnames(idr)=c('int.g','diff.sd','rate.sd.g','log.i','log|d|','log|r|',
    'sbn','wgt','li+lr')
#------------------------------------------------ Rates
SP23[1:3,]
nr=nrow(SP23);nr						#number of rows
nn3=.5*(nr-1)*nr;nn3
zc3=0							#zero count
nc=0
cn=3
  for (i in 1:(nr-1)){				#increment
    for (sr in 1:(nr-i)){			#start row
      id<-numeric(length=9)
      nc=nc+1
      iyr=1000000*(SP23[sr+i,1]-SP23[sr,1])		#interval in years
      w1g=exp(1.5147*SP23[sr,cn]+3.6049)		#Ongul. regr. of Legendre 1989
      w2g=exp(1.5147*SP23[sr+1,cn]+3.6049)		#GSR body weight in grams
      wg=exp(.5*(log(w1g)+log(w2g)))		#pooled weight
      gt=10^(.266*log10(wg)-.553)			#Eq 3.2 gen. time in years
      id[1]=iyr/gt					#interval in generations
      meandiff=SP23[sr+i,cn]-SP23[sr,cn]		#mean diff.
      if (!is.na(meandiff)){if (meandiff==0){zc3=zc3+1}}
      psd=PoolSD(SP23[sr+i,cn-1],SP23[sr,cn-1],SP23[sr+i,cn+1],SP23[sr,cn+1])#n1,n2,sd1,sd2
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
idrx3=idr[!rowSums(!is.finite(idr)),]  #remove rows that have -Inf
nn3;nrow(idr);nrow(idrx3);zc3
min(idrx3[,9]);max(idrx3[,9])
#=============================================== Rates: E. osbornianus
idr=matrix(nrow=0,ncol=9)
  colnames(idr)=c('int.g','diff.sd','rate.sd.g','log.i','log|d|','log|r|',
    'sbn','wgt','li+lr')
#------------------------------------------------ Rates
SP3[1:3,]
nr=nrow(SP3);nr						#number of rows
nn4=.5*(nr-1)*nr;nn4
zc4=0							#zero count
nc=0
cn=3
  for (i in 1:(nr-1)){				#increment
    for (sr in 1:(nr-i)){			#start row
      id<-numeric(length=9)
      nc=nc+1
      iyr=1000000*(SP3[sr+i,1]-SP3[sr,1])		#interval in years
      w1g=exp(1.5147*SP3[sr,cn]+3.6049)		#Ongul. regr. of Legendre 1989
      w2g=exp(1.5147*SP3[sr+1,cn]+3.6049)		#GSR body weight in grams
      wg=exp(.5*(log(w1g)+log(w2g)))		#pooled weight
      gt=10^(.266*log10(wg)-.553)			#Eq 3.2 gen. time in years
      id[1]=iyr/gt					#interval in generations
      meandiff=SP3[sr+i,cn]-SP3[sr,cn]		#mean diff.
      if (!is.na(meandiff)){if (meandiff==0){zc4=zc4+1}}
      psd=PoolSD(SP3[sr+i,cn-1],SP3[sr,cn-1],SP3[sr+i,cn+1],SP3[sr,cn+1])#n1,n2,sd1,sd2
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
idrx4=idr[!rowSums(!is.finite(idr)),]  #remove rows that have -Inf
nn4;nrow(idr);nrow(idrx4);zc4
min(idrx4[,9]);max(idrx4[,9])

#--------------------------------------------------------------------Summary
idrxa=rbind(idrx1,idrx2,idrx3,idrx4)
nrow(idrx1);nrow(idrx2);nrow(idrx3);nrow(idrx4);nrow(idrxa)
zc1;zc2;zc3;zc4;zca=zc1+zc2+zc3+zc4
zca;nrow(idrxa);zca+nrow(idrxa);nn1+nn2+nn3+nn4
writefile<-"C://R_aaROEVchapt09//9.2.26_Gingerich2003_Ectocion_out.csv"
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
#send idrx matrix, n, mode(diff/rate), panel placement coordinates, mode
bootresultr=PalPanelBC(idrx4[,1:8],"r",13,4.5,bootn,mode,psize,"normal",0)
proc.time()-ptm

text(9.2,5,expression(paste(italic("E. ralstonensis"))),pos=2,cex=1,col=1)
text(21.2,5,expression(paste(italic("E. osbornianus"))),pos=2,cex=1,col=1)

lines(c(idrx3[4]+1,idrx3[4]+1.5),c(idrx3[6]+6.5,idrx3[6]+7),
  lty=1,lwd=2,col='white')
lines(c(idrx3[4]+1,idrx3[4]+1.5),c(idrx3[6]+6.5,idrx3[6]+7),
  lty=1,lwd=1,col=1)
points(idrx2[4]+1,idrx2[6]+6.5,pch=19,cex=1.1,col=1)
points(idrx3[4]+1,idrx3[6]+6.5,pch=19,cex=1.1,col=1)
text(idrx3[4]+1.5-.6,idrx3[6]+7+.3,'transition rates',pos=4,cex=1,col=1)

lines(c(idrx3[4]+13,idrx3[4]+13.5),c(idrx3[6]+6.5,idrx3[6]+7),
  lty=1,lwd=2,col='white')
lines(c(idrx3[4]+13,idrx3[4]+13.5),c(idrx3[6]+6.5,idrx3[6]+7),
  lty=1,lwd=1,col=1)
points(idrx2[4]+13,idrx2[6]+6.5,pch=19,cex=1.1,col=1)
points(idrx3[4]+13,idrx3[6]+6.5,pch=19,cex=1.1,col=1)
text(idrx3[4]+13.5-.6,idrx3[6]+7+.3,'transition rates',pos=4,cex=1,col=1)

#----------------------------------------------------------
text(7,-2.6,paste("Processing time: ",round(proc.time()[3]-ptm[3],digits=2),
  "seconds"),pos=4,cex=1,col=4)



