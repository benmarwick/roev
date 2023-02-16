##Philip D. Gingerich RATES OF EVOLUTION (Cambridge University Press 2019)
##9.2.16_Malmgren&al1983_Globorotalia_size
#==========================================================================
cat(rep("\n",50))							#clear console
print (date()) #Ctr-s to save, Ctr-a to select all, Ctr-r to run
rm(list=ls(all=TRUE))#remove/clear all prev. variables
assign('last.warning',NULL,envir=baseenv())
ptm<-proc.time()
##======================================================================##
##library("devtools");library(roxygen2)
##setwd("c:/R_aaPackages/ROEV");document();setwd("..");install("ROEV")
#=====================================================================Setup
library(ROEV);library(MASS)
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
file1<-"C://R_aaROEVchapt09//9.2.16_Malmgren&al1983_Globorotalia.csv"
M<-read.csv(file1,header=TRUE,sep=",",quote="\"",dec=".",fill=TRUE)
M[1:5,];nrow(M)
#------------------------------------------------------Logged matrix
MM=apply(as.matrix.noquote(M[1:nrow(M),]),2,as.numeric)
MM[1:5,]
LM<-matrix(nrow=nrow(MM),ncol=7)
  colnames(LM)=c('age','Nsize','Xsize','Ssize','Nshap','Xshap','Sshap')
LM[,1]=MM[,1]
LM[,c(2,5)]=MM[,2]
LM[,3]=log(MM[,3])
LM[,4]=MM[,4]/MM[,3]
LM[,6]=MM[,5]
LM[,7]=MM[,6]

MM[1:5,]
LM[1:5,]
min(LM[,1]);max(LM[,1])
min(LM[,3]);max(LM[,3])

#======================================================= Plot panel a
xos=3;xoe=18;yos=9.5;yoe=17.5			#x,y original start and end
xfs=4;xfe=5.4;yfs=-7;yfe=1			#x,y foreground start and end
xfm=(xoe-xos)/(xfe-xfs)				#x foreground multiplier
yfm=(yoe-yos)/(yfe-yfs)				#y foreground multiplier
lines(c(xos,xoe),c(yos,yos),lwd=1,lty=1,col=1)
lines(c(xos,xos),c(yos,yoe),lwd=1,lty=1,col=1)
lines(c(xos,xoe),c(11.1,11.1),lwd=1,lty=2,col=1)
text(xos,11.5,'Pliocene',pos=4,cex=1,col=1)
text(xos,10.7,'Mio.',pos=4,cex=1,col=1)
#---------------------------------------------------------axis labels
for (i in seq(yos,yoe-.5,1)){				#y-axis
  lines(c(xos-.15,xos),c(i,i),lwd=1,lty=1,col=1)
}
for (i in seq(yos,yoe-.5,1)){
  text(xos+.1,i,format(round(abs(16-(i-.5)),digits=1),nsmall=1),pos=2,cex=1,col=1)
}
for (i in xos+10*(seq(xfs,xfe,.1)-xfs)){		#x-axis
  lines(c(i,i),c(yos,yos-.15),lwd=1,lty=1,col=1)
}
for (i in xos+10*(seq(xfs,xfe,.1)-xfs)){
  text(i,yos,format(round(-4.2+.2*i,digits=1),nsmall=1),
    pos=1,cex=1,col=1)
}
text(xos,yos+8.0,expression(paste(italic('Globorotalia tumida')~
	'lineage:  size')),pos=4,cex=1.3,col=1)
text(xos-1.5,yos+4,'Geological age (Ma)',
	srt=90,cex=1.3,col=1)
text(xos+7.5,yos-1.2,
  expression(paste('Ln cross-sectional area ('*mm^2*')')),cex=1.3,col=1)
#---------------------------------------------------------------- Plot ranges
LM[1:5,]
min(na.omit(LM[,3]));max(na.omit(LM[,3]))
min(LM[,1]);max(LM[,1])
xpos=18;ypos=25
for (i in 1:nrow(LM)){
  lines(xos+5*c(LM[i,3]-LM[i,4],LM[i,3]+LM[i,4])+xpos,
    ypos-c(yos-1+LM[i,1],yos-1+LM[i,1]),lty=1,lwd=1,col=gray(6/10))
}
for (i in 1:nrow(LM)){
  points(xos+5*LM[i,3]+xpos,ypos-(yos-1+LM[i,1]),pch=20,cex=.8,col=1)
}
lines(xos+5*LM[,3]+xpos,ypos-(yos-1+LM[,1]),lty=1,lwd=1.5,col=1)

text(14,12,expression(paste(italic('G. tumida'))),pos=4,cex=1,col=1)
text(8.6,10.2,expression(paste(italic('G. plesiotumida'))),pos=4,cex=1,col=1)

#==========================================================================
SP1=LM[87:95,];SP1
SPt=LM[44:87,];SPt
SP2=LM[1:44,];SP2
#================================================ Rates SP1 G. plesiotumida
idr=matrix(nrow=0,ncol=9)
  colnames(idr)=c('int.g','diff.sd','rate.sd.g','log.i','log|d|','log|r|',
    'sbn','wgt','li+lr')
#------------------------------------------------ Rates
SP1[1:3,]
gt=1	#generation time assumed
nr=nrow(SP1)						#number of rows
nn1=.5*(nr-1)*nr;nn1
nz1=0							#zero n count
nc=0
cn=3
  for (i in 1:(nr-1)){				#increment
    for (sr in 1:(nr-i)){			#start row
      id<-numeric(length=9)
      nc=nc+1
      iyr=1000000*(SP1[sr+i,1]-SP1[sr,1])	#interval in years
      #w1g=1000*LZ[sr,5];w2g=1000*LZ[(sr+i),5];wg=exp(.5*(log(w1g)+log(w2g)))
      #10^(.266*log10(wg)-.553)		#Eq 3.2 generation time in years
      id[1]=iyr/gt					#interval in generations
      meandiff=SP1[sr,cn]-SP1[sr+i,cn]		#mean diff.
      if (!is.na(meandiff)){if (meandiff==0){nz1=nz1+1}}
      psd=PoolSD(SP1[sr+i,cn-1],SP1[sr,cn-1],SP1[sr+i,cn+1],SP1[sr,cn+1])#n1,n2,sd1,sd2
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
nn1;nrow(idr);nrow(idrx1);nz1
min(idrx1[,9]);max(idrx1[,9])

#====================================================== Rates in transition
idr=matrix(nrow=0,ncol=9)
  colnames(idr)=c('int.g','diff.sd','rate.sd.g','log.i','log|d|','log|r|',
    'sbn','wgt','li+lr')
#------------------------------------------------------ Rates
gt=1	#generation time assumed
nr=nrow(SPt)-1					#number of rows
nnt=.5*(nr-1)*nr;nnt
nzt=0							#zero n count
nc=0
cn=3
  for (i in 1:(nr-1)){				#increment
    for (sr in 1:(nr-i)){			#start row
      id<-numeric(length=9)
      nc=nc+1
      iyr=1000000*(SPt[sr+i,1]-SPt[sr,1])	#interval in years
      #w1g=1000*LZ[sr,5];w2g=1000*LZ[(sr+i),5];wg=exp(.5*(log(w1g)+log(w2g)))
      #10^(.266*log10(wg)-.553)		#Eq 3.2 generation time in years
      id[1]=iyr/gt					#interval in generations
      meandiff=SPt[sr,cn]-SPt[sr+i,cn]		#mean diff.
      if (!is.na(meandiff)){if (meandiff==0){nzt=nzt+1}}
      psd=PoolSD(SPt[sr+i,cn-1],SPt[sr,cn-1],SPt[sr+i,cn+1],SPt[sr,cn+1])#n1,n2,sd1,sd2
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
idrxt=idr[!rowSums(!is.finite(idr)),]  #remove rows that have -Inf
nnt;nrow(idr);nrow(idrxt);nzt
min(na.omit(idrxt[,9]));na.omit(max(idrxt[,9]))

#================================================ Rates SP2 G. tumida
idr=matrix(nrow=0,ncol=9)
  colnames(idr)=c('int.g','diff.sd','rate.sd.g','log.i','log|d|','log|r|',
    'sbn','wgt','li+lr')
#------------------------------------------------ Rates
gt=1	#generation time assumed
nr=nrow(SP2)						#number of rows
nn2=.5*(nr-1)*nr;nn2
nz2=0							#zero n count
nc=0
cn=3
  for (i in 1:(nr-1)){				#increment
    for (sr in 1:(nr-i)){			#start row
      id<-numeric(length=9)
      nc=nc+1
      iyr=1000000*(SP2[sr+i,1]-SP2[sr,1])	#interval in years
      #w1g=1000*LZ[sr,5];w2g=1000*LZ[(sr+i),5];wg=exp(.5*(log(w1g)+log(w2g)))
      #10^(.266*log10(wg)-.553)		#Eq 3.2 generation time in years
      id[1]=iyr/gt					#interval in generations
      meandiff=SP2[sr,cn]-SP2[sr+i,cn]		#mean diff.
      if (!is.na(meandiff)){if (meandiff==0){nz2=nz2+1}}
      psd=PoolSD(SP2[sr+i,cn-1],SP2[sr,cn-1],SP2[sr+i,cn+1],SP2[sr,cn+1])#n1,n2,sd1,sd2
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
idrx2=idr[!rowSums(!is.finite(idr)),]  #remove rows that have -Inf
nn2;nrow(idr);nrow(idrx2);nz2
min(idrx2[,9]);max(idrx2[,9])

#--------------------------------------------------------------------Summary
idrxa=rbind(idrx1,idrxt,idrx2);idrx12=rbind(idrx1,idrx2)
nrow(idrx1);nrow(idrxt);nrow(idrx2);nrow(idrx12);nrow(idrxa)
nz1;nzt;nz2;nz1+nzt+nz2
nza=nz1+nzt+nz2
nza;nrow(idrxa);nza+nrow(idrxa);nn1+nnt+nn2
writefile<-"C://R_aaROEVchapt09//9.2.16_Malmgren&al1983_Globorotalia_size_out.csv"
write.csv(idrxa,file=writefile,na="NA",row.names=FALSE)

#======================================================== Plot LRI panel (b)
idrx12=rbind(idrx1,idrx2)
#assign 'bootn' as boot number
bootn<-1000;text(-1.7,-2.6,paste("Boot n =",bootn),pos=4,cex=1,col=4)
#assign 'mode' as "medians","all","mixed"
mode<-"all";text(2,-2.6,paste("Mode: ",mode),pos=4,cex=1,col=4)
#assign circle size for points (1.5 or 2)
psize<-1.2	#1.5/2
#assign 'equation' position as "normal","lower","none" at end of each call
#send (1)idrx matrix, (2)mode(diff/rate), (3)panel placement coordinate x,
#  (4)panel placement coordinate y, (5)bootn, (6)mode, (7)psize, (8)equation
bootresultd=PalPanelBC(idrx12[,1:8],"r",1,4.5,bootn,mode,psize,"normal",1)
proc.time()-ptm
#======================================================== Plot LRI panel (c)
#send idrx matrix, n, mode(diff/rate), panel placement coordinates, mode
bootresultr=PalPanelBC(idrxt[,1:8],"r",13,4.5,bootn,mode,psize,"normal",1)
proc.time()-ptm

text(7,4.6,expression(paste(italic("G. plesio. / tumida")~'rates')),pos=1,cex=1,col=1)
text(19,4.6,expression(paste("Transition rates")),pos=1,cex=1,col=1)

#----------------------------------------------------------
text(7,-2.6,paste("Processing time: ",round(proc.time()[3]-ptm[3],digits=2),
  "seconds"),pos=4,cex=1,col=4)



