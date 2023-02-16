##Philip D. Gingerich RATES OF EVOLUTION (Cambridge University Press 2019)
##9.2.17_Hayami1984_Cryptopecten
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
file1<-"../Gingerich2019_AnalysisFiles_FossilStudies/9.2.17_Hayami1984_Cryptopecten.csv"
H<-read.csv(file1,header=TRUE,sep=",",quote="\"",dec=".",fill=TRUE)
nrow(H);H
#------------------------------------------------------Logged matrix
LH<-matrix(nrow=nrow(H),ncol=7)
  colnames(LH)=c('Age','N1','M1','S1','N2','M2','S2')
LH[,1]=H[,2]
for (i in c(2,5)){
  LH[,i]=H[,i+1]
  LH[,i+1]=log(H[,i+2])
  LH[,i+2]=.01*H[,i+3]
};LH
PanelA=TRUE #=========================================== Plot panel a
if (PanelA==TRUE){
xos=2;xoe=20;yos=9.5;yoe=17.5			#x,y original start and end
xfs=6.1;xfe=7.9;yfs=-160;yfe=-0			#x,y foreground start and end
xfm=(xoe-xos)/(xfe-xfs)				#x foreground multiplier
yfm=(yoe-yos)/(yfe-yfs)				#y foreground multiplier
lines(c(xos,xoe),c(yos,yos),lwd=1,lty=1,col=1)
lines(c(xos,xos),c(yos,yoe),lwd=1,lty=1,col=1)
#---------------------------------------------------------axis labels
for (i in seq(yos,yoe-1.5,1)){				#y-axis
  lines(c(xos-.15,xos),c(i,i),lwd=1,lty=1,col=1)
}
for (i in seq(yos,yoe-1.5,1)){
  text(xos+.1,i,format(round(abs(375-25*(i-.5)),digits=1),nsmall=0),pos=2,cex=1,col=1)
}

for (i in xos+10*(seq(xfs,xfe,.1)-xfs)){		#x-axis
  lines(c(i,i),c(yos,yos-.15),lwd=1,lty=1,col=1)
}
for (i in xos+10*(seq(xfs,xfe,.1)-xfs)){
  text(i,yos,format(round(5.9+.1*i,digits=1),nsmall=1),
    pos=1,cex=1,col=1)
}
text(xos,yos+8,expression(paste(italic('Bison')~
	'spp.:  width across horn cores')),pos=4,cex=1.3,col=1)
text(xos-1.8,yos+3,'Geological age (kyr)',
	srt=90,cex=1.3,col=1)
text(xos+7.5,yos-1.2,
  expression(paste('Ln measurement (mm)')),cex=1.3,col=1)
#---------------------------------------------------------------- Plot ranges
xpos=xos-61;ypos=yos+6
#------------------------------------------------------- mean lines
lines(xpos+10*c(LM[1,3],LM[1,7]),
  ypos-c(LM[1,1],LM[1,5])/25000,lty=1,lwd=2,col=1)
lines(xpos+10*c(LM[1,c(3,11,15)],LM[1,c(11,15,19)]),
  ypos-c(LM[1,c(1,9,13)],LM[1,c(9,13,17)])/25000,lty=1,lwd=2,col=1)
namexy<-matrix(nrow=5,ncol=2)
for (i in 5:1){
  xbar=xpos+10*LM[1,4*i-1]
  stdev=10*LM[1,4*i]
  base=ypos-LM[1,4*i-3]/25000
  hgt=.8
  DrawNorm(xbar,stdev,base,hgt)		#mean,sd,base,hgt
  namexy[i,1]=xbar;namexy[i,2]=base+1.3
}
text(namexy[1,1],namexy[1,2],expression(paste(italic('B. priscus'))),
  pos=4,cex=1,col=1)
text(namexy[2,1],namexy[2,2],expression(paste(italic('B. latifrons'))),
  pos=4,cex=1,col=1)
text(namexy[3,1]+.4,namexy[3,2]-.3,expression(paste(italic('B. antiquus'))),
  pos=4,cex=1,col=1)
text(namexy[4,1],namexy[4,2],expression(paste(italic('B. occidentalis'))),
  pos=4,cex=1,col=1)
text(namexy[5,1],namexy[5,2],expression(paste(italic('B. bison'))),
  pos=4,cex=1,col=1)
rect(namexy[1,1]-2,namexy[1,2]-1.5,namexy[1,1]+2,namexy[1,2]-1.29,
  col='white',lwd=1,border=NA)
}
#=============================================================== Rates
idr=matrix(nrow=0,ncol=9)
  colnames(idr)=c('int.g','diff.sd','rate.sd.g','log.i','log|d|','log|r|',
    'sbn','wgt','li+lr')
#------------------------------------------------------
LH[1:5,]
gt=3	#generation time as median from age of reproduction and max age
nr=nrow(LH);nr					#number of rows
nn=.5*(nr-1)*nr;nn
zc=0							#zero count
nc=0
for (i in 1:(nr-1)){				#increment
  for (sr in 1:(nr-i)){			#start row
    for (cn in seq(3,3,3)){		#columns of means
      id<-numeric(length=9)
      nc=nc+1
      iyr=LH[sr,1]-LH[sr+i,1]	#interval in years
      id[1]=iyr/gt					#interval in generations
      meandiff=LH[sr,cn]-LH[sr+i,cn]		#mean diff.
      if (!is.na(meandiff)){if (meandiff==0){zc=zc+1}}
      psd=PoolSD(LH[sr+i,cn-1],LH[sr,cn-1],LH[sr+i,cn+1],LH[sr,cn+1])#n1,n2,sd1,sd2
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
idrs=idr[!rowSums(!is.finite(idr)),]  #remove rows that have -Inf
nn;nrow(idr);nrow(idrs);zc
min(idrs[,9]);max(idrs[,9])
#-------------------------------------------------------------------------
idr=matrix(nrow=0,ncol=9)
  colnames(idr)=c('int.g','diff.sd','rate.sd.g','log.i','log|d|','log|r|',
    'sbn','wgt','li+lr')
#------------------------------------------------------
LH[1:5,]
gt=3	#generation time as median from age of reproduction and max age
nr=nrow(LH);nr					#number of rows
nn=.5*(nr-1)*nr;nn
zc=0							#zero count
nc=0
for (i in 1:(nr-1)){				#increment
  for (sr in 1:(nr-i)){			#start row
    for (cn in seq(6,6,6)){		#columns of means
      id<-numeric(length=9)
      nc=nc+1
      iyr=LH[sr,1]-LH[sr+i,1]	#interval in years
      id[1]=iyr/gt					#interval in generations
      meandiff=LH[sr,cn]-LH[sr+i,cn]		#mean diff.
      if (!is.na(meandiff)){if (meandiff==0){zc=zc+1}}
      psd=PoolSD(LH[sr+i,cn-1],LH[sr,cn-1],LH[sr+i,cn+1],LH[sr,cn+1])#n1,n2,sd1,sd2
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
idrr=idr[!rowSums(!is.finite(idr)),]  #remove rows that have -Inf
nn;nrow(idr);nrow(idrr);zc
min(idrr[,9]);max(idrr[,9])

idrx=rbind(idrs,idrr)
#--------------------------------------------------------------------------
writefile<-"C://R_aaROEVchapt09//9.2.17_Hayami1984_Cryptopecten_out.csv"
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
bootresultd=PalPanelBC(idrs[,1:8],"r",1,4.5,bootn,mode,psize,"normal",1)
proc.time()-ptm
#======================================================== Plot LRI panel (c)
#send idrx matrix, n, mode(diff/rate), panel placement coordinates, mode
bootresultr=PalPanelBC(idrr[,1:8],"r",13,4.5,bootn,mode,psize,"normal",1)
proc.time()-ptm
#----------------------------------------------------------
text(7,-2.6,paste("Processing time: ",round(proc.time()[3]-ptm[3],digits=2),
  "seconds"),pos=4,cex=1,col=4)


