##Philip D. Gingerich RATES OF EVOLUTION (Cambridge University Press 2019) 
##8.1.17_Milner&al1999_Ovisaries_1wgt
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
xos=3.5;xoe=17;yos=9;yoe=17			#x,y original start and end
xfs=-1;xfe=12;yfs=2.4;yfe=4.0			#x,y foreground start and end
xfm=(xoe-xos)/(xfe-xfs)				#x foreground multiplier
yfm=(yoe-yos)/(yfe-yfs)				#y foreground multiplier
lines(c(xos,xoe),c(yos,yos),lwd=1,lty=1,col=1)
lines(c(xos,xos),c(yos,yoe),lwd=1,lty=1,col=1)

#---------------------------------------------------------axis labels
for (i in seq(yos,yoe-.5,1)){
  lines(c(xos-.12,xos),c(i,i),lwd=1,lty=1,col=1)
}
for (i in seq(yos,yoe-.5,1)){
  text(xos,i,format(round(.2*(i-7)+2.0,digits=2),nsmall=2),pos=2,cex=1,col=1)
}
for (i in xos+c(xfs:xfe)*1*xfm){
  lines(c(i,i),c(yos,yos-.1),lwd=1,lty=1,col=1)
}
for (i in c(1986,1988,1990,1992,1994,1996)){
  text(xos+(i-1986+2)*xfm,yos,i,pos=1,cex=1,col=1)
}
text(xos,yos+8.0,expression(paste(italic('Ovis aries')*
	': adult body weight')),pos=4,cex=1.3,col=1)
text(xos-1.9,yos+4,'Ln weight (kg)',
	srt=90,cex=1.3,col=1)
text(xos+7,yos-1.2,'Year',
	cex=1.3,col=1)
#============================================================ Load file
file1<-"C://R_aaROEVchapt08//8.1.17_Milner&al1999_Ovisaries_1wgt.csv"
M<-read.csv(file1,header=TRUE,sep=",",quote="\"",dec=".",fill=TRUE);M
#----------------------------Fill MF female matrix for 'all' and 'survivors'
MF<-matrix(nrow=6,ncol=12); rownames(MF)=c('Na','Ma','Va','Ns','Ms','Vs')
  colnames(MF)=as.character(seq(1985,1996))
for (i in 1:12){
  MF[1,i]=M[1,(i+3)]+M[4,(i+3)]
  if (M[4,(i+3)]>0){
    MF[2,i]=(M[1,(i+3)]*M[2,(i+3)]+M[4,(i+3)]*M[5,(i+3)])/MF[1,i]
  }else{
    MF[2,i]=M[2,(i+3)]
  }
  if (M[4,(i+3)]>1){
    MF[3,i]=(M[1,(i+3)]*M[3,(i+3)]+M[4,(i+3)]*M[6,(i+3)])/MF[1,i]
  }else{
    MF[3,i]=M[3,(i+3)]
  }
  MF[4,i]=M[1,(i+3)];MF[5,i]=M[2,(i+3)];MF[6,i]=M[3,(i+3)]
};MF=t(MF);MF
#----------------------------Fill MM male matrix for 'all' and 'survivors'
MM<-matrix(nrow=6,ncol=12); rownames(MM)=c('Na','Ma','Va','Ns','Ms','Vs')
  colnames(MM)=as.character(seq(1985,1996))
for (i in 1:12){
  MM[1,i]=M[7,(i+3)]+M[10,(i+3)]
  if (M[10,(i+3)]>0){
    MM[2,i]=(M[7,(i+3)]*M[8,(i+3)]+M[10,(i+3)]*M[11,(i+3)])/MM[1,i]
  }else{
    MM[2,i]=M[8,(i+3)]
  }
  if (M[10,(i+3)]>1){
    MM[3,i]=(M[7,(i+3)]*M[9,(i+3)]+M[10,(i+3)]*M[12,(i+3)])/MM[1,i]
  }else{
    MM[3,i]=M[9,(i+3)]
  }
  MM[4,i]=M[7,(i+3)];MM[5,i]=M[8,(i+3)];MM[6,i]=M[9,(i+3)]
};MM=t(MM);MM
#----------------------------------------------------------------------
ML<-matrix(nrow=12,ncol=7)
  colnames(ML)=c('Gen.','Nf','Ln.fx','Ln.fsd','Nm','Ln.mx','Ln.msd')
for (i in 1:length(ML[,1])){
  ML[i,1]=i						#interval
  ML[i,2]=MF[i,1]					#number of females
  ML[i,3]=log(MF[i,2])				#ln mean for females
  ML[i,4]=sqrt(MF[i,3])/MF[i,2]		#ln stdev for females
  ML[i,5]=MM[i,1]					#number of males	
  ML[i,6]=log(MM[i,2])				#ln mean for males
  ML[i,7]=sqrt(MM[i,3])/MM[i,2]		#ln stdev for males
}
ML[c(1:2,11:12),]
#------------------------------------------------Plot male ln weights
for (i in 1:12){
  lines(xos+.03+c(ML[i,1],ML[i,1])*xfm,
    5*c(ML[i,6]-ML[i,7],ML[i,6]+ML[i,7])-3,
    lty=1,lwd=1,col=gray(6/10))
}
for (i in 2:12){
  lines(xos+c(ML[i-1,1],ML[i,1])*xfm,5*c(ML[i-1,6],ML[i,6])-3,
	lty=1,lwd=1,col=1)
}
for (i in 1:12){
  points(xos+ML[i,1]*xfm,5*ML[i,6]-3,pch=19,cex=1.1,col=1)
}
#------------------------------------------------Plot female ln weights
for (i in 1:12){
  lines(xos-.03+c(ML[i,1],ML[i,1])*xfm,
    5*c(ML[i,3]-ML[i,4],ML[i,3]+ML[i,4])-3,
    lty=1,lwd=1,col=gray(6/10))
}
for (i in 2:12){
  lines(xos+c(ML[i-1,1],ML[i,1])*xfm,5*c(ML[i-1,3],ML[i,3])-3,
	lty=1,lwd=1,col=1)
}
for (i in 1:12){
  points(xos+ML[i,1]*xfm,5*ML[i,3]-3,pch=21,cex=1.2,bg='white',col=1)
}
#============================================= Rate calc: female weights
gentime=4								#generation time
rect(3.6,9.2,9,9.6,col='white',border=NA)
text(3.5,9.4,"Generation time: 4 years",pos=4,cex=.9,col=1)
idrsum=0
ML[1:3,]
n=nrow(ML);n
nn=.5*(n-1)*n;nn
idr=matrix(nrow=nn,ncol=9)
  colnames(idr)=c('int','diff.sd','rate.sd.g','log.i','log|d|','log|r|',
    'sbn','wgt','fgen')		#fgen is fraction of a generation
nc=0
for (k in 1:(n-1)){	#run length
  for (i in 1:(n-k)){  						#starting position
    nc=nc+1
    if(k<=gentime){idr[nc,1]=1} else {idr[nc,1]=k/gentime}	#interval
    meandiff=ML[(i+k),3]-ML[i,3]					#mean diff.
    poolsd=PoolSD(ML[i+k,2],ML[i,2],ML[i+k,4],ML[i,4])	#n1,n2,sd1,sd2
    idr[nc,2]=meandiff/poolsd						#diff.sd
    idr[nc,3]=abs(idr[nc,2])/idr[nc,1]				#rate.sd.gen
    idr[nc,4]=log10(idr[nc,1])					#log.i
    idr[nc,5]=log10(abs(idr[nc,2]))					#log.d
    idr[nc,6]=log10(idr[nc,3])					#log.r
    if(k==1){idr[nc,7]=1}else{idr[nc,7]=3}			#sbn
    idr[nc,8]=1/idr[nc,1]						#wgt
    idr[nc,9]=k/gentime							#fract of gen.
  }
}
idr[1:3,]
idrsum=idrsum+nrow(idr)
idrx=idr[!rowSums(!is.finite(idr)),]  #remove rows that have -Inf
idrxF=idrx[,1:8]
idrc=idrx[idrx[,9]>=1,]  #keep only rows where interval >= gentime
idrcF=idrc[,1:8]
nrow(idr);nrow(idrx);nrow(idrxF);nrow(idrcF)
#============================================= Rate calc: male weights
ML[1:3,]
n=nrow(ML);n
nn=.5*(n-1)*n;nn
idr=matrix(nrow=nn,ncol=9)
  colnames(idr)=c('int','diff.sd','rate.sd.g','log.i','log|d|','log|r|',
    'sbn','wgt','fgen')		#fgen is fraction of a generation
nc=0
for (k in 1:(n-1)){    #run length
  for (i in 1:(n-k)){  #starting position
    nc=nc+1
    if(k<=gentime){idr[nc,1]=1} else {idr[nc,1]=k/gentime}	#interval
    meandiff=ML[(i+k),6]-ML[i,6]					#mean diff.
    poolsd=PoolSD(ML[i+k,5],ML[i,5],ML[i+k,7],ML[i,7])	#n1,n2,sd1,sd2
    idr[nc,2]=meandiff/poolsd						#diff.sd
    idr[nc,3]=abs(idr[nc,2])/idr[nc,1]				#rate.sd
    idr[nc,4]=log10(idr[nc,1])					#log.i
    idr[nc,5]=log10(abs(idr[nc,2]))					#log.d
    idr[nc,6]=log10(idr[nc,3])					#log.r
    if(k==1){idr[nc,7]=1}else{idr[nc,7]=3}			#sbn
    idr[nc,8]=1/idr[nc,1]						#wgt
    idr[nc,9]=k/gentime							#fract of gen.
  }
}
idr[1:3,]
idrsum=idrsum+nrow(idr)
idrx=idr[!rowSums(!is.finite(idr)),]  #remove rows that have -Inf
idrxM=idrx[,1:8]
idrc=idrx[idrx[,9]>=1,]  #keep only rows where interval >= gentime
idrcM=idrc[,1:8]
nrow(idr);nrow(idrx);nrow(idrxM);nrow(idrcM)
#========================================================== Write 'all' file
idrxA=rbind(idrxF,idrxM)
writefile<-"C://R_aaROEVchapt08//8.1.17_Milner&al1999_Ovisaries_1wgt_out.csv"
write.csv(idrxA,file=writefile,na="NA",row.names=FALSE) 

idrsum;nrow(idrxA)
#========================================================== Plot LRI panel b
#assign 'bootn' as boot number
bootn<-1000;text(-1.7,-2.6,paste("Boot n =",bootn),pos=4,cex=1,col=4)
#assign 'mode' as "medians","all","mixed"
mode<-"all";text(2,-2.6,paste("Mode: ",mode),pos=4,cex=1,col=4)
#assign circle size for points (1.5 or 2)
psize<-1.5	#2	
#assign 'equation' position as "normal","lower","none" at end of each call
#send (1)idrx matrix, (2)mode(diff/rate), (3)panel placement coordinate x, 
#  (4)panel placement coordinate y, (5)bootn, (6)mode, (7)psize, (8)equation
bootresultd=TriPanelBC(idrcF,"r",1,4.5,bootn,mode,psize,"lower")
wrlb.d=bootresultd[[1]];bootmat.d=bootresultd[[2]]
proc.time()-ptm
#========================================================== Plot LRI panel c
#send idrx matrix, n, mode(diff/rate), panel placement coordinates, mode
bootresultr=TriPanelBC(idrcM,"r",13,4.5,bootn,mode,psize,"lower")
wrlb.r=bootresultr[[1]];bootmat.r=bootresultr[[2]]
proc.time()-ptm
#========================================================== Label B and C
rect(-.5,6.7,21.5,7.5,border=NA,col='white')
rect(.1,6,9,6.7,border=NA,col='white')
rect(12.1,6,21,6.7,border=NA,col='white')
text(0,6.5,"Female weight",pos=4,cex=1.3,col=1)
text(12,6.5,"Male weight",pos=4,cex=1.3,col=1)

#----------------------------------------------------------
text(7,-2.6,paste("Processing time: ",round(proc.time()[3]-ptm[3],digits=2),
  "seconds"),pos=4,cex=1,col=4)


