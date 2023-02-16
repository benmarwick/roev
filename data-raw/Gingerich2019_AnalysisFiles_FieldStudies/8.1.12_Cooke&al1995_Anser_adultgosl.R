##Philip D. Gingerich RATES OF EVOLUTION (Cambridge University Press 2019)
##8.1.12_Cooke&al1995_Anser_adultgosl
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
	xaxp=c(xr[1],xr[2],xr[1]+xr[2]),
	yaxp=c(yr[1],yr[2],abs(yr[1]-yr[2])),		#y-axis extreme ticks and number of ticks
	cex.axis=.8,
	cex.lab=1.1)			#label size
dev.size('in')
#======================================================= Plot panel a
xos=3.5;xoe=17;yos=9;yoe=17			#x,y original start and end
xfs=0;xfe=18;yfs=5;yfe=5.4			#x,y foreground start and end
xfm=(xoe-xos)/(xfe-xfs)				#x foreground multiplier
yfm=(yoe-yos)/(yfe-yfs)				#y foreground multiplier
lines(c(xos,xoe),c(yos,yos),lwd=1,lty=1,col=1)
lines(c(xos,xos),c(yos,yoe),lwd=1,lty=1,col=1)

#---------------------------------------------------------axis labels
for (i in seq(yos,yoe-.5,1)){
  lines(c(xos-.12,xos),c(i,i),lwd=1,lty=1,col=1)
}
for (i in seq(yos,yoe-.5,2)){
  text(xos,i+1,format(round(.05*(i-7)+7.3,digits=2),nsmall=2),pos=2,cex=1,col=1)
}
for (i in xos+c(xfs:xfe)*1*xfm){
  lines(c(i,i),c(yos,yos-.1),lwd=1,lty=1,col=1)
}
for (i in seq(1970,1986,2)){
  text(xos+(i-1970+2)*xfm,yos,i,pos=1,cex=1,col=1)
}
text(xos,yos+8.0,expression(paste(italic('Anser caerulescens')*
	': female body weight')),pos=4,cex=1.3,col=1)
text(xos-1.9,yos+4,'Ln weight (g)',
	srt=90,cex=1.3,col=1)
text(xos+7,yos-1.2,'Year',
	cex=1.3,col=1)
#============================================================ Load file
file1<-"../Gingerich2019_AnalysisFiles_FieldStudies/8.1.12_Cooke&al1995_Anser_1adultwgt.csv"
file2<-"../Gingerich2019_AnalysisFiles_FieldStudies/8.1.12_Cooke&al1995_Anser_2goslall.csv"
A<-read.csv(file1,header=TRUE,sep=",",quote="\"",dec=".",fill=TRUE);A
G<-read.csv(file2,header=TRUE,sep=",",quote="\"",dec=".",fill=TRUE);G
#-------------------------------------------------Adult female weight
LA<-matrix(nrow=nrow(A),ncol=7)
  colnames(LA)=c('Gen','Yr','N','Mean','Stdev','Ln.M','Ln.SD')
for (i in 1:nrow(A)){
  LA[i,1]=i						#Gen
  LA[i,2]=A[i,1]					#Year
  LA[i,3]=A[i,2]					#N
  LA[i,4]=A[i,3]					#Mean
  LA[i,5]=sqrt(A[i,2])*A[i,4]			#Stdev
  LA[i,6]=log(LA[i,4])				#Ln mean
  LA[i,7]=LA[i,5]/LA[i,4]			#Ln stdev
}
LA[c(1:2,(nrow(LA)-1):nrow(LA)),]
#-------------------------------------------------Gosling weight
LW<-matrix(nrow=nrow(G),ncol=7)
  colnames(LW)=c('Gen','Yr','N','Mean','Stdev','Ln.M','Ln.SD')
for (i in 1:nrow(G)){
  LW[i,1]=i						#Gen
  LW[i,2]=G[i,1]					#Year
  LW[i,3]=G[i,2]					#N
  LW[i,4]=G[i,3]					#Mean
  LW[i,5]=sqrt(G[i,2])*G[i,4]			#Stdev
  LW[i,6]=log(LW[i,4])				#Ln mean
  LW[i,7]=LW[i,5]/LW[i,4]			#Ln stdev
}
LW[c(1:2,(nrow(LW)-1):nrow(LW)),]
#-------------------------------------------------Gosling tarsus length
LT<-matrix(nrow=nrow(G),ncol=7)
  colnames(LT)=c('Gen','Yr','N','Mean','Stdev','Ln.M','Ln.SD')
for (i in 1:nrow(G)){
  LT[i,1]=i						#Gen
  LT[i,2]=G[i,1]					#Year
  LT[i,3]=G[i,5]					#N
  LT[i,4]=G[i,6]					#Mean
  LT[i,5]=G[i,7]					#Stdev
  LT[i,6]=log(LT[i,4])				#Ln mean
  LT[i,7]=LT[i,5]/LT[i,4]			#Ln stdev
}
LT[c(1:2,(nrow(LT)-1):nrow(LT)),]
#-------------------------------------------------Gosling beak length
LB<-matrix(nrow=nrow(G),ncol=7)
  colnames(LB)=c('Gen','Yr','N','Mean','Stdev','Ln.M','Ln.SD')
for (i in 1:nrow(G)){
  LB[i,1]=i						#Gen
  LB[i,2]=G[i,1]					#Year
  LB[i,3]=G[i,8]					#N
  LB[i,4]=G[i,9]					#Mean
  LB[i,5]=G[i,10]					#Stdev
  LB[i,6]=log(LB[i,4])				#Ln mean
  LB[i,7]=LB[i,5]/LB[i,4]			#Ln stdev
}
LB[c(1:2,(nrow(LB)-1):nrow(LB)),]
#------------------------------------------------Plot adult female weights
for (i in 1:nrow(LA)){
  lines(xos+c(LA[i,1],LA[i,1])*xfm,
    20*c(LA[i,6]-LA[i,7],LA[i,6]+LA[i,7])-138,
    lty=1,lwd=1,col=gray(6/10))
}
for (i in 2:nrow(LA)){
  lines(xos+c(LA[i-1,1],LA[i,1])*xfm,20*c(LA[i-1,6],LA[i,6])-138,
	lty=1,lwd=1,col=1)
}
for (i in 1:nrow(LA)){
  points(xos+(LA[i,1])*xfm,20*LA[i,6]-138,pch=19,cex=1.1,col=1)
}
#========================================== Rate calc: adult female weights
gentime=5				#generation time fide Niel and Lebreton 2005
rect(3.6,9.2,9,9.6,col='white',border=NA)
text(3.5,9.4,"Generation time: 5 years",pos=4,cex=.9,col=1)
LA[1:3,]
n=nrow(LA);n
nn=.5*(n-1)*n;nn
idr=matrix(nrow=nn,ncol=9)
  colnames(idr)=c('int','diff.sd','rate.sd.g','log.i','log|d|','log|r|',
    'sbn','wgt','fgen')		#fgen is fraction of a generation
nc=0
for (k in 1:(n-1)){	#run length
  for (i in 1:(n-k)){  						#starting position
    nc=nc+1
    if(k<=gentime){idr[nc,1]=1} else {idr[nc,1]=k/gentime}	#interval
    meandiff=LA[(i+k),6]-LA[i,6]				                  	#mean diff.
    poolsd=PoolSD(LA[i+k,3],LA[i,3],LA[i+k,7],LA[i,7])     	#n1,n2,sd1,sd2
    idr[nc,2]=meandiff/poolsd					                    	#diff.sd
    idr[nc,3]=abs(idr[nc,2])/idr[nc,1]	              			#rate.sd.gen
    idr[nc,4]=log10(idr[nc,1])			                    		#log.i
    idr[nc,5]=log10(abs(idr[nc,2]))	                				#log.d
    idr[nc,6]=log10(idr[nc,3])			                    		#log.r
    if(k==1){idr[nc,7]=1}else{idr[nc,7]=3}	            		#sbn
    idr[nc,8]=1/idr[nc,1]				                        		#wgt
    idr[nc,9]=k/gentime				                        			#fract of gen.
  }
}
idr[1:3,]
idrx=idr[!rowSums(!is.finite(idr)),]  #remove rows that have -Inf
idrxA=idrx[,1:8]
idrc=idrx[idrx[,9]>=1,]  #keep only rows where interval >= gentime
idrcA=idrc[,1:8]
nrow(idr);nrow(idrx);nrow(idrxA);nrow(idrcA)
writefile<-"C://R_aaROEVchapt08//8.1.12_Cooke&al1995_Anser_1adultwgt_out.csv"
write.csv(idrxA,file=writefile,na="NA",row.names=FALSE)
#============================================= Rate calc: gosling weights
LW[1:3,]
idrsum=0
n=nrow(LW);n
nn=.5*(n-1)*n;nn
idr=matrix(nrow=nn,ncol=9)
  colnames(idr)=c('int','diff.sd','rate.sd.g','log.i','log|d|','log|r|',
    'sbn','wgt','fgen')		#fgen is fraction of a generation
nc=0
for (k in 1:(n-1)){	#run length
  for (i in 1:(n-k)){  						#starting position
    nc=nc+1
    if(k<=gentime){idr[nc,1]=1} else {idr[nc,1]=k/gentime}	#interval
    meandiff=LW[(i+k),6]-LW[i,6]					#mean diff.
    poolsd=PoolSD(LW[i+k,3],LW[i,3],LW[i+k,7],LW[i,7])	#n1,n2,sd1,sd2
    idr[nc,2]=meandiff/poolsd						#diff.sd
    idr[nc,3]=abs(idr[nc,2])/idr[nc,1]				#rate.sd.gen
    idr[nc,4]=log10(idr[nc,1])					#log.i
    idr[nc,5]=log10(abs(idr[nc,2]))					#log.d
    idr[nc,6]=log10(idr[nc,3])					#log.r
    if(k==1){idr[nc,7]=1}else{idr[nc,7]=3}			#sbn
    idr[nc,8]=1/idr[nc,1]						#wgt
    idr[nc,9]=k/gentime						#fract of gen.
  }
}
idr[1:3,]
idrsum=idrsum+nrow(idr)
idrx=idr[!rowSums(!is.finite(idr)),]  #remove rows that have -Inf
idrxW=idrx[,1:8]
idrc=idrx[idrx[,9]>=1,]  #keep only rows where interval >= gentime
idrcW=idrc[,1:8]
nrow(idr);nrow(idrx);nrow(idrxW);nrow(idrcW)
#========================================= Rate calc: gosling tarsus length
LT=LT[!is.na(LT[,4]),];LT		#remove rows with no mean
n=nrow(LT);n
nn=.5*(n-1)*n;nn
idr=matrix(nrow=nn,ncol=9)
  colnames(idr)=c('int','diff.sd','rate.sd.g','log.i','log|d|','log|r|',
    'sbn','wgt','fgen')		#fgen is fraction of a generation
nc=0
for (k in 1:(n-1)){	#run length
  for (i in 1:(n-k)){  						#starting position
    nc=nc+1
    if(k<=gentime){idr[nc,1]=1} else {idr[nc,1]=k/gentime}	#interval
    meandiff=LT[(i+k),6]-LT[i,6]					#mean diff.
    poolsd=PoolSD(LT[i+k,3],LT[i,3],LT[i+k,7],LT[i,7])	#n1,n2,sd1,sd2
    idr[nc,2]=meandiff/poolsd						#diff.sd
    idr[nc,3]=abs(idr[nc,2])/idr[nc,1]				#rate.sd.gen
    idr[nc,4]=log10(idr[nc,1])					#log.i
    idr[nc,5]=log10(abs(idr[nc,2]))					#log.d
    idr[nc,6]=log10(idr[nc,3])					#log.r
    if(k==1){idr[nc,7]=1}else{idr[nc,7]=3}			#sbn
    idr[nc,8]=1/idr[nc,1]						#wgt
    idr[nc,9]=k/gentime						#fract of gen.
  }
}
idr[1:3,]
idrsum=idrsum+nrow(idr)
idrx=idr[!rowSums(!is.finite(idr)),]  #remove rows that have -Inf
idrxT=idrx[,1:8]
idrc=idrx[idrx[,9]>=1,]  #keep only rows where interval >= gentime
idrcT=idrc[,1:8]
nrow(idr);nrow(idrx);nrow(idrxT);nrow(idrcT)
#=========================================== Rate calc: gosling beak length
LB[1:3,]
n=nrow(LB);n
nn=.5*(n-1)*n;nn
idr=matrix(nrow=nn,ncol=9)
  colnames(idr)=c('int','diff.sd','rate.sd.g','log.i','log|d|','log|r|',
    'sbn','wgt','fgen')		#fgen is fraction of a generation
nc=0
for (k in 1:(n-1)){	#run length
  for (i in 1:(n-k)){  						#starting position
    nc=nc+1
    #idr[nc,1]=k/gentime						#interval
    if(k<=gentime){idr[nc,1]=1} else {idr[nc,1]=k/gentime}	#interval
    meandiff=LB[(i+k),6]-LB[i,6]				#mean diff.
    poolsd=PoolSD(LB[i+k,3],LB[i,3],LB[i+k,7],LB[i,7])#n1,n2,sd1,sd2
    idr[nc,2]=meandiff/poolsd					#diff.sd
    idr[nc,3]=abs(idr[nc,2])/idr[nc,1]			#rate.sd.gen
    idr[nc,4]=log10(idr[nc,1])				#log.i
    idr[nc,5]=log10(abs(idr[nc,2]))				#log.d
    idr[nc,6]=log10(idr[nc,3])				#log.r
    if(k==1){idr[nc,7]=1}else{idr[nc,7]=3}			#sbn
    idr[nc,8]=1/idr[nc,1]						#wgt
    idr[nc,9]=k/gentime						#fract of gen.
  }
}
idr[1:3,]
idrsum=idrsum+nrow(idr)
idrx=idr[!rowSums(!is.finite(idr)),]  #remove rows that have -Inf
idrxB=idrx[,1:8]
idrc=idrx[idrx[,9]>=1,]  #keep only rows where interval >= gentime
idrcB=idrc[,1:8]
nrow(idr);nrow(idrx);nrow(idrxB);nrow(idrcB)
#========================================================== Write gosling file
idrxG=rbind(idrxW,idrxT,idrxB)
  colnames(idrxG)=c('int','diff.sd','rate.sd.g','log.i','log|d|','log|r|','wgt')
idrsum
nrow(idrxG)
writefile<-"C://R_aaROEVchapt08//8.1.12_Cooke&al1995_Anser_2goslall_out.csv"
write.csv(idrxG,file=writefile,na="NA",row.names=FALSE)

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
bootresultd=TriPanelBC(idrcA,"r",1,4.5,bootn,mode,psize,"lower")
wrlb.d=bootresultd[[1]];bootmat.d=bootresultd[[2]]
proc.time()-ptm
#========================================================== Plot LRI panel c
#send idrx matrix, n, mode(diff/rate), panel placement coordinates, mode
idrcG=rbind(idrcW,idrcT,idrcB)
bootresultr=TriPanelBC(idrcG,"r",13,4.5,bootn,mode,psize,"lower")
wrlb.r=bootresultr[[1]];bootmat.r=bootresultr[[2]]
proc.time()-ptm
#========================================================== Label B and C
rect(-.5,6.7,21.5,7.5,border=NA,col='white')
rect(.1,6,9,6.7,border=NA,col='white')
rect(12.1,6,21,6.7,border=NA,col='white')
text(0,6.5,"Adult female weight",pos=4,cex=1.3,col=1)
text(12,6.5,"Gosling weight, etc.",pos=4,cex=1.3,col=1)

#----------------------------------------------------------
text(7,-2.6,paste("Processing time: ",round(proc.time()[3]-ptm[3],digits=2),
  "seconds"),pos=4,cex=1,col=4)
warnings()

