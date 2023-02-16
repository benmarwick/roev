##Philip D. Gingerich RATES OF EVOLUTION (Cambridge University Press 2019)
##8.1.39_Geiger&al2018_Musmusculus
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
xos=3.5;xoe=18;yos=9;yoe=17			#x,y original start and end
xfs=0;xfe=13;yfs=2.9;yfe=3.3			#x,y foreground start and end
xfm=(xoe-xos)/(xfe-xfs)				#x foreground multiplier
yfm=(yoe-yos)/(yfe-yfs)				#y foreground multiplier
lines(c(xos,xoe),c(yos,yos),lwd=1,lty=1,col=1)
lines(c(xos,xos),c(yos,yoe),lwd=1,lty=1,col=1)

#---------------------------------------------------------axis labels
for (i in seq(yos,yoe-.05,1)){
  lines(c(xos-.12,xos),c(i,i),lwd=1,lty=1,col=1)
}
for (i in seq(yos,yoe-.05,1)){
  text(xos,i,format(round(.1*(i-7)+1.2,digits=2),nsmall=2),pos=2,cex=1,col=1)
}
#-----------------------------------------------
for (i in xos+c(xfs:xfe)*1*xfm){
  lines(c(i,i),c(yos,yos-.1),lwd=1,lty=1,col=1)
}
for (i in c(1:13)){
  text(xos+(i)*xfm,yos,i,pos=1,cex=1,col=1)
}
#-----------------------------------------------
text(xos,yos+8.0,expression(paste(italic('Mus musculus')*
	': mouse pup body weight')),pos=4,cex=1.3,col=1)
text(xos-1.9,yos+4,'Ln weight (g)',
	srt=90,cex=1.3,col=1)
text(xos+7.5,yos-1.2,'Generation',
	cex=1.3,col=1)

#============================================================ Load file
file1<-"../Gingerich2019_AnalysisFiles_FieldStudies/8.1.39_Geiger&al2018_Musmusculus.csv"
G<-read.csv(file1,header=TRUE,sep=",",quote="\"",dec=".",fill=TRUE);G[c(1:3,nrow(G)),]
dayno=as.integer(difftime(G[,2],G[1,2],units="days"))+1;max(dayno);max(dayno)/365
G=cbind(G,as.integer(dayno));G[1:3,]
gendays=264
cutpt=seq(0,as.integer(max(dayno)),gendays);cutpt
LG=matrix(nrow=0,ncol=7);colnames(LG)=c('gen','Nw','ln Mw','ln Sw','Nh','ln Mh','ln Sh')
for (i in 2:14){    #14){
  nc=0;temp<-matrix(nrow=0,ncol=2);colnames(temp)=c('wgt','hlen')
  for(j in 1:nrow(G)){
    if (G[j,10]>cutpt[i-1]&&G[j,10]<=cutpt[i]){
      temp=rbind(temp, c(log(G[j,8]),log(G[j,9])) )
    }
  }
  LG=rbind(LG, c(i-1,length(temp[,1]),mean(temp[,1]),sd(temp[,1]),
    length(temp[,2]),mean(temp[,2]),sd(temp[,2])) )
}
LG;sum(LG[,2])
#------------------------------------------------Plot body weight time series
for (i in 1:nrow(LG)){
  lines(xos-.04+c(LG[i,1],LG[i,1])*xfm,
    10*c(LG[i,3]-LG[i,4],LG[i,3]+LG[i,4])-5,
    lty=1,lwd=1,col=gray(6/10))
}
for (i in 2:nrow(LG)){
  lines(xos+c(LG[i-1,1],LG[i,1])*xfm,10*c(LG[i-1,3],LG[i,3])-5,
	lty=1,lwd=1,col=1)
}
for (i in 1:nrow(LG)){
  points(xos+LG[i,1]*xfm,10*LG[i,3]-5,pch=19,col=1)
}

#============================================= Rate calculation
n=nrow(LG);n
nn=.5*(n-1)*n;nn
#---------------------------------------------------------Rates for pup weights
idr=matrix(nrow=2*nn,ncol=9)
  colnames(idr)=c('int','diff.sd','rate.sd.g','log.i','log|d|','log|r|',
    'sbn','wgt','sum')		#fgen is fraction of a generation
nc=0
for (k in 1:(n-1)){		#run length
  for (i in 1:(n-k)){  		#starting position
    nc=nc+1
    idr[nc,1]=k
    meandiff=LG[(i+k),3]-LG[i,3]					#mean diff.
    poolsd=PoolSD(LG[i+k,2],LG[i,2],LG[i+k,4],LG[i,4])	#n1,n2,sd1,sd2
    idr[nc,2]=meandiff/poolsd						#diff.sd
    idr[nc,3]=abs(idr[nc,2])/idr[nc,1]				#rate.sd.gen
    idr[nc,4]=log10(idr[nc,1])					#log.i
    idr[nc,5]=log10(abs(idr[nc,2]))					#log.d
    idr[nc,6]=log10(idr[nc,3])					#log.r
    if(idr[nc,1]==1){idr[nc,7]=1}else{idr[nc,7]=3}		#sbn
    idr[nc,8]=1/idr[nc,1]						#wgt
    idr[nc,9]=idr[nc,4]+idr[nc,6]					#sum
  }
}
idr[1:3,]
idrxW=idr[!rowSums(!is.finite(idr)),]  #remove rows that have -Inf
nrow(idr);nrow(idrxW)
#-----------------------------------------------------Rates for pup head lengths
idr=matrix(nrow=2*nn,ncol=9)
  colnames(idr)=c('int','diff.sd','rate.sd.g','log.i','log|d|','log|r|',
    'sbn','wgt','sum')		#fgen is fraction of a generation
nc=0
for (k in 1:(n-1)){		#run length
  for (i in 1:(n-k)){  		#starting position
    nc=nc+1
    idr[nc,1]=k
    meandiff=LG[(i+k),6]-LG[i,6]					#mean diff.
    poolsd=PoolSD(LG[i+k,5],LG[i,5],LG[i+k,7],LG[i,7])	#n1,n2,sd1,sd2
    idr[nc,2]=meandiff/poolsd						#diff.sd
    idr[nc,3]=abs(idr[nc,2])/idr[nc,1]				#rate.sd.gen
    idr[nc,4]=log10(idr[nc,1])					#log.i
    idr[nc,5]=log10(abs(idr[nc,2]))					#log.d
    idr[nc,6]=log10(idr[nc,3])					#log.r
    if(idr[nc,1]==1){idr[nc,7]=1}else{idr[nc,7]=3}		#sbn
    idr[nc,8]=1/idr[nc,1]						#wgt
    idr[nc,9]=idr[nc,4]+idr[nc,6]					#sum
  }
}
idr[1:3,]
idrxH=idr[!rowSums(!is.finite(idr)),]  #remove rows that have -Inf
nrow(idr);nrow(idrxH)
idrxA=rbind(idrxW,idrxH)
writefile<-"C://R_aaROEVchapt08//8.1.39_Geiger&al2018_Musmusculus_out.csv"
write.csv(idrxA[,1:8],file=writefile,na="NA",row.names=FALSE)

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
bootresultd=TriPanelBC(idrxW,"r",1,4.5,bootn,mode,psize,"lower")
wrlb.d=bootresultd[[1]];bootmat.d=bootresultd[[2]]
proc.time()-ptm
#========================================================== Plot LRI panel c
#send idrx matrix, n, mode(diff/rate), panel placement coordinates, mode
bootresultr=TriPanelBC(idrxH,"r",13,4.5,bootn,mode,psize,"lower")
wrlb.r=bootresultr[[1]];bootmat.r=bootresultr[[2]]
proc.time()-ptm
#========================================================== Label B and C
rect(-.5,6.7,21.5,7.5,border=NA,col='white')
rect(.1,6,9,6.7,border=NA,col='white')
rect(12.1,6,21,6.7,border=NA,col='white')
text(0,6.5,"Mouse pup body weight",pos=4,cex=1.3,col=1)
text(12,6.5,"Mouse pup head length",pos=4,cex=1.3,col=1)

#----------------------------------------------------------
text(7,-2.6,paste("Processing time: ",round(proc.time()[3]-ptm[3],digits=2),
  "seconds"),pos=4,cex=1,col=4)
warnings()


