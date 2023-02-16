##Philip D. Gingerich RATES OF EVOLUTION (Cambridge University Press 2019) 
##8.1.16_Koontz&al2001_Dipodomys
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
xfs=0;xfe=8;yfs=5;yfe=5.4			#x,y foreground start and end
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
for (i in seq(1989,1996,1)){
  text(xos+(i-1989+1)*xfm,yos,i,pos=1,cex=1,col=1)
}
text(xos,yos+8.0,expression(paste(italic('Dipodomys merriami')*
	': body weight')),pos=4,cex=1.3,col=1)
text(xos-1.9,yos+4,'Ln weight (g)',
	srt=90,cex=1.3,col=1)
text(xos+7,yos-1.2,'Year',
	cex=1.3,col=1)
#============================================================ Load file
file1<-"C://R_aaROEVchapt08//8.1.16_Koontz&al2001_Dipodomys.csv"
K<-read.csv(file1,header=TRUE,sep=",",quote="\"",dec=".",fill=TRUE)
K[1:3,]
#-------------------------------------------------Body weight locality A
KA<-matrix(nrow=nrow(K[1:8,]),ncol=12)
  colnames(KA)=c('Gen','Yr','N','Mean','Stdev','N','Mean','Stdev',
    'Ln.Mf','Ln.SDf','Ln.Mm','Ln.SDm')
for (i in 1:nrow(K[1:8,])){
  KA[i,1]=i						#Gen
  KA[i,2]=K[i,1]					#Year
  KA[i,3]=K[i,4]					#N
  KA[i,4]=K[i,5]					#Mean
  KA[i,5]=K[i,6]					#Stdev
  KA[i,6]=K[i,7]					#N
  KA[i,7]=K[i,8]					#Mean
  KA[i,8]=K[i,9]					#Stdev
  KA[i,9]=log(KA[i,4])				#Ln mean
  KA[i,10]=KA[i,5]/KA[i,4]			#Ln stdev
  KA[i,11]=log(KA[i,7])				#Ln mean
  KA[i,12]=KA[i,8]/KA[i,7]			#Ln stdev
}
KA
#-------------------------------------------------Body weight locality B
KB<-matrix(nrow=nrow(K[1:8,]),ncol=12)
  colnames(KB)=c('Gen','Yr','N','Mean','Stdev','N','Mean','Stdev',
    'Ln.Mf','Ln.SDf','Ln.Mm','Ln.SDm')
for (i in 1:nrow(K[1:8,])){
  KB[i,1]=i						#Gen
  KB[i,2]=K[i,1]					#Year
  KB[i,3]=K[i,4]					#N
  KB[i,4]=K[i,5]					#Mean
  KB[i,5]=K[i,6]					#Stdev
  KB[i,6]=K[i,7]					#N
  KB[i,7]=K[i,8]					#Mean
  KB[i,8]=K[i,9]					#Stdev
  KB[i,9]=log(KB[i,4])				#Ln mean
  KB[i,10]=KB[i,5]/KB[i,4]			#Ln stdev
  KB[i,11]=log(KB[i,7])				#Ln mean
  KB[i,12]=KB[i,8]/KB[i,7]			#Ln stdev
}
KB
#-------------------------------------------------Body weight locality C
KC<-matrix(nrow=nrow(K[1:8,]),ncol=12)
  colnames(KC)=c('Gen','Yr','N','Mean','Stdev','N','Mean','Stdev',
    'Ln.Mf','Ln.SDf','Ln.Mm','Ln.SDm')
for (i in 1:nrow(K[1:8,])){
  KC[i,1]=i						#Gen
  KC[i,2]=K[i,1]					#Year
  KC[i,3]=K[i,4]					#N
  KC[i,4]=K[i,5]					#Mean
  KC[i,5]=K[i,6]					#Stdev
  KC[i,6]=K[i,7]					#N
  KC[i,7]=K[i,8]					#Mean
  KC[i,8]=K[i,9]					#Stdev
  KC[i,9]=log(KC[i,4])				#Ln mean
  KC[i,10]=KC[i,5]/KC[i,4]			#Ln stdev
  KC[i,11]=log(KC[i,7])				#Ln mean
  KC[i,12]=KC[i,8]/KC[i,7]			#Ln stdev
}
KC
#------------------------------------------------Plot KA weights
for (i in seq(yos,yoe-.5,1)){
  lines(xos+c(0,xfe)*xfm,c(i,i),lty=2,lwd=1,col=gray(6/10))
}
for (i in 1:nrow(KA)){
  lines(xos+c(KA[i,1],KA[i,1])*xfm,
    20*c(KA[i,11]-KA[i,12],KA[i,11]+KA[i,12])-62,
    lty=1,lwd=1,col=gray(6/10))
  lines(xos+c(KA[i,1],KA[i,1])*xfm,
    20*c(KA[i,9]-KA[i,10],KA[i,9]+KA[i,10])-62,
    lty=1,lwd=1,col=gray(6/10))
}
for (i in 2:nrow(KA)){
  lines(xos+c(KA[i-1,1],KA[i,1])*xfm,20*c(KA[i-1,11],KA[i,11])-62,
	lty=1,lwd=1,col=1)
  lines(xos+c(KA[i-1,1],KA[i,1])*xfm,20*c(KA[i-1,9],KA[i,9])-62,
	lty=1,lwd=1,col=1)
}
for (i in 1:nrow(KA)){
  points(xos+(KA[i,1])*xfm,20*KA[i,11]-62,pch=19,cex=1.1,col=1)
  points(xos+(KA[i,1])*xfm,20*KA[i,9]-62,pch=21,bg='white',cex=1.1,col=1)
}
#=========================================================== Rate calc: KAfm
gentime=1				#estimated
rect(3.6,9.2,9,9.6,col='white',border=NA)
text(3.5,9.4,"Generation time: 1 year",pos=4,cex=.9,col=1)
KA[1:3,]
n=nrow(KA);n
nn=.5*(n-1)*n;nn
idr=matrix(nrow=nn,ncol=7)
  colnames(idr)=c('int','diff.sd','rate.sd.g','log.i','log|d|','log|r|','wgt')
nc=0
stdev<-numeric(length=n-1)
for (k in gentime:(n-1)){	#run length
  for (i in 1:(n-k)){  						#starting position
    nc=nc+1
    idr[nc,1]=k/gentime						#intercept
    meandiff=KA[(i+k),9]-KA[i,9]				#mean diff.
    poolsd=PoolSD(KA[i+k,3],KA[i,3],KA[i+k,10],KA[i,10])#n1,n2,sd1,sd2
    idr[nc,2]=meandiff/poolsd					#diff.sd
    idr[nc,3]=abs(idr[nc,2])/idr[nc,1]			#rate.sd.gen
    idr[nc,4]=log10(idr[nc,1])				#log.i
    idr[nc,5]=log10(abs(idr[nc,2]))				#log.d
    idr[nc,6]=log10(idr[nc,3])				#log.r
    idr[nc,7]=1/idr[nc,1]					#wgt
    stdev[i]=poolsd
  }
}
vartot=sum(stdev^2,na.rm=TRUE);varcnt=length(table(stdev,useNA='no'))
varmean=vartot/varcnt; avestdev=sqrt(varmean);avestdev
idr[1:3,]
idrxKAf=idr[!rowSums(!is.finite(idr)),]  #remove rows that have -Inf
n=nrow(KA);n	#---------------------------------------------------------
nn=.5*(n-1)*n;nn
idr=matrix(nrow=nn,ncol=7)
  colnames(idr)=c('int','diff.sd','rate.sd.g','log.i','log|d|','log|r|','wgt')
nc=0
stdev<-numeric(length=n-1)
for (k in gentime:(n-1)){	#run length
  for (i in 1:(n-k)){  						#starting position
    nc=nc+1
    idr[nc,1]=k/gentime						#intercept
    meandiff=KA[(i+k),11]-KA[i,11]				#mean diff.
    poolsd=PoolSD(KA[i+k,3],KA[i,3],KA[i+k,12],KA[i,12])#n1,n2,sd1,sd2
    idr[nc,2]=meandiff/poolsd					#diff.sd
    idr[nc,3]=abs(idr[nc,2])/idr[nc,1]			#rate.sd.gen
    idr[nc,4]=log10(idr[nc,1])				#log.i
    idr[nc,5]=log10(abs(idr[nc,2]))				#log.d
    idr[nc,6]=log10(idr[nc,3])				#log.r
    idr[nc,7]=1/idr[nc,1]					#wgt
    stdev[i]=poolsd
  }
}
vartot=sum(stdev^2,na.rm=TRUE);varcnt=length(table(stdev,useNA='no'))
varmean=vartot/varcnt; avestdev=sqrt(varmean);avestdev
idr[1:3,]
idrxKAm=idr[!rowSums(!is.finite(idr)),]  #remove rows that have -Inf
#=========================================================== Rate calc: KBfm
KB[1:3,]
n=nrow(KB);n
nn=.5*(n-1)*n;nn
idr=matrix(nrow=nn,ncol=7)
  colnames(idr)=c('int','diff.sd','rate.sd.g','log.i','log|d|','log|r|','wgt')
nc=0
stdev<-numeric(length=n-1)
for (k in gentime:(n-1)){	#run length
  for (i in 1:(n-k)){  						#starting position
    nc=nc+1
    idr[nc,1]=k/gentime						#intercept
    meandiff=KB[(i+k),9]-KB[i,9]				#mean diff.
    poolsd=PoolSD(KB[i+k,3],KB[i,3],KB[i+k,10],KB[i,10])#n1,n2,sd1,sd2
    idr[nc,2]=meandiff/poolsd					#diff.sd
    idr[nc,3]=abs(idr[nc,2])/idr[nc,1]			#rate.sd.gen
    idr[nc,4]=log10(idr[nc,1])				#log.i
    idr[nc,5]=log10(abs(idr[nc,2]))				#log.d
    idr[nc,6]=log10(idr[nc,3])				#log.r
    idr[nc,7]=1/idr[nc,1]					#wgt
    stdev[i]=poolsd
  }
}
vartot=sum(stdev^2,na.rm=TRUE);varcnt=length(table(stdev,useNA='no'))
varmean=vartot/varcnt; avestdev=sqrt(varmean);avestdev
idr[1:3,]
idrxKBf=idr[!rowSums(!is.finite(idr)),]  #remove rows that have -Inf
n=nrow(KB);n	#---------------------------------------------------------
nn=.5*(n-1)*n;nn
idr=matrix(nrow=nn,ncol=7)
  colnames(idr)=c('int','diff.sd','rate.sd.g','log.i','log|d|','log|r|','wgt')
nc=0
stdev<-numeric(length=n-1)
for (k in gentime:(n-1)){	#run length
  for (i in 1:(n-k)){  						#starting position
    nc=nc+1
    idr[nc,1]=k/gentime						#intercept
    meandiff=KB[(i+k),11]-KB[i,11]				#mean diff.
    poolsd=PoolSD(KB[i+k,3],KB[i,3],KB[i+k,12],KB[i,12])#n1,n2,sd1,sd2
    idr[nc,2]=meandiff/poolsd					#diff.sd
    idr[nc,3]=abs(idr[nc,2])/idr[nc,1]			#rate.sd.gen
    idr[nc,4]=log10(idr[nc,1])				#log.i
    idr[nc,5]=log10(abs(idr[nc,2]))				#log.d
    idr[nc,6]=log10(idr[nc,3])				#log.r
    idr[nc,7]=1/idr[nc,1]					#wgt
    stdev[i]=poolsd
  }
}
vartot=sum(stdev^2,na.rm=TRUE);varcnt=length(table(stdev,useNA='no'))
varmean=vartot/varcnt; avestdev=sqrt(varmean);avestdev
idr[1:3,]
idrxKBm=idr[!rowSums(!is.finite(idr)),]  #remove rows that have -Inf
#=========================================================== Rate calc: KCfm
gentime=1				#estimated
rect(3.6,9.2,9,9.6,col='white',border=NA)
text(3.5,9.4,"Generation time: 1 year",pos=4,cex=.9,col=1)
KC[1:3,]
n=nrow(KC);n
nn=.5*(n-1)*n;nn
idr=matrix(nrow=nn,ncol=7)
  colnames(idr)=c('int','diff.sd','rate.sd.g','log.i','log|d|','log|r|','wgt')
nc=0
stdev<-numeric(length=n-1)
for (k in gentime:(n-1)){	#run length
  for (i in 1:(n-k)){  						#starting position
    nc=nc+1
    idr[nc,1]=k/gentime						#intercept
    meandiff=KC[(i+k),9]-KC[i,9]				#mean diff.
    poolsd=PoolSD(KC[i+k,3],KC[i,3],KC[i+k,10],KC[i,10])#n1,n2,sd1,sd2
    idr[nc,2]=meandiff/poolsd					#diff.sd
    idr[nc,3]=abs(idr[nc,2])/idr[nc,1]			#rate.sd.gen
    idr[nc,4]=log10(idr[nc,1])				#log.i
    idr[nc,5]=log10(abs(idr[nc,2]))				#log.d
    idr[nc,6]=log10(idr[nc,3])				#log.r
    idr[nc,7]=1/idr[nc,1]					#wgt
    stdev[i]=poolsd
  }
}
vartot=sum(stdev^2,na.rm=TRUE);varcnt=length(table(stdev,useNA='no'))
varmean=vartot/varcnt; avestdev=sqrt(varmean);avestdev
idr[1:3,]
idrxKCf=idr[!rowSums(!is.finite(idr)),]  #remove rows that have -Inf
n=nrow(KC);n	#---------------------------------------------------------
nn=.5*(n-1)*n;nn
idr=matrix(nrow=nn,ncol=7)
  colnames(idr)=c('int','diff.sd','rate.sd.g','log.i','log|d|','log|r|','wgt')
nc=0
stdev<-numeric(length=n-1)
for (k in gentime:(n-1)){	#run length
  for (i in 1:(n-k)){  						#starting position
    nc=nc+1
    idr[nc,1]=k/gentime						#intercept
    meandiff=KC[(i+k),11]-KC[i,11]				#mean diff.
    poolsd=PoolSD(KC[i+k,3],KC[i,3],KC[i+k,12],KC[i,12])#n1,n2,sd1,sd2
    idr[nc,2]=meandiff/poolsd					#diff.sd
    idr[nc,3]=abs(idr[nc,2])/idr[nc,1]			#rate.sd.gen
    idr[nc,4]=log10(idr[nc,1])				#log.i
    idr[nc,5]=log10(abs(idr[nc,2]))				#log.d
    idr[nc,6]=log10(idr[nc,3])				#log.r
    idr[nc,7]=1/idr[nc,1]					#wgt
    stdev[i]=poolsd
  }
}
vartot=sum(stdev^2,na.rm=TRUE);varcnt=length(table(stdev,useNA='no'))
varmean=vartot/varcnt; avestdev=sqrt(varmean);avestdev
idr[1:3,]
idrxKCm=idr[!rowSums(!is.finite(idr)),]  #remove rows that have -Inf
#========================================================== Write 'all' file
idrxAll=rbind(idrxKAf,idrxKAm,idrxKBf,idrxKBm,idrxKCf,idrxKCm)
#  colnames(idrxAll)=c('int','diff.sd','rate.sd.g','log.i','log|d|','log|r|','wgt')
nn;nrow(idrxAll)
writefile<-"C://R_aaROEVchapt08//8.1.16_Koontz&al2001_Dipodomys_out.csv"
write.csv(idrxAll,file=writefile,na="NA",row.names=FALSE) 

#----------------------------------------------------------
text(7,-2.6,paste("Processing time: ",round(proc.time()[3]-ptm[3],digits=2),
  "seconds"),pos=4,cex=1,col=4)
warnings()

