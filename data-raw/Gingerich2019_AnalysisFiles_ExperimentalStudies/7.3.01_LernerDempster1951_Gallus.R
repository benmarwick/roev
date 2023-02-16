##Philip D. Gingerich RATES OF EVOLUTION (Cambridge University Press 2019)
##7.3.01_LernerDempster1951_Gallus
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
xos=2.5;xoe=17.5;yos=9;yoe=18	#x,y original start and end
xfs=0;xfe=14;yfs=2.1;yfe=2.55		#x,y foreground start and end
xfm=(xoe-xos)/(xfe-xfs)				#x foreground multiplier
yfm=(yoe-yos)/(yfe-yfs)				#y foreground multiplier
lines(c(xos,xoe),c(yos,yos),lwd=1,lty=1,col=1)
lines(c(xos,xos),c(yos,yoe),lwd=1,lty=1,col=1)

#---------------------------------------------------------axis labels
for (i in seq(yos,yoe-1,1)){
	lines(c(xos-.12,xos),c(i,i),lwd=1,lty=1,col=1)
}
for (i in seq(yos,yoe-1,1)){
	text(xos,i,format(round(.05*(i-8)+2.05,digits=2),nsmall=2),pos=2,cex=1,col=1)
}
for (i in xos+c(xfs:(xfe/1))*1*xfm){
	lines(c(i,i),c(yos,yos-.1),lwd=1,lty=1,col=1)
}
for (i in xos+c(xfs:(xfe/1))*1*xfm){
	text(i,yos,(i-2.5)/xfm,pos=1,cex=1,col=1)
}
text(xos,yos+8.8,expression(paste(italic('Gallus gallus')*': shank length')),
	pos=4,cex=1.3,col=1)
text(xos-1.9,yos+4,'Ln shank length (cm)',
	srt=90,cex=1.3,col=1)
text(xos+7.5,yos-1.2,'Generation',
	cex=1.3,col=1)
#============================================================ Load file
file1<-"../Gingerich2019_AnalysisFiles_ExperimentalStudies/7.3.01_LernerDempster1951_Gallus.csv"
LD<-read.csv(file1,header=TRUE,sep=",",quote="\"",dec=".",fill=TRUE)
LD
#---------------------------------------------------------------------
print(c(min(LD[,2]),max(LD[,2]))) # generation
print(c(min(LD[,5]-LD[,6]),max(LD[,5]+LD[,6]))) # Cmean -/+  Csd
print(c(min(LD[,10]-LD[,11]),max(LD[,10]+LD[,11]))) # Smean  -/+ Ssd
lines(c(2.5,17.5),20*c(LD[1,5],LD[1,5])-33,lty=2,lwd=1.4,col=gray(6/10))
#----------------------------------------------- Plot control line [,3&5-6]
# Gen by
for (i in 2:15){
  lines(2.47+c(LD[i,2],LD[i,2])*xfm,
    20*c(LD[i,5]-LD[i,6],LD[i,5]+LD[i,6])-33,
    lty=1,lwd=1,col=gray(6/10))
  lines(2.5+c(LD[i-1,2],LD[i,2])*xfm,
    20*c(LD[i-1,5],LD[i,5])-33,
    lty=1,lwd=1,col=1)
}
#------------------------------------------- Plot selection line [,3&10-11]
for (i in 2:15){
  lines(2.53+c(LD[i,2],LD[i,2])*xfm,
    20*c(LD[i,10]-LD[i,11],LD[i,10]+LD[i,11])-33,
    lty=1,lwd=1,col=gray(6/10))
  lines(2.5+c(LD[i-1,2],LD[i,2])*xfm,
    20*c(LD[i-1,10],LD[i,10])-33,
    lty=1,lwd=1,col=1)
}
#------------------------------------------------Overlay points
points(2.5+LD[,2]*xfm,
  20*LD[,10]-33,
  pch=21,bg=1,cex=1.2,col=1)
points(2.5+LD[,2]*xfm,
  20*LD[,5]-33,
  pch=21,bg='white',cex=1.2,col=1)
#================================================ Rate calc LDselection
LD[1:3,]
n=length(LD[,2]);n
nn=.5*(n-1)*n;nn
idr=matrix(nrow=nn,ncol=9)
  colnames(idr)=c('int','diff.sd','rate.sd','log.i','log.d','log.r',
    'sbn','wgt','sum')
nc=0
stdev<-numeric(length=n-1)
for (k in 1:(n-1)){    #run length
  for (i in 1:(n-k)){  #starting position
    nc=nc+1
    idr[nc,1]=k							#intercept
    meandiff=LD[(i+k),10]-LD[i,10]				#mean diff.
    poolsd=PoolSD(LD[i+k,7],LD[i,7],LD[i+k,11],LD[i,11])
    idr[nc,2]=meandiff/poolsd					#diff.sd
    idr[nc,3]=idr[nc,2]/idr[nc,1]				#rate.sd
    idr[nc,4]=log10(idr[nc,1])				#log.i
    idr[nc,5]=log10(abs(idr[nc,2]))				#log.d
    idr[nc,6]=log10(abs(idr[nc,3]))				#log.r
    if(idr[nc,1]==1){idr[nc,7]=1}else{idr[nc,7]=3}	#sbn
    idr[nc,8]=1/idr[nc,1]					#wgt
    idr[nc,9]=idr[nc,4]+idr[nc,6]				#sum
    stdev[i]=poolsd
  }
}
vartot=sum(stdev^2,na.rm=TRUE);varcnt=length(table(stdev,useNA='no'))
varmean=vartot/varcnt; avestdev=sqrt(varmean);avestdev
idr[1:3,]
idrx1=idr[!rowSums(!is.finite(idr)),]  #remove rows that have -Inf
#writefile<-"C://R_aaRateBook//RB_Chapter07//LernerDempster_selection.csv"
#write.csv(idrx1,file=writefile,na="NA",row.names=FALSE)

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
bootresultd=TriPanelBC(idrx1,"r",1,4.5,bootn,mode,psize,"lower")
wrlb.d=bootresultd[[1]];bootmat.d=bootresultd[[2]]
proc.time()-ptm
#================================================ Rate calc LDcontrol
LD[1:3,]
n=length(LD[,2]);n
nn=.5*(n-1)*n;nn
idr=matrix(nrow=nn,ncol=9)
  colnames(idr)=c('int','diff.sd','rate.sd','log.i','log.d','log.r',
    'sbn','wgt','sum')
nc=0
stdev<-numeric(length=n-1)
for (k in 1:(n-1)){    #run length
  for (i in 1:(n-k)){  #starting position
    nc=nc+1
    idr[nc,1]=k							#intercept
    meandiff=LD[(i+k),5]-LD[i,5]				#mean diff.
    poolsd=PoolSD(LD[i+k,3],LD[i,3],LD[i+k,6],LD[i,6])
    idr[nc,2]=meandiff/poolsd					#diff.sd
    idr[nc,3]=idr[nc,2]/idr[nc,1]				#rate.sd
    idr[nc,4]=log10(idr[nc,1])				#log.i
    idr[nc,5]=log10(abs(idr[nc,2]))				#log.d
    idr[nc,6]=log10(abs(idr[nc,3]))				#log.r
    if(idr[nc,1]==1){idr[nc,7]=1}else{idr[nc,7]=3}	#sbn
    idr[nc,8]=1/idr[nc,1]					#wgt
    idr[nc,9]=idr[nc,4]+idr[nc,6]				#sum
    stdev[i]=poolsd
  }
}
vartot=sum(stdev^2,na.rm=TRUE);varcnt=length(table(stdev,useNA='no'))
varmean=vartot/varcnt; avestdev=sqrt(varmean);avestdev
idr[1:3,]
idrx2=idr[!rowSums(!is.finite(idr)),]  #remove rows that have -Inf
#writefile<-"C://R_aaRateBook//RB_Chapter07//LernerDempster_control.csv"
#write.csv(idrx2,file=writefile,na="NA",row.names=FALSE)
#========================================================== Write rate files
idrxA=rbind(idrx1,idrx2)
writefile<-"C://R_aaROEVchapt07//7.3.01_LernerDempster1951_Gallus_out.csv"
write.csv(idrxA,file=writefile,na="NA",row.names=FALSE)

#========================================================== Plot LRI panel c
#send idrx matrix, n, mode(diff/rate), panel placement coordinates, mode
bootresultr=TriPanelBC(idrx2,"r",13,4.5,bootn,mode,psize,"lower")
wrlb.r=bootresultr[[1]];bootmat.r=bootresultr[[2]]
proc.time()-ptm
#========================================================== Label B and C
rect(-.5,6.7,21.5,7.5,border=NA,col='white')
rect(.1,6,9,6.7,border=NA,col='white')
rect(12.1,6,21,6.7,border=NA,col='white')
text(0,6.5,"Selected line",pos=4,cex=1.3,col=1)
text(12,6.5,"Control population",pos=4,cex=1.3,col=1)

#----------------------------------------------------------
text(7,-2.6,paste("Processing time: ",round(proc.time()[3]-ptm[3],digits=2),
  "seconds"),pos=4,cex=1,col=4)
print(bootcountd);print(bootcountr)
#======================================================== Likelihoods
#Last Draw.Likelihoods number is a vertical scale multiplier relative to 800
#Draw.Likelihoods(bootmat.d,bootmat.r,1)


