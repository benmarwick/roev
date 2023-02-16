##Philip D. Gingerich RATES OF EVOLUTION (Cambridge University Press 2019) 
##8.1.28_Babin-Fenske&al2008_Gyretes
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
xos=1;xoe=xos+13.5;yos=9;yoe=17			#x,y original start and end
xfs=0;xfe=13;yfs=5;yfe=5.4			#x,y foreground start and end
xfm=(xoe-xos)/(xfe-xfs)				#x foreground multiplier
yfm=(yoe-yos)/(yfe-yfs)				#y foreground multiplier
lines(c(xos,xoe),c(yos,yos),lwd=1,lty=1,col=1)
lines(c(xos,xos),c(yos,yoe),lwd=1,lty=1,col=1)

#---------------------------------------------------------axis labels
for (i in seq(yos,yoe,1)){
  lines(c(xos-.12,xos),c(i,i),lwd=1,lty=1,col=1)
}
for (i in seq(yos,yoe,2)){
  #text(xos,i+1,format(round(.05*(i-7)+7.3,digits=2),nsmall=2),pos=2,cex=1,col=1)
  text(xos,i,format(round(.025*(i-7)+1.5,digits=2),nsmall=2),pos=2,cex=1,col=1)
}
for (i in xos+c(xfs:xfe)*1*xfm){
  lines(c(i,i),c(yos,yos-.1),lwd=1,lty=1,col=1)
}
for (i in seq(1930,1990,10)){
  text(xos+(i-1930+5)*xfm/5,yos,i,pos=1,cex=1,col=1)
}
text(xos,yos+8.5,expression(paste(italic('Gyretes sinuatus')*
	': body length')),pos=4,cex=1.3,col=1)
text(xos-1.9,yos+4,'Ln length (mm)',
	srt=90,cex=1.3,col=1)
text(xos+7,yos-1.2,'Year',
	cex=1.3,col=1)
#============================================================ Load file
list.files("C://R_aaRateBook//RB_Chapter08",pattern=".csv")
file1<-"C://R_aaROEVchapt08//8.1.28_Babin-Fenske&al2008_Gyretes.csv"
B<-read.csv(file1,header=TRUE,sep=",",quote="\"",dec=".",fill=TRUE)
B[1:3,]
#-------------------------------------------------Stats - five groups
Gr<-matrix(nrow=5,ncol=5)
  colnames(Gr)=c('start','end','syear','eyear','year')
Gr[1,1:5]=c( 1,11,1928,1932,1930)
Gr[2,1:5]=c(12,23,1950,1953,1951)
Gr[3,1:5]=c(30,39,1971,1973,1972)
Gr[4,1:5]=c(42,69,1979,1982,1981)
Gr[5,1:5]=c(70,79,1986,1988,1987);Gr
LL<-matrix(nrow=5,ncol=7)		#---------------------Length meas.
  colnames(LL)=c('Gen','Yr','N','Mean','Stdev','Ln.M','Ln.SD')
for (i in 1:5){
  LL[i,1]=Gr[i,5]-Gr[1,5]
  LL[i,2]=Gr[i,5]
  LL[i,3]=length(B[Gr[i,1]:Gr[i,2],2][!is.na(B[Gr[i,1]:Gr[i,2],2])])
  LL[i,4]=mean(B[Gr[i,1]:Gr[i,2],2],na.rm=TRUE)
  LL[i,5]=sd(B[Gr[i,1]:Gr[i,2],2],na.rm=TRUE)
  LL[i,6]=log(LL[i,4])
  LL[i,7]=LL[i,5]/LL[i,4]
};LL
LW<-matrix(nrow=5,ncol=7)		#---------------------Width meas.
  colnames(LW)=c('Gen','Yr','N','Mean','Stdev','Ln.M','Ln.SD')
for (i in 1:5){
  LW[i,1]=Gr[i,5]-Gr[1,5]
  LW[i,2]=Gr[i,5]
  LW[i,3]=length(B[Gr[i,1]:Gr[i,2],3][!is.na(B[Gr[i,1]:Gr[i,2],2])])
  LW[i,4]=mean(B[Gr[i,1]:Gr[i,2],3],na.rm=TRUE)
  LW[i,5]=sd(B[Gr[i,1]:Gr[i,2],3],na.rm=TRUE)
  LW[i,6]=log(LW[i,4])
  LW[i,7]=LW[i,5]/LW[i,4]
};LW
LF<-matrix(nrow=5,ncol=7)		#---------------------Fineness
  colnames(LF)=c('Gen','Yr','N','Mean','Stdev','Ln.M','Ln.SD')
for (i in 1:5){
  LF[i,1]=Gr[i,5]-Gr[1,5]
  LF[i,2]=Gr[i,5]
  LF[i,3]=length(B[Gr[i,1]:Gr[i,2],4][!is.na(B[Gr[i,1]:Gr[i,2],2])])
  LF[i,4]=mean(B[Gr[i,1]:Gr[i,2],4],na.rm=TRUE)
  LF[i,5]=sd(B[Gr[i,1]:Gr[i,2],4],na.rm=TRUE)
  LF[i,6]=log(LF[i,4])
  LF[i,7]=LF[i,5]/LF[i,4]
};LF
#------------------------------------------------Plot measurements
bx=xos+(B[,1]-1930+5)*.2*xfm;by=40*log(B[,2])-53
points(bx,by,pch=19,cex=.8,col=gray(5/10))
min(bx);max(bx);min(by);max(by)
min(B[,1]);max(B[,1]);min(B[,2]);max(B[,2])
lmB2=lm(log(B[,2])~B[,1])
b=coef(lmB2)[1];m=coef(lmB2)[2]	#slope m and intercept b
rsx1=xos+(min(B[,1])-1930+5)*.2*xfm	#regression line coordinates
rsx2=xos+(max(B[,1])-1930+5)*.2*xfm
rsy1=40*(m*min(B[,1])+b)-53
rsy2=40*(m*max(B[,1])+b)-53
lines(c(rsx1,rsx2),c(rsy1,rsy2),lty=2,lwd=2,col=gray(5/10))
text(xos+1,yos+7,bquote(italic(Y)==.(format(m,nsmall=2,digits=2))~italic(X)~-~
  .(format(abs(b),nsmall=2,digits=4))),pos=4,cex=1.1,col=gray(4/10))

text(15,16,"Years:",pos=4,cex=1,col=1)
text(21.5,16,paste(format(max(B[,1])-min(B[,1]),digits=3,nsmall=3)),
  pos=2,cex=1,col=1)
text(15,15.3,"Interval (gen.):",pos=4,cex=1,col=1)
text(21.5,15.3,paste(format(max(B[,1])-min(B[,1]),digits=3,nsmall=3)),
  pos=2,cex=1,col=1)

text(15,14.3,bquote(italic(Y)[1988]~-~italic(Y)[1928]~":"),pos=4,cex=1,col=1)
text(21.5,14.3,format((m*max(B[,1])+b)-(m*min(B[,1])+b),digits=3,nsmall=4),
  pos=2,cex=1,col=1)
text(15,13.6,"Unexpl. MS:",pos=4,cex=1,col=1)
text(21.5,13.6,paste(format(anova(lmB2)[2,3],digits=2,nsmall=3)),
  pos=2,cex=1,col=1)
text(15,12.9,"Std. dev.:",pos=4,cex=1,col=1)
text(21.5,12.9,format(sqrt(anova(lmB2)[2,3]),digits=3,nsmall=3),
  pos=2,cex=1,col=1)
D=((m*max(B[,1])+b)-(m*min(B[,1])+b))/sqrt(anova(lmB2)[2,3]);D
text(15,12.2,"Diff. (std. dev.):",pos=4,cex=1,col=1)
text(21.5,12.2,format(D,digits=3,nsmall=4),pos=2,cex=1,col=1)

R=abs(D)/(max(B[,1])-min(B[,1]));R
text(15,11.2,"Rate:",pos=4,cex=1,col=1)
subscript=format(log10(max(B[,1])-min(B[,1])),digits=3,nsmall=3)
text(21.5,11.2,bquote(italic(h)[.(subscript)]==.(format(R,digits=3,nsmall=3))),
  pos=2,cex=1,col=1)

#Standard deviation is square root of mean square error unexplained by the
#  regression.  This is also the standard deviation of the residuals.

#------------------------------------------------Plot lengths
for (i in 1:nrow(LL)){
  lines(xos+c(LL[i,1]+5,LL[i,1]+5)*.2*xfm,
    40*c(LL[i,6]-LL[i,7],LL[i,6]+LL[i,7])-53,
    lty=1,lwd=1,col=1)
}
for (i in 2:nrow(LL)){
  lines(xos+c(LL[i-1,1]+5,LL[i,1]+5)*.2*xfm,40*c(LL[i-1,6],LL[i,6])-53,
	lty=1,lwd=1,col=1)
}
for (i in 1:nrow(LL)){
  points(xos+(LL[i,1]+5)*.2*xfm,40*LL[i,6]-53,pch=21,bg='white',cex=1.4,col=1)
}
#========================================== Rate calc: lengths
gentime=1
rect(xos+.1,9.2,9,9.6,col='white',border=NA)
text(xos,9.4,"Generation time: 1 year",pos=4,cex=.9,col=1)
LL
n=nrow(LL);n
nn=.5*(n-1)*n;nn
idr=matrix(nrow=nn,ncol=7)
  colnames(idr)=c('int','diff.sd','rate.sd.g','log.i','log|d|','log|r|','wgt')
nc=0
stdev<-numeric(length=n-1)
for (k in gentime:(n-1)){	#run length
  for (i in 1:(n-k)){  						#starting position
    nc=nc+1
    idr[nc,1]=(LL[(i+k),2]-LL[i,2])/gentime		#interval
    meandiff=LL[(i+k),6]-LL[i,6]				#mean diff.
    poolsd=PoolSD(LL[i+k,3],LL[i,3],LL[i+k,7],LL[i,7])#n1,n2,sd1,sd2
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
idrxL=idr[!rowSums(!is.finite(idr)),]  #remove rows that have -Inf
#============================================= Rate calc: widths
LW[1:3,]
n=nrow(LW);n
nn=.5*(n-1)*n;nn
idr=matrix(nrow=nn,ncol=7)
  colnames(idr)=c('int','diff.sd','rate.sd.g','log.i','log|d|','log|r|','wgt')
nc=0
stdev<-numeric(length=n-1)
for (k in gentime:(n-1)){	#run length
  for (i in 1:(n-k)){  						#starting position
    nc=nc+1
    idr[nc,1]=(LW[(i+k),2]-LW[i,2])/gentime		#intercept
    meandiff=LW[(i+k),6]-LW[i,6]				#mean diff.
    poolsd=PoolSD(LW[i+k,3],LW[i,3],LW[i+k,7],LW[i,7])#n1,n2,sd1,sd2
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
idrxW=idr[!rowSums(!is.finite(idr)),]  #remove rows that have -Inf
#========================================= Rate calc: fineness
LF
n=nrow(LF);n
nn=.5*(n-1)*n;nn
idr=matrix(nrow=nn,ncol=7)
  colnames(idr)=c('int','diff.sd','rate.sd.g','log.i','log|d|','log|r|','wgt')
nc=0
stdev<-numeric(length=n-1)
for (k in gentime:(n-1)){	#run length
  for (i in 1:(n-k)){  						#starting position
    nc=nc+1
    idr[nc,1]=(LF[(i+k),2]-LF[i,2])/gentime		#intercept
    meandiff=LF[(i+k),6]-LF[i,6]				#mean diff.
    poolsd=PoolSD(LF[i+k,3],LF[i,3],LF[i+k,7],LF[i,7])#n1,n2,sd1,sd2
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
idrxF=idr[!rowSums(!is.finite(idr)),]  #remove rows that have -Inf
#========================================================== Write 'all' file
idrxA=rbind(idrxL,idrxW,idrxF)
  colnames(idrxA)=c('int','diff.sd','rate.sd.g','log.i','log|d|','log|r|','wgt')
writefile<-"C://R_aaROEVchapt08//8.1.28_Babin-Fenske&al2008_Gyretes_out.csv"
write.csv(idrxA,file=writefile,na="NA",row.names=FALSE) 

#----------------------------------------------------------
text(7,-2.6,paste("Processing time: ",round(proc.time()[3]-ptm[3],digits=2),
  "seconds"),pos=4,cex=1,col=4)
warnings()

