##Philip D. Gingerich RATES OF EVOLUTION (Cambridge University Press 2019) 
##8.2.05_Baker&al1990_Fringilla
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
text(xos,yos+8.0,expression(paste(italic('Fringilla coelebs')*
	': body weight')),pos=4,cex=1.3,col=1)
text(xos-1.9,yos+4,'Ln weight (g)',
	srt=90,cex=1.3,col=1)
text(xos+7,yos-1.2,'Year',
	cex=1.3,col=1)
#============================================================ Load file
list.files("C://R_aaRateBook//RB_Chapter08",pattern=".csv")
file1<-"C://R_aaROEVchapt08//8.2.05_Baker&al1990_Fringilla.csv"
B<-read.csv(file1,header=TRUE,sep=",",quote="\"",dec=".",fill=TRUE)
B[1:3,]
LnB=B
for (i in 1:8){
  c1=3*i-1;c2=c1+1;c3=c1+2;print(c(c1,c2,c3))
  LnB[,c2]=log(B[,c2])
  LnB[,c3]=(sqrt(B[,c1])*B[,c3])/B[,c2]
}
LnB[1:3,]
#======================================================= Calculate rates
idr=matrix(nrow=12*28,ncol=8)
  colnames(idr)=c('int','diff.sd','rate.sd.g','log.i','log|d|','log|r|',
    'S-B-N','wgt')
nr=0
for (i in 1:7){					#initial locality
  for (j in (i+1):8){				#second locality
    if (j==8){int=90}else{int=120}	#2-year generation time;cross-sectional
    for (k in 1:12){
       nr=nr+1
       idr[nr,1]=int				#interval in gen
       diff=LnB[k,(3*i)]-LnB[k,(3*j)]
       wgtsumvar=LnB[k,(3*i-1)]*LnB[k,(3*i+1)]^2+LnB[k,(3*j-1)]*LnB[k,(3*j+1)]^2
       poolsd=sqrt(wgtsumvar/(LnB[k,(3*i-1)]+LnB[k,(3*j-1)]))
       idr[nr,2]=diff/poolsd			#diff.sd		
       idr[nr,3]=abs(idr[nr,2])/idr[nr,1]	#rate.sd.g
       idr[nr,4]=log10(idr[nr,1])
       idr[nr,5]=log10(abs(idr[nr,2]))
       idr[nr,6]=log10(idr[nr,3])
       idr[nr,7]=2
       idr[nr,8]=1/idr[nr,1]
    }
  }
}
nrow(idr)
idrx=idr[!rowSums(!is.finite(idr)),]  #remove rows that have -Inf
nrow(idrx)
writefile<-"C://R_aaROEVchapt08//8.2.05_Baker&al1990_Fringilla_out.csv"
write.csv(idrx,file=writefile,na="NA",row.names=FALSE) 

#----------------------------------------------------------
text(7,-2.6,paste("Processing time: ",round(proc.time()[3]-ptm[3],digits=2),
  "seconds"),pos=4,cex=1,col=4)
warnings()

