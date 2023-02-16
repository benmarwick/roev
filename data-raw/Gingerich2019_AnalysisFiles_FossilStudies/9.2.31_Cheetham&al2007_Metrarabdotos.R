##Philip D. Gingerich RATES OF EVOLUTION (Cambridge University Press 2019) 
##9.2.31_Cheetham&al2007_Metrarabdotos
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
	#asp=16*.5/1400,			#aspect ratio asp=x/y)
	col=1,				#color black
	las=1,				#axis labels always horizontal
	#mgp=c(2.2,.3,0),			#margin for axis title/labels/tickline
	mgp=c(2,.3,0),			#margin for axis title/labels/tickline
	tck=-0.01,				#tick-mark length
	xaxp=c(xr[1],xr[2],xr[1]+xr[2]),		#x-axis ticks, no. ticks
	yaxp=c(yr[1],yr[2],abs(yr[1]-yr[2])),	#y-axis ticks, no. ticks
	cex.axis=.8,
	cex.lab=1.1)			#label size
dev.size('in')
#============================================================ Load file
file1<-"C://R_aaROEVchapt09//9.2.31_Cheetham&al2007_Metrarabdotos.csv" 
C<-read.csv(file1,header=TRUE,sep=",",quote="\"",dec=".",fill=TRUE)
C[1:5,];nrow(C);ncol(C)
#======================================================= Plot panel a
xos=1;xoe=19;yos=9;yoe=18			#x,y original start and end
xfs=1.1;xfe=1.7;yfs=-56.2;yfe=53.8		#x,y foreground start and end
xfm=(xoe-xos)/(xfe-xfs)				#x foreground multiplier
yfm=(yoe-yos)/(yfe-yfs)				#y foreground multiplier

for (i in seq(xos,xos+15,5)){
  print(i)
  lines(c(i,i+3.5),c(yos,yos),lwd=1,lty=1,col=1)	#x lines
  for (j in 1:3){
    lines(c(i+j,i+j),c(yos,yos-.15),lwd=1,lty=1,col=1)
    text(i+j,yos,round(10*j-10,digits=0),pos=1,cex=1,col=1)
    text(i+2,yos-.5,"Mahalanobis D",pos=1,cex=1,col=1)
  }
  lines(c(i,i),c(yos,yoe),lwd=1,lty=1,col=1)		#y lines
  for (j in 1:9){
    lines(c(i,i-.15),c(9+j,9+j),lwd=1,lty=1,col=1)
    text(i,9+j,round(9-j,digits=0),pos=2,cex=1,col=1)
    text(i-1,yos+4.5,"Age (Ma)",srt=90,cex=1,col=1)
  }
}
#------------------------------------------------------ Plot saundersi-boldi
C[,1:6];sp=3;X12=18.6		#sp is initial species (2nd species add 3)
L1=sum(!is.na(C[,sp]));L1;L2=sum(!is.na(C[,sp+3]));L2
Dx1=mean(C[1:L1,sp]);Dx1;Dx2=mean(C[1:L2,sp+3]);Dx2
D1<-matrix(nrow=L1,ncol=2);D2<-matrix(nrow=L2,ncol=2)
D1[1:L1,1]=C[1:L1,sp-1];D1[1:L1,2]=C[1:L1,sp]-Dx1;mean(D1[,2])
D2[1:L2,1]=C[1:L2,sp+2];D2[1:L2,2]=-C[1:L2,sp+3]+Dx2+X12;mean(D2[,2])
lines(.1*c(D1[L1,2],D2[L2,2])+2,18-c(D1[L1,1],D2[L2,1]),lty=2,lwd=1,col=1)
for (i in 2:L1){
  lines(c(.1*D1[i-1,2]+2,.1*D1[i,2]+2),
  c(18-D1[i-1,1],18-D1[i,1]),lty=1,lwd=1,col=1)
}
points(.1*D1[,2]+2,18-D1[,1],pch=20,cex=.8,col=1)
text(.1*D1[1,2]+2-.3,18-D1[1,1]+.2,"M. saundersi",
  srt=90,pos=4,font=3,cex=.9,col=1)
for (i in 2:L2){
  lines(c(.1*D2[i-1,2]+2,.1*D2[i,2]+2),
  c(18-D2[i-1,1],18-D2[i,1]),lty=1,lwd=1,col=1)
}
points(.1*D2[,2]+2,18-D2[,1],pch=20,cex=.8,col=1)
text(.1*D2[1,2]+2-.3,18-D2[1,1]+.2,"M. boldi",
  srt=90,pos=4,font=3,cex=.9,col=1)

#------------------------------------------------------ Plot boldi-coatesi
C[,7:12];sp=9;X12=6.9		#sp is initial species (2nd species add 3)
L1=sum(!is.na(C[,sp]));L1;L2=sum(!is.na(C[,sp+3]));L2
Dx1=mean(C[1:L1,sp]);Dx1;Dx2=mean(C[1:L2,sp+3]);Dx2
D1<-matrix(nrow=L1,ncol=2);D2<-matrix(nrow=L2,ncol=2)
D1[1:L1,1]=C[1:L1,sp-1];D1[1:L1,2]=C[1:L1,sp]-Dx1;mean(D1[,2])
D2[1:L2,1]=C[1:L2,sp+2];D2[1:L2,2]=C[1:L2,sp+3]-Dx2+X12;mean(D2[,2])
lines(.1*c(D1[L1,2],D2[L2,2])+7,18-c(D1[L1,1],D2[L2,1]),lty=2,lwd=1,col=1)
for (i in 2:L1){
  lines(c(.1*D1[i-1,2]+7,.1*D1[i,2]+7),
  c(18-D1[i-1,1],18-D1[i,1]),lty=1,lwd=1,col=1)
}
points(.1*D1[,2]+7,18-D1[,1],pch=20,cex=.8,col=1)
text(.1*D1[1,2]+7-.3,18-D1[1,1]+.2,"M. boldi",
  srt=90,pos=4,font=3,cex=.9,col=1)
for (i in 2:L2){
  lines(c(.1*D2[i-1,2]+7,.1*D2[i,2]+7),
  c(18-D2[i-1,1],18-D2[i,1]),lty=1,lwd=1,col=1)
}
points(.1*D2[,2]+7,18-D2[,1],pch=20,cex=.8,col=1)
text(.1*D2[1,2]+7-.3,18-D2[1,1]+.2,"M. coatesi",
  srt=90,pos=4,font=3,cex=.9,col=1)
#----------------------------------------------- Plot colligatum-auriculatum
C[,7:12];sp=15;X12=17.2		#sp is initial species (2nd species add 3)
L1=sum(!is.na(C[,sp]));L1;L2=sum(!is.na(C[,sp+3]));L2
Dx1=mean(C[1:L1,sp]);Dx1;Dx2=mean(C[1:L2,sp+3]);Dx2
D1<-matrix(nrow=L1,ncol=2);D2<-matrix(nrow=L2,ncol=2)
D1[1:L1,1]=C[1:L1,sp-1];D1[1:L1,2]=-C[1:L1,sp]+Dx1;mean(D1[,2])
D2[1:L2,1]=C[1:L2,sp+2];D2[1:L2,2]=-C[1:L2,sp+3]+Dx2+X12;mean(D2[,2])
lines(.1*c(D1[L1,2],D2[L2,2])+12,18-c(D1[L1,1],D2[L2,1]),lty=2,lwd=1,col=1)
for (i in 2:L1){
  lines(.1*c(D1[i-1,2],D1[i,2])+12,
  18-c(D1[i-1,1],D1[i,1]),lty=1,lwd=1,col=1)
}
points(.1*D1[,2]+12,18-D1[,1],pch=20,cex=.8,col=1)
text(.1*D1[1,2]+12-.3,18-D1[1,1]+.2,"M. colligatum",
  srt=90,pos=4,font=3,cex=.9,col=1)
for (i in 2:L2){
  lines(.1*c(D2[i-1,2],D2[i,2])+12,
  18-c(D2[i-1,1],D2[i,1]),lty=1,lwd=1,col=1)
}
points(.1*D2[,2]+12,18-D2[,1],pch=20,cex=.8,col=1)
text(.1*D2[5,2]+12-.6,18-D2[5,1]+.2,"M. auriculatum",
  srt=90,pos=4,font=3,cex=.9,col=1)
#------------------------------------------------- Plot cubaguaense-vokesorum
C[,7:12];sp=21;X12=20.4		#sp is initial species (2nd species add 3)
L1=sum(!is.na(C[,sp]));L1;L2=sum(!is.na(C[,sp+3]));L2
Dx1=mean(C[1:L1,sp]);Dx1;Dx2=mean(C[1:L2,sp+3]);Dx2
D1<-matrix(nrow=L1,ncol=2);D2<-matrix(nrow=L2,ncol=2)
D1[1:L1,1]=C[1:L1,sp-1];D1[1:L1,2]=C[1:L1,sp]-Dx1;mean(D1[,2])
D2[1:L2,1]=C[1:L2,sp+2];D2[1:L2,2]=-C[1:L2,sp+3]+Dx2+X12;mean(D2[,2])
lines(.1*c(D1[L1,2],D2[L2,2])+17,18-c(D1[L1,1],D2[L2,1]),lty=2,lwd=1,col=1)
for (i in 2:L1){
  lines(.1*c(D1[i-1,2],D1[i,2])+17,
  18-c(D1[i-1,1],D1[i,1]),lty=1,lwd=1,col=1)
}
points(.1*D1[,2]+17,18-D1[,1],pch=20,cex=.8,col=1)
text(.1*D1[1,2]+17-.3,18-D1[1,1]+.2,"M. vokesorum",
  srt=90,pos=4,font=3,cex=.9,col=1)
points(.1*D2[,2]+17,18-D2[,1],pch=20,cex=.8,col=1)
text(.1*D2[1,2]+17-.3,18-D2[1,1]+.2,"M. cubaguaense",
  srt=90,pos=4,font=3,cex=.9,col=1)

#===================================================== Rates: seven lineages
idr=matrix(nrow=0,ncol=9)
  colnames(idr)=c('int.g','diff.sd','rate.sd.g','log.i','log|d|','log|r|',
    'sbn','wgt','li+lr')
#------------------------------------------------ Rates
C[1:3,]
gt=10			#generation time from Cheetham and Jackson 1995
zc=0							#zero count
nc=0
for (cn in seq(3,21,3)){		#columns of means
  nr=sum(!is.na(C[,cn]))
  for (i in 1:(nr-1)){				#increment
    for (sr in 1:(nr-i)){			#start row
      id<-numeric(length=9)
      nc=nc+1
      iyr=1000000*(C[sr+i,cn-1]-C[sr,cn-1])		#interval in years
      id[1]=iyr/gt					#interval in generations
      meandiff=C[sr+i,cn]-C[sr,cn]		#mean diff.
      if (!is.na(meandiff)){if (meandiff==0){zc=zc+1}}
      id[2]=meandiff					#diff.sd
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
idrx=idr[!rowSums(!is.finite(idr)),]  #remove rows that have -Inf
nrow(idr);nrow(idrx);zc
min(idrx[,9]);max(idrx[,9])
writefile<-"C://R_aaROEVchapt09//9.2.31_Cheetham&al2007_Metrarabdotos_out.csv"
write.csv(idrx,file=writefile,na="NA",row.names=FALSE) 

#======================================================== Plot LRI panel (e)
#assign 'bootn' as boot number
bootn<-1000;text(-1.7,-2.6,paste("Boot n =",bootn),pos=4,cex=1,col=4)
#assign 'mode' as "medians","all","mixed"
mode<-"all";text(2,-2.6,paste("Mode: ",mode),pos=4,cex=1,col=4)
#assign circle size for points (1.5 or 2)
psize<-1.2	#1.5/2	
#assign 'equation' position as "normal","lower","none" at end of each call
#send (1)idrx matrix, (2)mode(diff/rate), (3)panel placement coordinate x, 
#  (4)panel placement coordinate y, (5)bootn, (6)mode, (7)psize, (8)equation
bootresultd=PalPanelBC(idrx[,1:8],"d",1,4.5,bootn,mode,psize,"normal",0)
proc.time()-ptm
#======================================================== Plot LRI panel (f)
#send idrx matrix, n, mode(diff/rate), panel placement coordinates, mode
bootresultr=PalPanelBC(idrx[,1:8],"r",13,4.5,bootn,mode,psize,"normal",0)
proc.time()-ptm

#==================================================== Rates between lineages
C[1:3,]
nc=0
rab<-matrix(nrow=4,ncol=4)						#rate a to b
  colnames(rab)=c('age1','D1','age2','D2')
for (cn in seq(3,21,6)){
  nc=nc+1
  nr1=sum(!is.na(C[,cn]));nr1
  rab[nc,1]=C[nr1,cn-1];rab[nc,2]=C[nr1,cn]
  nr2=sum(!is.na(C[,cn+3]));nr2
  rab[nc,3]=C[nr2,cn+2];rab[nc,4]=C[nr2,cn+3]
};rab
rab[,2]=c(0,0,0,0)
rab[,4]=c(18.6,6.9,17.2,20.4);rab
bet=matrix(nrow=4,ncol=9)			#Between species matrix to append
  colnames(bet)=c('int.g','diff.sd','rate.sd.g','log.i','log|d|','log|r|',
    'sbn','wgt','li+lr')
for (i in 1:4){
  bet[i,1]=100000*(rab[i,1]-rab[i,3])		#int.g
  bet[i,2]=rab[i,4]				#diff.sd
  bet[i,3]=abs(bet[i,2])/bet[i,1]		#rate.sd.g
  bet[i,4]=log10(bet[i,1])			#log int
  bet[i,5]=log10(abs(bet[i,2]))			#log diff
  bet[i,6]=log10(bet[i,3])			#log rate
  bet[i,7]=2					#sbn
  bet[i,8]=1/bet[i,1]				#wgt
  bet[i,9]=bet[i,4]+bet[i,6]			#LI+LR
};bet
points(bet[,4]+1,bet[,5]+4.5,pch=19,cex=1,col=1)
text(9.8,5.6,expression(paste("Between spp.")),pos=2,cex=.9,col=1)
points(bet[,4]+13,bet[,6]+6.5,pch=19,cex=1,col=1)
text(20.9,3.9,expression(paste("Between species")),pos=2,cex=.9,col=1)
writefile<-"C://R_aaROEVchapt09//9.2.31_Cheetham&al2007_Metrarabdotos2_out.csv"
write.csv(bet[,1:8],file=writefile,na="NA",row.names=FALSE) 
 
#----------------------------------------------------------
text(7,-2.6,paste("Processing time: ",round(proc.time()[3]-ptm[3],digits=2),
  "seconds"),pos=4,cex=1,col=4)
#print(bootcountd);print(bootcountr)



