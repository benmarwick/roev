##Philip D. Gingerich RATES OF EVOLUTION (Cambridge University Press 2019)
##9.2.32_Secord&al2012_Hyracotherium
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
	#mgp=c(2.2,.3,0),			#margin for axis title/labels/tickline
	mgp=c(2,.3,0),			#margin for axis title/labels/tickline
	tck=-0.01,				#tick-mark length
	xaxp=c(xr[1],xr[2],xr[1]+xr[2]),		#x-axis ticks, no. ticks
	yaxp=c(yr[1],yr[2],abs(yr[1]-yr[2])),	#y-axis ticks, no. ticks
	cex.axis=.8,
	cex.lab=1.1)			#label size
dev.size('in')
#============================================================ Load file
file1<-"../Gingerich2019_AnalysisFiles_FossilStudies/9.2.32_Secord&al2012_Hyracotherium.csv"
H<-read.csv(file1,header=TRUE,sep=",",quote="\"",dec=".",fill=TRUE)
H[1:5,];nrow(H);ncol(H)
#------------------------------------------------ Recognize samples
H[1:7,]
SH<-matrix(nrow=0,ncol=4);colnames(SH)=c('age','N','lnM','lnS')	#Stat H
st=H[1,12];nc=0		#strat level
samp<-matrix(nrow=0,ncol=4);colnames(samp)=c('age','n','lnm','lns')
for (i in 1:nrow(H)){
  if (H[i,12]==st){
    nc=nc+1;samp=rbind(samp,c(st,nc,H[i,10],H[i,10]))
    st=H[i,12]
  }else{
    #print(samp)
    if (nc>1){
      sh<-numeric(length=4)
      sh[1]=st;sh[2]=nrow(samp);sh[3]=mean(samp[,3]);sh[4]=sd(samp[,3])
      SH=rbind(SH,sh)
    }
    if (nc==1){
      #samp[i,4]=0
      SH=rbind(SH,samp)

    }
    samp<-matrix(nrow=0,ncol=4);colnames(samp)=c('age','n','lnm','lns')
    nc=1;samp=c(H[i,12],nc,H[i,10],H[i,10])
    st=H[i,12]
  }
}
samp<-matrix(nrow=0,ncol=4);colnames(samp)=c('age','n','lnm','lns')
nc=1;samp=c(st,nc,H[i,10],H[i,10])#;samp[i,4]=0
SH=rbind(SH,samp)
for (i in 1:nrow(SH)){
  if(SH[i,2]==1){SH[i,4]=NA}
};SH
wgtsd=sqrt(weighted.mean(SH[,4]^2,SH[,2],na.rm=TRUE))
SH[,4]=wgtsd;SH
#================================================================= Rates
idr=matrix(nrow=0,ncol=9)
  colnames(idr)=c('int.g','diff.sd','rate.sd.g','log.i','log|d|','log|r|',
    'sbn','wgt','li+lr')
#------------------------------------------------ Base rates only
SH[1:5,]
nr=nrow(SH)
zc=0							#zero count
nc=0
#i=1
for (i in 1:(nr-1)){				#increment
  id<-numeric(length=9)
  nc=nc+1
  iyr=175000/(48.04-14.45)*(SH[i,1]-SH[i+1,1])		#interval in years
  Mwm=(SH[i,2]*SH[i,3]+SH[i+1,2]*SH[i+1,3])/(SH[i,2]+SH[i+1,2])#molar wgt.mean
  wg=1000*exp(1.515*Mwm-3.261)
  gt=10^(.266*log10(wg)-.553)			#Eq 3.2 generation time in years
  id[1]=iyr/gt					#interval in generations
  meandiff=SH[i,3]-SH[i+1,3]			#mean diff.
  if (!is.na(meandiff)){if (meandiff==0){zc=zc+1}}
  psd=SH[i,4] #PoolSD(LH[sr+i,cn-1],LH[sr,cn-1],LH[sr+i,cn+1],LH[sr,cn+1])#n1,n2,sd1,sd2
  id[2]=meandiff/psd				#diff.sd
  id[3]=abs(id[2])/id[1]			#rate.sd.gen
  id[4]=log10(id[1])				#log.i
  id[5]=log10(abs(id[2]))			#log.d
  id[6]=log10(id[3])				#log.r
  id[7]=2						#sbn
  id[8]=1/id[1]					#wgt
  id[9]=id[4]+id[6]				#sum
  idr=rbind(idr,id)
}

#--------------------------------------------------------------------------
nrow(idr)
#--------------------------------------------------------------------------
idrx=idr[!rowSums(!is.finite(idr)),]  #remove rows that have -Inf
nrow(idr);nrow(idrx);zc
min(idrx[,9]);max(idrx[,9])

writefile<-"C://R_aaROEVchapt09//9.2.32_Secord&al2012_Hyracotherium_out.csv"
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
bootresultd=PalPanelBC(idrx[,1:8],"d",1,4.5,bootn,mode,psize,"normal",0)
proc.time()-ptm
#======================================================== Plot LRI panel (c)
#send idrx matrix, n, mode(diff/rate), panel placement coordinates, mode
bootresultr=PalPanelBC(idrx[,1:8],"r",13,4.5,bootn,mode,psize,"normal",0)
proc.time()-ptm

#----------------------------------------------------------
text(7,-2.6,paste("Processing time: ",round(proc.time()[3]-ptm[3],digits=2),
  "seconds"),pos=4,cex=1,col=4)
#print(bootcountd);print(bootcountr)



