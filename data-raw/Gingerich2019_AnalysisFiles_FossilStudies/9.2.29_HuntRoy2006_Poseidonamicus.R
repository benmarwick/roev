##Philip D. Gingerich RATES OF EVOLUTION (Cambridge University Press 2019) 
##9.2.29_HuntRoy2006_Poseidonamicus
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
file1<-"C://R_aaROEVchapt09//9.2.29_HuntRoy2006_Poseidonamicus.csv" 
H<-read.csv(file1,header=TRUE,sep=",",quote="\"",dec=".",fill=TRUE)
H[1:5,];nrow(H);ncol(H)
#------------------------------------------------ Recognize species
H[1:5,]
sp=H[1,2]
spn=1
LH<-matrix(nrow=nrow(H),ncol=4);colnames(LH)=c('sp','age','lnM','lnS')
for (i in 1:nrow(H)){
  if (H[i,2]==sp){
    LH[i,1]=spn
    LH[i,2]=H[i,3]
    LH[i,3]=log(H[i,7])
    LH[i,4]=.03
  }else{
    sp=H[i,2]
    spn=spn+1
    LH[i,1]=spn
    LH[i,2]=H[i,3]
    LH[i,3]=log(H[i,7])
    LH[i,4]=.03
  }
}
LH

#================================================ Rates 
idr=matrix(nrow=0,ncol=9)
  colnames(idr)=c('int.g','diff.sd','rate.sd.g','log.i','log|d|','log|r|',
    'sbn','wgt','li+lr')
#------------------------------------------------ 
LH[1:5,]
gt=5	#generation time from Hunt 2007
zc=0							#zero count
nc=0
for (sp in 1:max(LH[,1])){
  nr=nrow(LH[LH[,1]==sp,]);nr				#number of rows
  TLH<-matrix(nrow=nr,ncol=ncol(LH))		#temp. subset
  TLH=LH[LH[,1]==sp,]
  nn=.5*(nr-1)*nr;nn
  for (i in 1:(nr-1)){					#increment
    for (sr in 1:(nr-i)){				#start row
      id<-numeric(length=9)
      nc=nc+1
      iyr=1000000*(TLH[sr,2]-TLH[sr+i,2])		#interval in years
      id[1]=iyr/gt					#interval in generations
      meandiff=TLH[sr+i,3]-TLH[sr,3]			#mean diff.
      if (!is.na(meandiff)){if (meandiff==0){zc=zc+1}}
      psd=.03 #PoolSD(TLH[sr+i,cn-1],TLH[sr,cn-1],TLH[sr+i,cn+1],TLH[sr,cn+1])#n1,n2,sd1,sd2
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
idrx=idr[!rowSums(!is.finite(idr)),]  #remove rows that have -Inf
nrow(idr);nrow(idrx);zc
min(idrx[,9]);max(idrx[,9])

writefile<-"C://R_aaROEVchapt09//9.2.29_HuntRoy2006_Poseidonamicus_out.csv"
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
bootresultr=PalPanelBC(idrx[,1:8],"r",13,4.5,bootn,mode,psize,"normal",0)
proc.time()-ptm

#----------------------------------------------------------
text(7,-2.6,paste("Processing time: ",round(proc.time()[3]-ptm[3],digits=2),
  "seconds"),pos=4,cex=1,col=4)
