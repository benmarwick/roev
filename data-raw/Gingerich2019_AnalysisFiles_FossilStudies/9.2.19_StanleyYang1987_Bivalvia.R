##Philip D. Gingerich RATES OF EVOLUTION (Cambridge University Press 2019) 
##9.2.19_StanleyYang1987_Bivalvia
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
file1<-"C://R_aaROEVchapt09//9.2.19_StanleyYang1987_Bivalvia.csv" 
S<-read.csv(file1,header=TRUE,sep=",",quote="\"",dec=".",fill=TRUE)
S[1:5,];nrow(S)
#================================================ Rates 
idr=matrix(nrow=0,ncol=9)
  colnames(idr)=c('int.g','diff.sd','rate.sd.g','log.i','log|d|','log|r|',
    'sbn','wgt','li+lr')
#------------------------------------------------ Rates
gt=5	#generation time assumed
nr=nrow(S);nr						#number of rows
nn=.5*(nr-1)*nr;nn
nz=0							#zero n count
nc=0
for (i in 1:nrow(S)){
  id<-numeric(length=9)
  nc=nc+1
  id[1]=S[i,1]*1000000/gt	#int_g
  id[2]=S[i,2]			#diff_sd
  id[3]=abs(id[2])/id[1]
  id[4]=log10(id[1])
  id[5]=log10(abs(id[2]))
  id[6]=log10(id[3])
  id[7]=NA
  id[8]=1/id[1]
  id[9]=id[4]+id[6]
  idr=rbind(idr,id)
}
#--------------------------------------------------------------------------
nrow(idr)
#--------------------------------------------------------------------------
idrx=idr
nn;nrow(idr);nrow(idrx);nz
min(idrx[,9]);max(idrx[,9])
writefile<-"C://R_aaROEVchapt09//9.2.19_StanleyYang1987_Bivalvia_out.csv"
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
xos=1;yos=-.5
text(xos-2.8,yos+8.2,expression(paste('Neogene Bivalvia: Mahalanobis D')),
  pos=4,cex=1.3,col=1)
#----------------------------------------------------------
text(7,-2.6,paste("Processing time: ",round(proc.time()[3]-ptm[3],digits=2),
  "seconds"),pos=4,cex=1,col=4)



