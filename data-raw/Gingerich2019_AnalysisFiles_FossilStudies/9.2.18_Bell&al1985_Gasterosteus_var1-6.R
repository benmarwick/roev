##Philip D. Gingerich RATES OF EVOLUTION (Cambridge University Press 2019) 
##9.2.18_Bell&al1985_Gasterosteus_var1-6
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
file1<-"C://R_aaROEVchapt09//9.2.18_Bell&al1985_Gasterosteus.csv" 
B<-read.csv(file1,header=TRUE,sep=",",quote="\"",dec=".",fill=TRUE)
B[1:5,];nrow(B);ncol(B)
#------------------------------------------------------Logged matrix
BB=apply(as.matrix.noquote(B[1:nrow(B),]),2,as.numeric)
BB[1:5,]
LB<-matrix(nrow=nrow(BB),ncol=19)
  colnames(LB)=colnames(BB)[2:20]
LB[,1]=BB[,2]						#ages
LB[,seq(2,17,3)]=BB[,seq(3,18,3)]			#Ns
LB[,3]=BB[,4]						#untransformed pelvic score
LB[,4]=BB[,5]						#untransformed score stdev
LB[,seq(6,18,3)]=log(BB[,seq(7,19,3)])			#log transformed means
LB[,seq(7,19,3)]=BB[,seq(8,20,3)]/BB[,seq(7,19,3)]	#log trans. stdevs
BB[1:5,]
LB[1:5,]
min(LB[,1]);max(LB[,1])
min(LB[,3])-min(LB[,4])
max(LB[,3])+max(LB[,4])

#======================================================= Plot panel a
xos=3;xoe=18;yos=9.5;yoe=16.5			#x,y original start and end
xfs=.3;xfe=3.3;yfs=-6;yfe=1			#x,y foreground start and end
xfm=(xoe-xos)/(xfe-xfs)				#x foreground multiplier
yfm=(yoe-yos)/(yfe-yfs)				#y foreground multiplier
lines(c(xos,xoe),c(yos,yos),lwd=1,lty=1,col=1)
lines(c(xos,xos),c(yos,yoe),lwd=1,lty=1,col=1)
#---------------------------------------------------------axis labels
for (i in seq(yos+.5,yoe-1,1)){				#y-axis
  lines(c(xos-.15,xos),c(i,i),lwd=1,lty=1,col=1)
}
for (i in seq(yos+.5,yoe-1,1)){
  text(xos+.1,i,format(round(abs(20*(i)-200),digits=0),nsmall=0),pos=2,cex=1,col=1)
}
for (i in xos-1.5+5*(seq(xfs,xfe,.1))){		#x-axis
  lines(c(i,i),c(yos,yos-.15),lwd=1,lty=1,col=1)
}
for (i in xos+5*(seq(.5,3,.5))){		#x-axis
  text(i-1.5,yos,format(round(.2*i-.6,digits=1),nsmall=1),
    pos=1,cex=1,col=1)
}
text(xos,yos+7.0,expression(paste(italic('Gasterosteus doryssus')~
	':  pelvic score')),pos=4,cex=1.3,col=1)
text(xos-1.5,yos+3,'Time (kyr)',
	srt=90,cex=1.3,col=1)
text(xos+7.5,yos-1.2,
  expression(paste('Pelvic score (0-3)')),cex=1.3,col=1)
#---------------------------------------------------------------- Plot ranges
LB[1:5,]
min(na.omit(LB[,3]));max(na.omit(LB[,3]))
min(LB[,1]);max(LB[,1])
xpos=-1.5;ypos=20
for (i in 1:nrow(LB)){
  lines(xos+5*c(LB[i,3]-LB[i,4],LB[i,3]+LB[i,4])+xpos,
    yos+.5+.00005*c(LB[i,1],LB[i,1]),lty=1,lwd=1,col=gray(6/10))
}
for (i in 1:nrow(LB)){
  points(xos+5*LB[i,3]+xpos,yos+.5+.00005*LB[i,1],pch=20,cex=.8,col=1)
}
lines(xos+5*LB[,3]+xpos,yos+.5+.00005*LB[,1],lty=1,lwd=1.5,col=1)
#==========================================================================
LB[,1:4]
SP2=LB[21:26,];SP2
SPt=LB[20:21,];SPt
SP1=LB[1:20,];SP1
#================================================ Rates SP1 G. doryssus
idr=matrix(nrow=0,ncol=9)
  colnames(idr)=c('int.g','diff.sd','rate.sd.g','log.i','log|d|','log|r|',
    'sbn','wgt','li+lr')
#------------------------------------------------ Rates
gt=2	#generation time assumed
nr=nrow(SP1)						#number of rows
nn1=.5*(nr-1)*nr;nn1
zc1=0							#zero count
nc=0
for (cn in seq(3,ncol(SP1),3)){		#columns of means
  for (i in 1:(nr-1)){				#increment
    for (sr in 1:(nr-i)){			#start row
      id<-numeric(length=9)
      nc=nc+1
      iyr=SP1[sr+i,1]-SP1[sr,1]		#interval in years
      id[1]=iyr/gt					#interval in generations
      meandiff=SP1[sr+i,cn]-SP1[sr,cn]		#mean diff.
      if (!is.na(meandiff)){if (meandiff==0){zc1=zc1+1}}
      psd=PoolSD(SP1[sr+i,cn-1],SP1[sr,cn-1],SP1[sr+i,cn+1],SP1[sr,cn+1])#n1,n2,sd1,sd2
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
idrx1=idr[!rowSums(!is.finite(idr)),]  #remove rows that have -Inf
nn1;nrow(idr);nrow(idrx1);zc1
min(idrx1[,9]);max(idrx1[,9])

#====================================================== Rates in transition
idr=matrix(nrow=0,ncol=9)
  colnames(idr)=c('int.g','diff.sd','rate.sd.g','log.i','log|d|','log|r|',
    'sbn','wgt','li+lr')
#------------------------------------------------------ Rates
gt=2	#generation time assumed
nr=nrow(SPt)					#number of rows
nnt=.5*(nr-1)*nr;nnt
zct=0							#zero n count
nc=0
for (cn in seq(3,ncol(SPt),3)){		#columns of means
  for (i in 1:(nr-1)){				#increment
    for (sr in 1:(nr-i)){			#start row
      id<-numeric(length=9)
      nc=nc+1
      iyr=SPt[sr+i,1]-SPt[sr,1]		#interval in years
      id[1]=iyr/gt;tr_int_g=id[1]			#interval in generations
      meandiff=SPt[sr+i,cn]-SPt[sr,cn]		#mean diff.
      if (!is.na(meandiff)){if (meandiff==0){zct=zct+1}}
      psd=PoolSD(SPt[sr+i,cn-1],SPt[sr,cn-1],SPt[sr+i,cn+1],SPt[sr,cn+1])#n1,n2,sd1,sd2
      id[2]=meandiff/psd;tr_psd_sd=id[2]		#diff.sd
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
idrxt=idr[!rowSums(!is.finite(idr)),]  #remove rows that have -Inf
nnt;nrow(idr);nrow(idrxt);zct
min(na.omit(idrxt[,9]));na.omit(max(idrxt[,9]))

#================================================ Rates SP2 
idr=matrix(nrow=0,ncol=9)
  colnames(idr)=c('int.g','diff.sd','rate.sd.g','log.i','log|d|','log|r|',
    'sbn','wgt','li+lr')
#------------------------------------------------ Rates
gt=2	#generation time assumed
nr=nrow(SP2)						#number of rows
nn2=.5*(nr-1)*nr;nn2
zc2=0							#zero n count
nc=0
for (cn in seq(3,ncol(SP2),3)){		#columns of means
  for (i in 1:(nr-1)){				#increment
    for (sr in 1:(nr-i)){			#start row
      id<-numeric(length=9)
      nc=nc+1
      iyr=SP2[sr+i,1]-SP2[sr,1]		#interval in years
      id[1]=iyr/gt					#interval in generations
      meandiff=SP2[sr+i,cn]-SP2[sr,cn]		#mean diff.
      if (!is.na(meandiff)){if (meandiff==0){zc2=zc2+1}}
      psd=PoolSD(SP2[sr+i,cn-1],SP2[sr,cn-1],SP2[sr+i,cn+1],SP2[sr,cn+1])#n1,n2,sd1,sd2
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
idrx2=idr[!rowSums(!is.finite(idr)),]  #remove rows that have -Inf
nn2;nrow(idr);nrow(idrx2);zc2
min(idrx2[,9]);max(idrx2[,9])

#--------------------------------------------------------------------Summary
idrxa=rbind(idrx1,idrxt,idrx2);idrx12=rbind(idrx1,idrx2)
nrow(idrx1);nrow(idrxt);nrow(idrx2);nrow(idrx12);nrow(idrxa)
zc1;zct;zc2;zc1+zct+zc2
zca=zc1+zct+zc2
zca;nrow(idrxa);zca+nrow(idrxa);nn1+nnt+nn2
writefile<-"C://R_aaROEVchapt09//9.2.18_Bell&al1985_Gasterosteus_var1-6_out.csv"
write.csv(idrxa,file=writefile,na="NA",row.names=FALSE) 

#======================================================== Add idrxt to plots
text(8.5,14.3,paste('highest rate: ',round(tr_psd_sd,digits=2),
  'st. dev. in',round(tr_int_g,digits=0),'gen.'),
  pos=4,cex=1,col=1)

#----------------------------------------------------------
text(7,-2.6,paste("Processing time: ",round(proc.time()[3]-ptm[3],digits=2),
  "seconds"),pos=4,cex=1,col=4)



