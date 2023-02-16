

#=====================================================================Setup
library(RATES);library(MASS)
library(ggplot2)
#----------------------------------------------------------------------Plot


#============================================================ Load file
file1<-"../Gingerich2019_AnalysisFiles_FieldStudies/eutotal.csv"
LD<-read.csv(file1,header=TRUE,sep=",",quote="\"",dec=".",fill=TRUE)
LD <- LD[complete.cases(LD), ]
LD
#---------------------------------------------------------------------
print(c(min(LD[,1]),max(LD[,1]))) # year
print(c(min(LD[,3]-LD[,4]),max(LD[,3]+LD[,4]))) # mean -/+  sd

#----------------------------------------------- Plot
#  Reference lines at
# 768 (coronation of Charlemagne),
# 1000 (start of second millennium),
# 1140 (start of Gothic architecture)
# 1315 (beginning of the Great Northern European Famine) and
# 1348 (advent of the Black Death).

# breakpoints at 1140 and 1380

ggplot(LD) +
  aes(decade,
      mean) +
  geom_point() +
  geom_line() +
  geom_vline(xintercept = c(768, 1140, 1315, 1348),
             colour = "red") +
  theme_minimal()

LD1 <- LD[LD$decade < 1140, ]
LD2 <- LD[LD$decade >= 1140 & LD$decade <= 1348, ]
LD3 <- LD[LD$decade >= 1348, ]

#================================================ Rate calc LD1
LD1[1:3,]
n=length(LD1[,2]);n
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
    meandiff=LD1[(i+k),3]-LD1[i,3]				#mean diff.
    poolsd=PoolSD(LD1[i+k,5],LD1[i,5],LD1[i+k,4],LD1[i,4])
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

#================================================ Rate calc LD2
LD2[1:3,]
n=length(LD2[,2]);n
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
    meandiff=LD2[(i+k),3]-LD2[i,3]				#mean diff.
    poolsd=PoolSD(LD2[i+k,5],LD2[i,5],LD2[i+k,4],LD2[i,4])
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

#================================================ Rate calc LD3
LD3[1:3,]
n=length(LD3[,2]);n
nn=.5*(n-1)*n;nn
idr=matrix(nrow=nn,ncol=9)
colnames(idr)=c('int','diff.sd','rate.sd','log.i','log.d','log.r',
                'sbn','wgt','sum')
nc=0
stdev<-numeric(length=n-1)
for (k in 1:(n-1)){    #run length
  for (i in 1:(n-k)){  #starting position
    nc=nc+1
    idr[nc,1]=k						              	#intercept
    meandiff=LD3[(i+k),3]-LD3[i,3]				#mean diff.
    poolsd=PoolSD(LD3[i+k,5],LD3[i,5],LD3[i+k,4],LD3[i,4])
    idr[nc,2]=meandiff/poolsd			    		#diff.sd
    idr[nc,3]=idr[nc,2]/idr[nc,1]	  			#rate.sd
    idr[nc,4]=log10(idr[nc,1])		    		#log.i
    idr[nc,5]=log10(abs(idr[nc,2]))				#log.d
    idr[nc,6]=log10(abs(idr[nc,3]))				#log.r
    if(idr[nc,1]==1){idr[nc,7]=1}else{idr[nc,7]=3}	#sbn
    idr[nc,8]=1/idr[nc,1]			        		#wgt
    idr[nc,9]=idr[nc,4]+idr[nc,6]		  		#sum
    stdev[i]=poolsd
  }
}

vartot=sum(stdev^2,na.rm=TRUE);varcnt=length(table(stdev,useNA='no'))
varmean=vartot/varcnt; avestdev=sqrt(varmean);avestdev
idr[1:3,]
idrx3=idr[!rowSums(!is.finite(idr)),]  #remove rows that have -Inf


#========================================================== Plot LRI panel b
#assign 'bootn' as boot number
bootn<-1000
#assign 'mode' as "medians","all","mixed"
mode<-"all"
#assign circle size for points (1.5 or 2)
psize<-2	#2
#assign 'equation' position as "normal","lower","none" at end of each call
#send (1)idrx matrix, (2)mode(diff/rate), (3)panel placement coordinate x,
#  (4)panel placement coordinate y, (5)bootn, (6)mode, (7)psize, (8)equation

xr<-c(-20,20)
yr<-c(-10,10)		#xrange;yrange for plot axes
plot(xr,yr,					#set up plot
     type='n',
     xaxt='n',
     yaxt='n',
     axes = FALSE,
     ann=FALSE,
     asp = 1) 				#aspect ratio (y/x))

bootresultd=TriPanelBC(idrx1,"r",-18,0,bootn,mode,psize,"normal")
bootresultd=TriPanelBC(idrx2,"r",-2,0,bootn,mode,psize,"normal")
bootresultd=TriPanelBC(idrx3,"r",13,0,bootn,mode,psize,"normal")

# p 109 of ROEV

# Random time series have differences that scale with a slope at or near 0.500 on
# a log difference versus log interval or LDI plot. The corresponding rates scale
# with a slope at or near 0.500 on a log rate versus log interval or LRI plot.

# Stationary time series have differences that scale with a slope at or near 0.000
# on an LDI plot. The corresponding rates scale with a slope at or near 1.000
# on an LRI plot.

# Directional time series have differences that scale with a slope at or near 1.000
# on an LDI plot. The corresponding rates scale with a slope at or near 0.000 on
# an LRI plot.

# directional:    should include 0.000
# random:         should include 0.500 or -0.500
# stationary:     should include 1.000




