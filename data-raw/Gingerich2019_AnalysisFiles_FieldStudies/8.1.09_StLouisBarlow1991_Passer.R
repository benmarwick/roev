##Philip D. Gingerich RATES OF EVOLUTION (Cambridge University Press 2019) 
##8.1.09_StLouisBarlow1991_Passer
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
#============================================================ Load file
file1<-"C://R_aaROEVchapt08//8.1.09_StLouisBarlow1991_Passer.csv"
P<-read.csv(file1,header=TRUE,sep=",",quote="\"",dec=".",fill=TRUE)
P
#----------------------------------------------------------------------
P[1,]
LP<-matrix(nrow=nrow(P),ncol=12)		#Ln matrix for Passer montanus
  colnames(LP)=c('Ngf','Mgf','Sgf','Ngm','Mgm','Sgm','Naf','Maf','Saf',
    'Nam','Mam','Sam')
for (i in 1:nrow(P)){
  LP[i,1]=P[i,2]	#-----------------N
  LP[i,2]=log(P[i,3])			#ln mean
  sdgf=sqrt(P[i,2])*P[i,4]
  LP[i,3]=sdgf/P[i,3]			#ln sd
  LP[i,4]=P[i,5]	#-----------------N
  LP[i,5]=log(P[i,6])			#ln mean
  sdgm=sqrt(P[i,5])*P[i,7]
  LP[i,6]=sdgm/P[i,6]			#ln sd
  LP[i,7]=P[i,8]	#-----------------N
  LP[i,8]=log(P[i,9])			#ln mean
  sdaf=sqrt(P[i,8])*P[i,10]
  LP[i,9]=sdaf/P[i,9]			#ln sd
  LP[i,10]=P[i,11]#-----------------N
  LP[i,11]=log(P[i,12])			#ln mean
  sdam=sqrt(P[i,11])*P[i,13]
  LP[i,12]=sdam/P[i,12]			#ln sd
}
LP[1,]
#============================= Rate calc 
gentime=2
idr=matrix(nrow=2*nrow(P),ncol=9)
  colnames(idr)=c('int','diff.sd','rate.sd.g','log.i','log|d|','log|r|',
    'sbn','wgt','fgen')	#fgen is fraction of a generation
nc=0
dt=114
for (i in 1:nrow(P)){							#females	
      nc=nc+1
      idr[nc,1]=dt/gentime						#interval gen
      meandiff=LP[i,8]-LP[i,2]					#mean diff.
      poolsd=PoolSD(LP[i,1],LP[i,7],
        LP[i,3],LP[i,9])						#n1,n2,sd1,sd2
      idr[nc,2]=meandiff/poolsd					#diff.sd
      idr[nc,3]=abs(idr[nc,2])/idr[nc,1]				#rate.sd.gen
      idr[nc,4]=log10(idr[nc,1])					#log.i
      idr[nc,5]=log10(abs(idr[nc,2]))				#log.d
      idr[nc,6]=log10(idr[nc,3])					#log.r
      idr[nc,7]=2								#sbn
      idr[nc,8]=1/idr[nc,1]						#wgt
      idr[nc,9]=idr[nc,4]+idr[nc,6]
    }
  }
}
for (i in 1:nrow(P)){							#males	
      nc=nc+1
      idr[nc,1]=dt/gentime						#interval gen
      meandiff=LP[i,11]-LP[i,5]					#mean diff.
      poolsd=PoolSD(LP[i,4],LP[i,10],
        LP[i,6],LP[i,12])						#n1,n2,sd1,sd2
      idr[nc,2]=meandiff/poolsd					#diff.sd
      idr[nc,3]=abs(idr[nc,2])/idr[nc,1]				#rate.sd.gen
      idr[nc,4]=log10(idr[nc,1])					#log.i
      idr[nc,5]=log10(abs(idr[nc,2]))				#log.d
      idr[nc,6]=log10(idr[nc,3])					#log.r
      idr[nc,7]=2								#sbn
      idr[nc,8]=1/idr[nc,1]						#wgt
      idr[nc,9]=idr[nc,4]+idr[nc,6]
    }
  }
}
idr
idrx=idr[!rowSums(!is.finite(idr)),]  #remove rows that have -Inf
nrow(idr);nrow(idrx)
plot(idr[,4],idr[,6])
writefile<-"C://R_aaROEVchapt08//8.1.09_StLouisBarlow1991_Passer_out.csv"
write.csv(idrx,file=writefile,na="NA",row.names=FALSE) 


