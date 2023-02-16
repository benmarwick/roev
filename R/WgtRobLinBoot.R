#' An ROEV Function
#'
#' Bootstrap weighted robust line fit
#'
#' @param datmat something
#' @param bootn something
#'
#' @keywords WgtRobLinBoot
#' @import MASS
#' @export
#' @examples
#' # WgtRobLinBoot()

WgtRobLinBoot<-function(datmat,bootn){	#input as datmat/medmat nx4 matrix
wrlb<-numeric(length=7)
boot025=ceiling(.025*bootn)
boot975=floor(.975*bootn)
bootmat<-matrix(nrow=bootn,ncol=2);colnames(bootmat)=c("slope","intercept")
bootcount<-0
for (b in 1:bootn){
	bootcount=bootcount+1
	bootr<-sample(length(datmat[,1]),replace=T)#row nos to use in resampling
	simat=datmat[bootr,]
	bootmat[bootcount,1]=WgtRobLinFit(simat)[1]
	bootmat[bootcount,2]=WgtRobLinFit(simat)[2]
}
print(bootcount)
bootmats1=sort(bootmat[,1])	#slopes
bootmats2=sort(bootmat[,2])	#intercepts
wrlb[1]=bootmats1[boot975]
wrlb[2]=median(bootmats1)
wrlb[3]=bootmats1[boot025]
wrlb[4]=bootmats2[boot975]
wrlb[5]=median(bootmats2)
wrlb[6]=bootmats2[boot025]
wrlb[7]=bootcount
return(list(wrlb,bootmat))
}
