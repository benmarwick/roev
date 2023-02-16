#' An ROEV Function
#'
#' Weighted-robust-linear-model line fit
#'
#' @param datmat matrix of line nos., weights, indep.var, and dep.var
#'
#' @keywords WgtRobLinFit
#' @import MASS
#' @export
#' @examples
#' # WgtRobLinFit()

WgtRobLinFit<-function(datmat){	#input as datmat nx4 matrix
DAT<-data.frame(datmat)
wrlf<-numeric(length=2)
rlmw<-rlm(datmat[,4]~datmat[,3],DAT,datmat[,2],maxit=200)
wrlf[1]=coef(rlmw)[2]	#slope
wrlf[2]=coef(rlmw)[1]	#intercept
return(wrlf)
}
