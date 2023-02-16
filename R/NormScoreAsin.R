#' An ROEV Function
#'
#' Normalize score mean and standard deviation for a score range
#' @param rangemin range minimum
#' @param rangemax range maximum
#' @param scoremean score mean
#' @param scoresd score standard deviation
#' @keywords range, mean, standard deviation, arcsin
#' @export
#' @examples
#' # NormScoreAsin()

NormScoreAsin<-function(rangemin,rangemax,scoremean,scoresd){
  #rangemin=1;rangemax=16;scoremean=8.70;scoresd=sqrt(7.50)
  #rangemin=1;rangemax=16;scoremean=3.34;scoresd=sqrt(7.55)
  rangemid=.5*(rangemin+rangemax)	#range midpoint 8.5
  rangediv=rangemax-rangemid		#range divisor 7.5
  newmean=asin((scoremean-rangemid)/rangediv)
  if (((scoremean-scoresd)-rangemid)/rangediv<(-1)){
    newsd.minus=asin(-1)
  }else{
    newsd.minus=asin(((scoremean-scoresd)-rangemid)/rangediv)
  }
  newsd.plus=asin(((scoremean+scoresd)-rangemid)/rangediv)
  newsd=.5*(newsd.plus-newsd.minus)
  #print(c(rangemid,rangediv,newmean,newsd.minus,newsd.plus,newsd))
  normscore=c(newmean,newsd)
  return(normscore)
}
