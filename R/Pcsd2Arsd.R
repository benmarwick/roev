#' An ROEV Function
#'
#' Conversion of a standard deviation from a percent to an arithmetic scale
#'
#' @param pc.m something
#' @param pc.sd something
#'
#' @keywords pc.m,pc.sd
#' @export
#' @examples
#' # Pcsd2Arsd()


Pcsd2Arsd <- function(pc.m, pc.sd) {
  arsd.a = (1 / pi) * acos(1 - .02 * (pc.m - pc.sd))
  arsd.b = (1 / pi) * acos(1 - .02 * (pc.m + pc.sd))
  ar.sd = .5 * (arsd.b - arsd.a)
  return(ar.sd)
}
