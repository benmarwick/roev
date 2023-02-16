#' An ROEV Function
#'
#' Conversion from a percent to an arithmetic scale
#'
#' @param percent something
#'
#' @keywords percent
#' @export
#' @examples
#' # Pc2Ar()

Pc2Ar <- function(percent) {
  arith = (1 / pi) * acos(1 - 2 * (.01 * percent))
  arith
  return(arith)
}
