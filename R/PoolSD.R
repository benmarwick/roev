#' An ROEV Function
#'
#' Pool standard deviations weighted by n
#'
#' @param n1 something
#' @param n2 something
#' @param sd1 something
#' @param sd2 something
#'
#' @keywords standard deviation
#' @export
#' @examples
#' # PoolSD()

PoolSD <- function(n1, n2, sd1, sd2) {
  if (is.na(sd1)) {
    poolsd = sd2
  }
  if (is.na(sd2)) {
    poolsd = sd1
  }
  if (!is.na(sd1) && !is.na(sd2)) {
    nm1 = n1 - 1
    nm2 = n2 - 1
    v1 = sd1 ^ 2
    v2 = sd2 ^ 2
    poolvar = (nm1 * v1 + nm2 * v2) / (nm1 + nm2)
    poolsd = sqrt(poolvar)
  }
  return(poolsd)
}
