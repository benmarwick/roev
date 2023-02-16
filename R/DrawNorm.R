#' An ROEV Function
#'
#' DrawNorm, histogram of step rates
#' @param mean Mean
#' @param sd Standard deviation
#' @param base Base
#' @param hgt Height
#' @keywords normal
#' @import stats graphics grDevices MASS
#' @export
#' @examples
#' # DrawNorm()

DrawNorm <- function(mean, sd, base, hgt) {
  xfit <- seq(mean - 3 * sd, mean + 3 * sd, length = 121)
  yfit <- dnorm(xfit, mean, sd)
  maxy = max(yfit) / hgt
  polygon(xfit,
          yfit / maxy + base,
          border = NA,
          col = gray(9 / 10))
  lines(xfit, yfit / maxy + base, lwd = 1.2, col = 1)
  lines(c(xfit[61], xfit[61]),
        c(0, yfit[61] / maxy) + base,
        lwd = 1,
        col = gray(4 / 10))
}
