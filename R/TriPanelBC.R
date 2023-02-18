#' An ROEV Function
#'
#' TriPanel B or C, drawing and analyzing LDI or LRI temporal scaling
#'
#' @param dr something
#' @param xac something
#' @param yac something
#' @param bootn something
#' @param mode something
#' @param psize something
#' @param equation something
#' @param idrx matrix with -Inf removed, n, diff or rate, xac, yac
#'
#' @keywords LDI LRI
#' @export
#' @examples
#' # TriPanelBC()


TriPanelBC <- function(idrx,
                       dr,
                       xac,
                       yac,
                       bootn,
                       mode,
                       psize,
                       equation) {
  #idrx is int. diff. rate matrix with -Inf rows removed
  # sometimes a little jitter may help, e.g. apply(idrx, 2, jitter, amount = 0.001)

  for (i in 1:ncol(idrx)) {
    idrx[, i] = rev(idrx[, i])
  }
  	#idrx is int. diff. rate matrix now reversed top to bottom
  #-------------------------------------------------------------------- y-axis
  lines(
    c(xac - 1, xac - 1),
    c(yac - 5, yac + 3),
    lty = 1,
    lwd = 1,
    col = 1
  )
  lines(
    c(xac, xac),
    c(yac - 5, yac + 3),
    lty = 1,
    lwd = 1,
    col = gray(6 / 10)
  )
  for (i in 0:7) {
    lines(
      c(xac - 1.15, xac - 1),
      c(yac - 4 + i, yac - 4 + i),
      lty = 1,
      lwd = 1,
      col = 1
    )
  }
  if (dr == "d") {
    text(
      xac - 2.2,
      yac - 1,
      expression(Log[10] ~ abs( ~ difference ~ italic(d) * ' ')),
      srt = 90,
      cex = 1.3,
      col = 1
    )
  } else{
    text(
      xac - 2.2,
      yac - 1,
      expression(Log[10] ~ abs( ~ rate ~ italic(r) * ' ')),
      srt = 90,
      cex = 1.3,
      col = 1
    )
  }
  for (i in 0:3) {
    text(xac - 1.5, yac - 4 + 2 * i, i - 2, cex = 1, col = 1)
  }
  #-------------------------------------------------------------------- x-axis
  lines(
    c(xac - 1, xac + 8),
    c(yac - 5, yac - 5),
    lty = 1,
    lwd = 1,
    col = 1
  )
  for (i in 0:7) {
    lines(
      c(xac + i, xac + i),
      c(yac - 5.15, yac - 5),
      lty = 1,
      lwd = 1,
      col = 1
    )
  }
  for (i in 0:3) {
    text(xac + 2 * i, yac - 5.5, i, cex = 1, col = 1)
  }
  text(
    xac + 3.5,
    yac - 6.2,
    expression(Log[10] ~ interval ~ italic(i)),
    cex = 1.3,
    col = 1
  )
  #-------------------------------- slope, intercept, and confidence intervals
  #TS.a1a0=TheilSen(idrx[,c(4,5)])	#call TheilSen function
  #---------------------------------------------------Full matrix
  if (dr == "d") {
    datmat = idrx[, c(1, 8, 4, 5)]
  }
  if (dr == "r") {
    datmat = idrx[, c(1, 8, 4, 6)]
  }
  #---------------------------------------------------Median matrix
  n = length(idrx[idrx[, 1] == 1, 1])
  medvvec <- numeric(length = n)	#median vector of dependent variable
  if (dr == "d") {
    for (i in 1:n) {
      tempstat = summary(idrx[idrx[, 1] == i, 5])
      medvvec[i] = tempstat[3]
    }
  }
  if (dr == "r") {
    for (i in 1:n) {
      tempstat = summary(idrx[idrx[, 1] == i, 6])
      medvvec[i] = tempstat[3]
    }
  }
  #--------------------------------------------------
  medmat <- matrix(nrow = n, ncol = 4)
  colnames(medmat) = c("int", "wt", "medint", "medval")
  medmat[, 1] = c(1:n)			#interval
  medmat[, 2] = 1 / medmat[, 1]		#weights
  medmat[, 3] = log10(medmat[, 1])		#intervals
  medmat[, 4] = medvvec[]		#values
  MED <- data.frame(medmat)
  #print(MED)
  #---------------------------------------------------Robust linear modeling
  if (mode == "medians") {
    si.cal = WgtRobLinFit(medmat)	#slope-intercept rlm
    output = WgtRobLinBoot(medmat, bootn)
    wrlb = output[[1]]
    bootmat = output[[2]]	#unpack
    #wrlb[2]=si.cal[1];wrlb[5]=si.cal[2]#replace with calculated medians
  }
  if (mode == "all") {
    si.cal = WgtRobLinFit(datmat)	#slope-intercept rlm
    output = WgtRobLinBoot(datmat, bootn)
    wrlb = output[[1]]
    bootmat = output[[2]]	#unpack
    wrlb[2] = si.cal[1]
    wrlb[5] = si.cal[2]#replace with calculated medians
  }
  if (mode == "mixed") {
    si.cal = WgtRobLinFit(datmat)	#slope-intercept rlm
    output = WgtRobLinBoot(medmat, bootn)
    wrlb = output[[1]]
    bootmat = output[[2]]	#unpack
    wrlb[2] = si.cal[1]
    wrlb[5] = si.cal[2]#replace with calculated medians
  }
  # print(wrlb)
  #-------------------------------------------------------------------------
  #med.li=median(idrx[,4])	#median of intervals
  #med.ld=median(idrx[,5])	#median of differences
  #rcnt.r=length(idrx[,5])
  #intercept=med.ld-wrlb[2]*med.li
  dx1 = min(idrx[, 4]) - .25
  dx2 = max(idrx[, 4]) + .25
  dy1.1 = wrlb[1] * dx1 + wrlb[6]
  dy2.1 = wrlb[1] * dx2 + wrlb[6]
  dy1.2 = wrlb[2] * dx1 + wrlb[5]
  dy2.2 = wrlb[2] * dx2 + wrlb[5]
  dy1.3 = wrlb[3] * dx1 + wrlb[4]
  dy2.3 = wrlb[3] * dx2 + wrlb[4]
  lines(
    xac + 2 * c(dx1, dx2),
    yac + 2 * c(dy1.1, dy2.1),
    lty = 2,
    lwd = 1,
    col = 1
  )
  lines(
    xac + 2 * c(dx1, dx2),
    yac + 2 * c(dy1.2, dy2.2),
    lty = 1,
    lwd = 3,
    col = 1
  )
  lines(
    xac + 2 * c(dx1, dx2),
    yac + 2 * c(dy1.3, dy2.3),
    lty = 2,
    lwd = 1,
    col = 1
  )
  #-------------------------------------------------------- points on panel b
  if (dr == "d") {
    if (mode == "all") {
      for (i in 1:length(idrx[, 5])) {
        if (idrx[i, 5] > -2.35) {
          if (idrx[i, 1] == 1) {
            points(
              xac + 2 * idrx[i, 4],
              yac + 2 * idrx[i, 5],
              pch = 21,
              cex = psize,
              bg = gray(6 / 10),
              col = 1
            )
          }
          if (idrx[i, 1] > 1 && idrx[i, 1] <= 2) {
            points(
              xac + 2 * idrx[i, 4],
              yac + 2 * idrx[i, 5],
              pch = 21,
              cex = psize,
              bg = gray(8 / 10),
              col = 1
            )
          }
          if (idrx[i, 1] > 2 && idrx[i, 1] <= 3) {
            points(
              xac + 2 * idrx[i, 4],
              yac + 2 * idrx[i, 5],
              pch = 21,
              cex = psize,
              bg = gray(9 / 10),
              col = 1
            )
          }
          if (idrx[i, 1] > 3) {
            points(
              xac + 2 * idrx[i, 4],
              yac + 2 * idrx[i, 5],
              pch = 21,
              cex = psize,
              bg = "white",
              col = gray(4 / 10)
            )
          }
        }
      }
    } else{
      for (i in 1:(n - 1)) {
        if (medmat[i, 4] > -2.35) {
          if (medmat[i, 1] == 1) {
            points(
              xac + 2 * medmat[i, 3],
              yac + 2 * medmat[i, 4],
              pch = 21,
              cex = psize,
              bg = gray(6 / 10),
              col = 1
            )
          }
          if (medmat[i, 1] > 1 && idrx[i, 1] <= 2) {
            points(
              xac + 2 * medmat[i, 3],
              yac + 2 * medmat[i, 4],
              pch = 21,
              cex = psize,
              bg = gray(8 / 10),
              col = 1
            )
          }
          if (medmat[i, 1] > 2 && idrx[i, 1] <= 3) {
            points(
              xac + 2 * medmat[i, 3],
              yac + 2 * medmat[i, 4],
              pch = 21,
              cex = psize,
              bg = gray(9 / 10),
              col = 1
            )
          }
          if (medmat[i, 1] > 3) {
            points(
              xac + 2 * medmat[i, 3],
              yac + 2 * medmat[i, 4],
              pch = 21,
              cex = psize,
              bg = "white",
              col = gray(4 / 10)
            )
          }
        }
      }
    }
  }
  if (wrlb[4] >= -2.5) {
    points(
      xac + 2 * 0,
      yac + 2 * wrlb[4],
      pch = 25,
      cex = 1,
      bg = 'white',
      col = 1
    )
  }
  if (wrlb[6] <= 1.5) {
    points(
      xac + 2 * 0,
      yac + 2 * wrlb[6],
      pch = 24,
      cex = 1,
      bg = 'white',
      col = 1
    )
  }
  #-------------------------------------------------------- points on panel c
  if (dr == "r") {
    if (mode == "all") {
      for (i in 1:length(idrx[, 6])) {
        if (idrx[i, 6] > -2.35) {
          if (idrx[i, 1] == 1) {
            points(
              xac + 2 * idrx[i, 4],
              yac + 2 * idrx[i, 6],
              pch = 21,
              cex = psize,
              bg = gray(6 / 10),
              col = 1
            )
          }
          if (idrx[i, 1] > 1 && idrx[i, 1] <= 2) {
            points(
              xac + 2 * idrx[i, 4],
              yac + 2 * idrx[i, 6],
              pch = 21,
              cex = psize,
              bg = gray(8 / 10),
              col = 1
            )
          }
          if (idrx[i, 1] > 2 && idrx[i, 1] <= 3) {
            points(
              xac + 2 * idrx[i, 4],
              yac + 2 * idrx[i, 6],
              pch = 21,
              cex = psize,
              bg = gray(9 / 10),
              col = 1
            )
          }
          if (idrx[i, 1] > 3) {
            points(
              xac + 2 * idrx[i, 4],
              yac + 2 * idrx[i, 6],
              pch = 21,
              cex = psize,
              bg = "white",
              col = gray(4 / 10)
            )
          }
        }
      }
    } else{
      for (i in 1:(n - 1)) {
        if (medmat[i, 4] > -2.35) {
          if (medmat[i, 1] == 1) {
            points(
              xac + 2 * medmat[i, 3],
              yac + 2 * medmat[i, 4],
              pch = 21,
              cex = psize,
              bg = gray(6 / 10),
              col = 1
            )
          }
          if (medmat[i, 1] > 1 && idrx[i, 1] <= 2) {
            points(
              xac + 2 * medmat[i, 3],
              yac + 2 * medmat[i, 4],
              pch = 21,
              cex = psize,
              bg = gray(8 / 10),
              col = 1
            )
          }
          if (medmat[i, 1] > 2 && idrx[i, 1] <= 3) {
            points(
              xac + 2 * medmat[i, 3],
              yac + 2 * medmat[i, 4],
              pch = 21,
              cex = psize,
              bg = gray(9 / 10),
              col = 1
            )
          }
          if (medmat[i, 1] > 3) {
            points(
              xac + 2 * medmat[i, 3],
              yac + 2 * medmat[i, 4],
              pch = 21,
              cex = psize,
              bg = "white",
              col = gray(4 / 10)
            )
          }
        }
      }
    }
  }
  if (wrlb[4] >= -2.5) {
    points(
      xac + 2 * 0,
      yac + 2 * wrlb[4],
      pch = 25,
      cex = 1,
      bg = 'white',
      col = 1
    )
  }
  if (wrlb[6] <= 1.5) {
    points(
      xac + 2 * 0,
      yac + 2 * wrlb[6],
      pch = 24,
      cex = 1,
      bg = 'white',
      col = 1
    )
  }
  #------------------------------------------------------------------------
  if (equation == "normal") {
    rect(
      xac + 5 - 2.8,
      yac + 3 - .8,
      xac + 5 + 2.8,
      yac + 3 - .2,
      col = "white",
      border = NA
    )
    if (wrlb[5] < 0) {
      text(
        xac + 5,
        yac + 3,
        bquote(Y == .(format(
          wrlb[2], digits = 1, nsmall = 3
        ))
        ~ X ~ - ~ .(format(
          abs(wrlb[5]), digits = 1, nsmall = 3
        ))),
        pos = 1,
        cex = 1.2,
        col = 1
      )
    } else{
      text(
        xac + 5,
        yac + 3,
        bquote(Y == .(format(
          wrlb[2], digits = 1, nsmall = 3
        ))
        ~ X ~ + ~ .(format(
          abs(wrlb[5]), digits = 1, nsmall = 3
        ))),
        pos = 1,
        cex = 1.2,
        col = 1
      )
    }
    if (mode == "all") {
      rect(
        xac + 5 - 1.9,
        yac + 2.1 - .7,
        xac + 4 + 2.6,
        yac + 2.1 - .2,
        col = "white",
        border = NA
      )
      text(
        xac + 5,
        yac + 2.1,
        bquote(#label n,r2,s
          italic(N) == .(
            format(length(idrx[, 1]), digits = 2, nsmall = 2)
          )),
        pos = 1,
        cex = 1,
        col = 1
      )
    }
    if (mode == "medians") {
      rect(
        xac + 5 - 1.4,
        yac + 2.1 - .7,
        xac + 4 + 2.1,
        yac + 2.1 - .2,
        col = "white",
        border = NA
      )
      text(
        xac + 5,
        yac + 2.1,
        bquote(#label n,r2,s
          italic(N) == .(format(
            n, digits = 0, nsmall = 0
          ))),
        pos = 1,
        cex = 1,
        col = 1
      )
    }
  }
  if (equation == "lower") {
    rect(xac + 5 - 2.8,
         yac - .3,
         xac + 5 + 2.8,
         yac + .3,
         col = "white",
         border = NA)
    if (wrlb[5] < 0) {
      text(
        xac + 5,
        yac + .5,
        bquote(Y == .(format(
          wrlb[2], digits = 1, nsmall = 3
        ))
        ~ X ~ - ~ .(format(
          abs(wrlb[5]), digits = 1, nsmall = 3
        ))),
        pos = 1,
        cex = 1.2,
        col = 1
      )
    } else{
      text(
        xac + 5,
        yac + .5,
        bquote(Y == .(format(
          wrlb[2], digits = 1, nsmall = 3
        ))
        ~ X ~ + ~ .(format(
          abs(wrlb[5]), digits = 1, nsmall = 3
        ))),
        pos = 1,
        cex = 1.2,
        col = 1
      )
    }
    if (mode == "all") {
      rect(
        xac + 5 - 1.6,
        yac - 1.1,
        xac + 5 + 1.6,
        yac - .6,
        col = "white",
        border = NA
      )
      text(
        xac + 5,
        yac - .4,
        bquote(#label n,r2,s
          italic(N) == .(
            format(length(idrx[, 1]), digits = 2, nsmall = 2)
          )),
        pos = 1,
        cex = 1,
        col = 1
      )
    }
    if (mode == "medians") {
      rect(
        xac + 5 - 1.4,
        yac - 1.1,
        xac + 4 + 2.1,
        yac - .6,
        col = "white",
        border = NA
      )
      text(
        xac + 5,
        yac - .4,
        bquote(#label n,r2,s
          italic(N) == .(format(
            n, digits = 0, nsmall = 0
          ))),
        pos = 1,
        cex = 1,
        col = 1
      )
    }
  }
  rect(xac + 3.7,
       yac - 4.8,
       xac + 8.55,
       yac - 1.7,
       col = 'white',
       border = NA)
  text(
    xac + 3.45,
    yac - 2,
    "Slope max.:",
    pos = 4,
    cex = .9,
    col = 1
  )
  text(
    xac + 8.4,
    yac - 2,
    format(round(wrlb[1], 3), nsmall = 3),
    pos = 2,
    cex = .9,
    col = 1
  )
  text(
    xac + 3.45,
    yac - 2.5,
    "Slope med.:",
    pos = 4,
    cex = .9,
    col = 1
  )
  text(
    xac + 8.4,
    yac - 2.5,
    format(round(wrlb[2], 3), nsmall = 3),
    pos = 2,
    cex = .9,
    col = 1
  )
  text(
    xac + 3.45,
    yac - 3,
    "Slope min.:",
    pos = 4,
    cex = .9,
    col = 1
  )
  text(
    xac + 8.4,
    yac - 3,
    format(round(wrlb[3], 3), nsmall = 3),
    pos = 2,
    cex = .9,
    col = 1
  )
  text(
    xac + 3.45,
    yac - 3.5,
    "Intcpt. max.:",
    pos = 4,
    cex = .9,
    col = 1
  )
  text(
    xac + 8.4,
    yac - 3.5,
    format(round(wrlb[4], 3), nsmall = 3),
    pos = 2,
    cex = .9,
    col = 1
  )
  text(
    xac + 3.45,
    yac - 4,
    "Intcpt. med.:",
    pos = 4,
    cex = .9,
    col = 1
  )
  text(
    xac + 8.4,
    yac - 4,
    format(round(wrlb[5], 3), nsmall = 3),
    pos = 2,
    cex = .9,
    col = 1
  )
  text(
    xac + 3.45,
    yac - 4.5,
    "Intcpt. min.:",
    pos = 4,
    cex = .9,
    col = 1
  )
  text(
    xac + 8.4,
    yac - 4.5,
    format(round(wrlb[6], 3), nsmall = 3),
    pos = 2,
    cex = .9,
    col = 1
  )
  #-------------------------------------------------------- significance
  if (dr == "d") {
    #slope maximum
    if (wrlb[1] < 1 && wrlb[1] > .5) {
      text(xac + 8.4 - .2, yac - 2 + .2, '*')
    }
    if (wrlb[1] < 1 && wrlb[1] < .5) {
      text(xac + 8.4 - .1, yac - 2 + .2, '**')
    }
    #slope minimum
    if (wrlb[3] > 0 && wrlb[3] < .5) {
      text(xac + 8.4 - .2, yac - 3 + .2, '*')
    }
    if (wrlb[3] > 0 && wrlb[3] > .5) {
      text(xac + 8.4 - .1, yac - 3 + .2, '**')
    }
  }
  if (dr == "r") {
    #slope maximum
    if (wrlb[1] < 0 && wrlb[1] > (-.5)) {
      text(xac + 8.4 - .2, yac - 2 + .2, '*')
    }
    if (wrlb[1] < 0 && wrlb[1] < (-.5)) {
      text(xac + 8.4 - .1, yac - 2 + .2, '**')
    }
    #slope minimum
    if (wrlb[3] > -1 && wrlb[3] < (-.5)) {
      text(xac + 8.4 - .2, yac - 3 + .2, '*')
    }
    if (wrlb[3] > -1 && wrlb[3] > (-.5)) {
      text(xac + 8.4 - .1, yac - 3 + .2, '**')
    }
  }
  return(output)	#includes wrlb with bootcount plus unsorted bootmat
}
