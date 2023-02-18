#' An ROEV Function
#'
#'  Prepare matrix with columns for Generation, Time step, Number of samples, Mean, Stdev, Ln Mean, Ln Stdev. Read to input into `TriPanelBC()`
#'
#' @param x data frame where the first four columns are time-mean-sd-n. The column names are not important, but the order of columns is.
#'
#' @keywords DataPrep
#' @export
#' @examples
#' # DataPrep()



DataPrep <- function(x){

  # our input data frame is called x

  n = length(x[,2])        # get the number of rows (measurements)
  nn = 0.5 * (n-1) * n     # get the number of pairwise comparisons

  idr = matrix(nrow = nn,
               ncol = 9)   # create an empty matrix with 9 columns

  # set the column names of our empty matrix
  colnames(idr) = c('int',
                    'diff.sd',
                    'rate.sd',
                    'log.i',
                    'log.d',
                    'log.r',
                    'sbn',
                    'wgt',
                    'sum')

  # populate the empty matrix with rate calculations
  nc = 0
  stdev <- numeric(length = n - 1)

  for (k in 1:(n - 1)){                           # run length
    for (i in 1:(n - k)){                         # starting position
      nc = nc + 1
      idr[nc, 1] = k						                  # intercept
      meandiff = x[(i + k), 2] - x[i, 2]			    # mean diff.
      poolsd = roev::PoolSD(x[i + k, 4],
                            x[i, 4],
                            x[i + k, 3],
                            x[i, 3])
      idr[nc, 2] = meandiff / poolsd			    		# diff.sd
      idr[nc, 3] = idr[nc, 2] / idr[nc, 1]		  	# rate.sd
      idr[nc, 4] = log10(idr[nc, 1])			      	# log.i
      idr[nc, 5] = log10(abs(idr[nc, 2]))		  		# log.d
      idr[nc, 6] = log10(abs(idr[nc, 3]))			  	# log.r
      if(idr[nc, 1] == 1){idr[nc, 7] = 1} else {idr[nc, 7] = 3}	# sbn
      idr[nc, 8] = 1 / idr[nc, 1]				        	# wgt
      idr[nc, 9] = idr[nc, 4] + idr[nc, 6]		  	# sum
      stdev[i] = poolsd
    }
  }

  idrx1 = idr[!rowSums(!is.finite(idr)), ]  # remove rows that have -Inf

  return(idrx1)

  }
