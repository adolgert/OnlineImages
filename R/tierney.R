
#' Retrieve median and confidence intervals.
#' @param base A Tierney quantile estimator.
#'
#' If every image had a NaN for a pixel, that pixel's quantiles will be all NaN.
#' If the pixel had some NaN values, the quantiles will be adjusted to account
#' for those fewer values.
#'
#' @export
quantiles <- function(base) {
  ci <- dim(base$buf)[1] %/% 2
  edge <- 0.5 * (100 - base$level) / 100
  probs <- c(edge, 0.5, 1 - edge)

  rot_param <- aperm(base$param, c(2, 3, 1))
  median <- rot_param[, , 1]
  N <- rot_param[, , 4]
  bds <- dim(median)
  lower <- array(0, dim = bds)
  upper <- array(0, dim = bds)
  for (j in seq(bds[2])) {
    for (i in seq(bds[1])) {
      m <- N[i, j]
      if (m == 0) {
        lower[i, j] <- NaN
        median[i, j] <- NaN
        upper[i, j] <- NaN

      } else if (m <= 2*ci) {
        qq <- unname(quantile(base$buf[1:m, i, j], probs = probs))
        lower[i, j] <- qq[1]
        median[i, j] <- qq[2]
        upper[i, j] <- qq[3]

      } else {
        # leave the median estimate alone.
        # Use the R quantile function on CI b/c it interpolates to account for
        # an entry that had some NaNs.
        low <- edge * m / ci # Expect this to be something like 0.96.
        lower[i, j] <- unname(quantile(base$buf[1:ci, i, j], probs = low))
        upper[i, j] <- unname(quantile(
          base$buf[(ci + 1):(2 * ci), i, j], probs = (1 - low)
          ))
      }
    }
  }
  list(lower = lower, median = median, upper = upper)
}


#' Create an object to gather quantiles.
#' @param extent The dimensions of the image as c(rows, cols).
#' @param total The total number of images from which to calculate quantiles.
#' @param level The confidence interval (C.I.), e.g. 95 or 90.
#'
#' Once you build the object, it expects to be called `total`
#' times with `quantile_add(obj, image)`, where image has
#' `dims=extent`.
#'
#' This returns a list, where `alpha` is the median quantile,
#' `n` is the number of images seen so far, `total` is the total
#' above, `param` is three-dimensional, with four main components,
#' 1: median, 2: fn from the paper, 3: d0 from the paper and 4: number of
#' non-NaN images seen for this pixel. Then the `buf` element is a stack of
#' images which are the lowest and highest in rank, to determine the C.I.
#' The `level` argument is stored, too.
#'
#' @export
build_tierney <- function(extent, total, level = 95) {
  edge = 0.5 * (100 - level) / 100
  # Add one so we can use quantile algorithm that interpolates the C.I.
  ci = as.integer(round(1000 * edge)) + 1
  list(
    alpha = 0.5,  # median quantile
    level = as.numeric(level), # quantile.
    n = 0L,  # Number of images seen.
    total = as.integer(total),  # Total number off images expected.
    param = array(0, dim = c(4, extent)),
    buf = array(0, dim = c(2 * ci, extent))
  )
}
