
#' Retrieve median and confidence intervals.
#' @param base A Tierney quantile estimator.
#' @export
quantiles <- function(base) {
  median <- base$param[1, , ]
  ci <- dim(base$buf)[1] %/% 2
  lower <- base$buf[ci, , ]
  upper <- base$buf[ci + 1, , ]
  list(lower = lower, median = median, upper = upper)
}


#' Create an object to gather quantiles.
#' @param extent The dimensions of the image as c(rows, cols).
#' @param total The total number of images from which to calculate quantiles.
#' @param level The confidence interval, e.g. 95 or 90.
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
#'
#' @export
build_tierney <- function(extent, total, level = 95) {
  edge = 0.5 * (100 - level) / 100
  ci = as.integer(round(1000 * edge))
  list(
    alpha = 0.5,  # median quantile
    n = 0L,  # Number of images seen.
    total = as.integer(total),  # Total number off images expected.
    param = array(0, dim = c(4, extent)),
    buf = array(0, dim = c(2 * ci, extent))
  )
}


#' Trial of the package.
#'
#' @export
trial <- function() {
  # Image dimensions
  i <- 10
  j <- 20
  base <- build_tierney(c(i, j), 1000L)
  img <- array(0, dim = c(i, j))
  base2 <- quantile_add(base, img)
}


trial2 <- function() {
  extent <- c(3, 5, 4)
  deep <- array(0, dim = extent);
  outbound <- sample_matrix(deep);
  outbound
}


trial3 <- function() {
  extent <- c(3, 5, 4)
  deep <- array(0, dim = extent);
  outbound <- iter_matrix(deep);
  outbound
}
