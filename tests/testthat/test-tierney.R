test_that("tierney builds", {
  i <- 10
  j <- 20
  base <- build_tierney(c(i, j), 1000L)

  expect_equal(base$total, 1000L)
  expect_equal(base$alpha, 0.5)
  expect_equal(base$n, 0L)
})


test_that("tierney adds an image", {
  i <- 10
  j <- 20
  base <- build_tierney(c(i, j), 1000L)
  img <- array(rnorm(i*j), dim = c(i, j))
  base2 <- quantile_add(base, img)
  expect_equal(base2$n, 1L)
})


test_that("tierney sorts the first set of images", {
  i <- 10
  j <- 20
  base <- build_tierney(c(i, j), 1000L)
  ci_cnt <- dim(base$buf)[1]

  for (ci_idx in 1:(ci_cnt - 1)) {
    img <- array(rnorm(i*j, mean = 3.7), dim = c(i, j))
    base <- quantile_add(base, img)
  }
  expect_equal(base$n, 49L)
  img <- array(rnorm(i*j, mean = 3.7), dim = c(i, j))
  base <- quantile_add(base, img)
  trials <- list(c(1, 3), c(7, 4), c(i, j))
  for (t in 1:length(trials)) {
    column = base$buf[, trials[[t]][1], trials[[t]][2]]
    colsort = sort(column)
    expect_equal(colsort, column)
  }
})


relerr <- function(a, b, median, sigma) {
  (b - a - median) / sigma
}


test_that("tierney gives final estimates", {
  i <- 10
  j <- 20
  base <- build_tierney(c(i, j), 1000L)
  ci_cnt <- dim(base$buf)[1]

  for (ci_idx in 1:1000L) {
    img <- array(rnorm(i*j, mean = 3.7), dim = c(i, j))
    base <- quantile_add(base, img)
  }
  expect_equal(base$n, 1000L)
  q <- quantiles(base)
  gold_low <- qnorm(0.025, mean = 3.7)
  gold_hi <- qnorm(0.025, mean = 3.7)
  gold_med <- 3.7
  trials <- list(c(1, 3), c(7, 4), c(i, j))
  for (t in 1:length(trials)) {
    ii <- trials[[t]][1]
    jj <- trials[[t]][2]
    low = q$lower[i, j]
    upper = q$upper[i, j]
    med = q$median[i, j]
    expect_lt(relerr(gold_low, low, 3.7, 1), 0.1)
    expect_lt(relerr(gold_hi, upper, 3.7, 1), 0.1)
    expect_lt(relerr(gold_med, med, 3.7, 1), 0.1)
  }
})
