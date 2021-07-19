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
  expect_equal(base$n, 51L)
  img <- array(rnorm(i*j, mean = 3.7), dim = c(i, j))
  base <- quantile_add(base, img)
  trials <- list(c(1, 3), c(7, 4), c(i, j))
  for (t in 1:length(trials)) {
    column = base$buf[, trials[[t]][1], trials[[t]][2]]
    colsort = sort(column)
    expect_equal(colsort, column)
  }
})


relerr <- function(a, b, sigma) {
  abs((b - a) / sigma)
}


test_that("tierney gives final estimates", {
  i <- 10
  j <- 20
  N <- 1000L
  stdev <- 1.0
  base <- build_tierney(c(i, j), N)
  ci_cnt <- dim(base$buf)[1]
  all <- array(0, dim = c(i, j, N))

  # Each point gets a different mean so that we know the i,j correspond.
  means <- ((1:N) * 19937) %% 20
  for (ci_idx in 1:N) {
    img <- array(rnorm(i*j, mean = means, sd = stdev), dim = c(i, j))
    base <- quantile_add(base, img)
    all[,,ci_idx] <- img
  }
  expect_equal(base$n, N)
  expect_true(all(base$param[4,,] == N))
  q <- quantiles(base)

  reord_all <- aperm(all, c(3, 1, 2))
  for (jr in 1:j) {
    for (ir in 1:i) {
      qq <- unname(quantile(reord_all[,ir, jr], probs = c(0.025, 0.5, 0.975)))

      low <- q$lower[ir, jr]
      expect_lt(relerr(qq[1], low, stdev), 0.1)
      med <- q$median[ir, jr]
      expect_lt(relerr(qq[2], med, stdev), 0.3)
      upper <- q$upper[ir, jr]
      expect_lt(relerr(qq[3], upper, stdev), 0.1)
    }
  }
})


test_that("tierney handles NaNs", {
  i <- 10
  j <- 20
  N <- 1000L
  stdev <- 1.0
  base <- build_tierney(c(i, j), N)
  ci_cnt <- dim(base$buf)[1]
  all <- array(0, dim = c(i, j, N))

  all_nan <- list(c(1, 2), c(7, 3))
  nan_percent <- list(c(1, 3), c(4, 2), c(9, 9), c(2, 7))
  nan_count <- integer(length = length(nan_percent))

  # Each point gets a different mean so that we know the i,j correspond.
  means <- ((1:N) * 19937) %% 20
  for (ci_idx in 1:N) {
    img <- array(rnorm(i*j, mean = means, sd = stdev), dim = c(i, j))
    for (anidx in seq_along(all_nan)) {
      coord <- all_nan[[anidx]]
      img[coord[0], coord[1]] <- NaN
    }
    for (npidx in seq_along(nan_percent)) {
      np <- nan_percent[[npidx]]
      if (rbinom(1, 1, 0.5) == 1L) {
        img[np[0], np[1]] <- NaN
        nan_count[npidx] = nan_count[npidx] + 1L
      }
    }
    base <- quantile_add(base, img)
    all[,,ci_idx] <- img
  }
  expect_equal(base$n, N)
  q <- quantiles(base)

  reord_all <- aperm(all, c(3, 1, 2))
  for (jr in 1:j) {
    for (ir in 1:i) {
      qq <- unname(quantile(reord_all[,ir, jr], probs = c(0.025, 0.5, 0.975)))

      low <- q$lower[ir, jr]
      expect_lt(relerr(qq[1], low, stdev), 0.1)
      med <- q$median[ir, jr]
      expect_lt(relerr(qq[2], med, stdev), 0.3)
      upper <- q$upper[ir, jr]
      expect_lt(relerr(qq[3], upper, stdev), 0.1)
    }
  }
})
