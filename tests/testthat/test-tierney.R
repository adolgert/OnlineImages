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

  for (z in 1:ci_cnt) {
    img <- array(rnorm(i*j, mean = 3.7), dim = c(i, j))
    base <- quantile_add(base, img)
  }
  column = base$buf[, 1, 1]
  colsort = sort(column)
  expect_equal(colsort, column)
})
