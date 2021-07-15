# OnlineImages

Calculates median and confidence intervals for a stack of images without storing the whole stack.

This solves a problem for generating data on maps. The input values are raster images of data. The calculation creates a set of draws taking that data as input. The output can be 1000 images of floating-point values, so it can use a lot of memory. A 2048x4096 image, with 1000 draws, is 62.5GB. This would reduce that to 3.3GB.

Maybe you were running out of memory with code that created a stack of images and then computed quantiles on that stack. Most of this code manipulates image data, but that creates copies of that data, which is a problem:
```R
i <- 5
j <- 10
N <- 12
data_stack <- lapply(1:N, function(idx) array(rnorm(i*j), dim = c(i, j, 1)))
together <- do.call(cbind, data_stack)
draws_first <- aperm(together, c(2, 1))
quantile_stack <- lapply(1:(i*j), function(idx) {
  quantile(draws_first[,idx], probs = c(0.025, 0.5, 0.975))
  })
do.call(rbind, quantile_stack)
```
This package stores a few images, in order to compute the confidence intervals, which have worse statistics. However, it uses a stochastic estimator to find the median, and that stores the equivalent of 3 images.
```R
i <- 5
j <- 10
N <- 12
sa_estimator <- build_tierney(c(i, j), N)
for (jidx in 1:j) {
  for (iidx in 1:i) {
    sa_estimator <- quantile_add(sa_estimator, array(rnorm(i*j), dim = c(i, j)))
  }
}
quantiles(sa_estimator)
```

This repository is an online quantile algorithm, applied to image processing. It assumes the exact situation we need, the calculation of a median and two quantile values. The central algorithm is Tierney's stochastic approximation
to the quantile.

Tierney, Luke. 1983. “A Space-Efficient Recursive Procedure for Estimating
a Quantile of an Unknown Distribution.” SIAM Journal on Scientific and
Statistical Computing. https://doi.org/10.1137/0904048.

I looked at other algorithms. The remedian doesn't converge quickly-enough to handle only
1000 draws. Liechty and Lin have a "Single-pass, low-storage arbitrary quantile estimation, for
massive datasets," Statistics and Computing 13: 91-100. 2003. It's fine, but it uses quite a bit
of storage to estimate the histogram of the PDF. Our target here is near 1000 draws, becuase that's
enough to get good error bounds on the CI, themselves.

The implementation is in C++ because it would be tough to do with block operations.
