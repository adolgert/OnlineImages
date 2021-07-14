# OnlineImages

Calculates median and confidence intervals for a stack of images without storing the whole stack.

This solves a problem for generating data on maps. The input values are raster images of data. The calculation creates a set of draws taking
that data as input. The output can be 1000 images of floating-point values, so it can use a lot of memory.

This repository is an online quantile algorithm, applied to image processing. It assumes the exact situation we need,
the calculation of a median and two quantile values. The central algorithm is Tierney's stochastic approximation
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
