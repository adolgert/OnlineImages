#include<algorithm>
#include<cmath>
#include<Rcpp.h>
using namespace Rcpp;

// [[Rcpp::plugins("cpp11")]]

// This file calculates quantiles for a set of images.
// The main feature of these functions is that they estimate the
// median and calculate confidence intervals without storing
// all of the images.
//
// The median is estimated using stochastic approximation. This algorithm
// is from: Tierney, Luke. 1983. “A Space-Efficient Recursive Procedure
// for Estimating a Quantile of an Unknown Distribution.” SIAM Journal on
// Scientific and Statistical Computing. https://doi.org/10.1137/0904048.
//
// The CI quantiles are calculated exactly by keeping all whose rank is
// outside the bounds of the confidence interval on each side.

//' Add an image to the stack from which to compute quantiles.
//'
//' @param base The object that `build_tierney` made.
//' @param image A two-dimensional array of numeric. Can include NaN values.
//' @return A modified base object, that you pass in the next time.
//'
// [[Rcpp::export]]
List quantile_add(List base, NumericMatrix image) {
  int irow = image.nrow();
  int icol = image.ncol();

  // Holds image slices on either end of confidence intervals.
  // dims = c(2 * ci, i, j)
  NumericVector buf = base["buf"];
  // dims = c(3, i, j). Params are median, fn, d0.
  NumericVector param = base["param"];
  IntegerVector param_dim = param.attr("dim");
  int param_dim_cnt = param_dim.length();
  if (param_dim_cnt != 3) {
    stop("parameter vector should have three dimensions");
  }
  int pc = 4; // Count of parameters, the first dimension of param.
  if (param_dim[0] != pc) {
    stop("parameter vector first dimension should be 4");
  }
  IntegerVector buf_dim = buf.attr("dim");
  if (buf_dim.length() != 3) {
    stop("buffer of images should have three dimensions");
  }
  int dci = buf_dim[0];
  int ci = dci >> 1;  // The number of values to keep at high and low.
  int icnt = buf_dim[1];
  int jcnt = buf_dim[2];
  if (icnt != irow) {
    stop("image rows should equal buffer rows");
  }
  if (jcnt != icol) {
    stop("image cols should equal buffer cols");
  }

  // This is the total images seen, but some pixels may have had NaN values.
  int ntotal = as<int>(base["n"]);
  int total = as<int>(base["total"]);  // Number of images expected by the end.
  double alpha = as<double>(base["alpha"]);

  for (R_xlen_t jidx = 0; jidx < jcnt; ++jidx) { // Fortran order from R.
    for (R_xlen_t iidx = 0; iidx < icnt; ++iidx) {
      R_xlen_t ioffset = iidx + icnt * jidx;  // image offset
      R_xlen_t poffset = pc * ioffset;   // param offset
      R_xlen_t boffset = dci * ioffset; // buf offset.

      int n = param(3 + poffset);
      double x = image(iidx, jidx);
      if (std::isnan(x)) {
        continue;
      }

      if (n < dci) {
        buf[n + boffset] = x;
        n++;
        param[3 + poffset] = n;

        if (n == dci) {
          NumericVector::iterator bufit = buf.begin() + boffset;
          NumericVector::iterator bufend = buf.begin() + boffset + dci;
          std::sort(bufit, bufend);
          // Set $\xi_0$ equal to the $10\alpha$th smallest of the first ten
          // observations. \xi is the median:
          double median = 0.5 * (buf(boffset + ci - 1), buf(boffset + ci));
          // Set $d_0^{-1}$ equal to the interquantile range of the first ten
          // observations (specifically to the difference between the eigth
          // and third smallest observations).
          R_xlen_t low = std::lround(dci * 0.3);
          R_xlen_t high = std::lround(dci * 0.8);
          // The initial estimate $\xi_0$ of $\xi$ and $d_0$, an initial
          // estimate of $f(\xi)^{-1}$, are treated as fixed; in practice they
          // can be obtained from a small preliminary sample.
          double d0 = 1.0 / (buf(boffset + high) - buf(boffset + low));
          // The paper says this can start as 0. 1/d0 might work, too.
          double fn = 0;
          param[0 + poffset] = median;
          param[1 + poffset] = fn;
          param[2 + poffset] = d0;
        } // else not equal to dci, so don't organize now.

      } else if (n < total) {
        double h = std::pow(n + 1, -0.5);
        // Tierney stochastic approximation of median to start.
        double xi = param(0 + poffset); // $\hat{\xi}$, the median estimate.
        double fn = param(1 + poffset); // $\hat{f}_n(\xi)$.
        double d0 = param(2 + poffset); // $d_0$.
        // \hat{f}_n can be calculated recursively by setting $\hat{f}_0(\xi)=0$
        // and using the identity
        // \hat{f}_{n+1}(\xi)=\frac{1}{n+1}(n\hat{f}_n(\xi) +
        //                         I(\hat{\xi}_n, X_{n+1},h_{n+1})/(2 h_{n+1}))
        // I(x,y,z) = 1 if |x-y| <= z, 0 if |x-y| > z.
        double I = (std::fabs(xi <= x) <= h) ? 1 : 0;
        double fn1 = (1 / (n + 1)) * (n * fn + I / (2 * h));
        double dn = std::min(1 / fn, d0 * std::pow(n, alpha));
        // Z(x,y) = 1 if x <= y, 0 if x > y.
        // \hat{\xi}_{n+1} = \hat{\xi}_n -
        //                     \frac{d_n}{n+1}(Z(X_{n+1}, \hat{\xi}_n) - \alpha)
        double Z = (x <= xi) ? 1 : 0;
        double xi1 = xi - (dn / (n + 1)) * (Z - alpha);
        param[0 + poffset] = xi1;
        param[1 + poffset] = fn1;
        param[3 + poffset] = n + 1;

        // Then do the CI at the boundaries.
        // Here we keep every value that is lower (or higher) rank.
        double max_low = buf(ci - 1 + boffset);
        double min_high = buf(ci + boffset);
        if (x < max_low) {
          // Add x to the lower C.I.
          buf[ci - 1 + boffset] = x;
          // Swap max value of lower C.I. to the last rank, but leave rest
          // unsorted.
          NumericVector::iterator bufit = buf.begin() + boffset;
          auto next_max_iter = std::max_element(bufit, bufit + ci);
          // Translate to integer b/c you can't write to Rcpp iterator.
          auto next_max_cnt = std::distance(bufit, next_max_iter);
          double item = buf(boffset + ci - 1);
          buf[boffset + ci - 1] = buf(boffset + next_max_cnt);
          buf[boffset + next_max_cnt] = item;

        } else if (x > min_high) {
          // Add x to the higher C.I.
          buf[ci + boffset] = x;
          // Swap min value of C.I. to the first rank, but leave rest alone.
          NumericVector::iterator bufit = buf.begin() + boffset + ci;
          auto next_min_iter = std::min_element(bufit, bufit + ci);
          auto next_min_cnt = std::distance(bufit, next_min_iter);
          double item = buf(boffset + ci);
          buf[boffset + ci] = buf(boffset + ci + next_min_cnt);
          buf[boffset + ci + next_min_cnt] = item;
        }
        // else x isn't relevant for rank of confidence intervals.
      } // else too many observations.
    }
  }
  ntotal++;

  return List::create(
    Named("alpha") = base["alpha"],
    Named("level") = base["level"],
    Named("n") = ntotal,
    Named("total") = base["total"],
    Named("param") = param,
    Named("buf") = buf
  );
}
