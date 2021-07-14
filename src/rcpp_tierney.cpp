#include<algorithm>
#include<cmath>
#include<Rcpp.h>
using namespace Rcpp;

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
// The quantiles are calculated exactly by keeping all whose rank is
// outside the bounds of the confidence interval on each side.
//
// [[Rcpp::plugins("cpp11")]]

// [[Rcpp::export]]
List quantile_add(List base, NumericMatrix image) {
  int irow = image.nrow();
  int icol = image.ncol();

  // Holds image slices on either end of confidence intervals.
  // dims = c(2 * ci, i, j)
  NumericVector buf = base["buf"];
  // dims = c(3, i, j). Params are median, fn, d0.
  NumericVector param = base["param"];
  IntegerVector buf_dim = buf.attr("dim");
  int dci = buf_dim[0];
  int ci = dci << 1;  // The number of values to keep at high and low.
  int icnt = buf_dim[1];
  int jcnt = buf_dim[2];
  int pc = 3; // Count of parameters

  int n = as<int>(base["n"]);
  double alpha = as<double>(base["alpha"]);

  if (n < 2 * ci) {
    for (R_xlen_t jidx = 0; jidx < jcnt; ++jidx) {
      for (R_xlen_t iidx = 0; iidx < icnt; ++iidx) {
        R_xlen_t offset = n + dci * (iidx + icnt * jidx);
        buf[offset] = image(iidx, jidx);
      }
    }
    n += 1;
    // Once buffer is filled, use it to initialize stochastic approximation.
    // Also set up the buffer for confidence intervals.
    if (n == dci) {
      for (R_xlen_t jinit = 0; jinit < jcnt; ++jinit) {
        for (R_xlen_t iinit = 0; iinit < icnt; ++iinit) {
          NumericVector::iterator bufit = buf.begin();
          R_xlen_t column = dci * (iinit + icnt * jinit);
          bufit += column;
          NumericVector::iterator bufend = buf.begin();
          bufend += dci + column;
          std::sort(bufit, bufend);  // Does this work with Rcpp pointers?
          double median = 0.5 * (buf(column + ci - 1), buf(column + ci));
          R_xlen_t low = std::lround(dci * 0.3);
          R_xlen_t high = std::lround(dci * 0.8);
          double fn = buf(column + high) - buf(column + low);
          double d0 = 1.0 / fn;
          param[0 + pc * iinit + pc * icnt * jinit] = median;
          param[1 + pc * iinit + pc * icnt * jinit] = fn;
          param[2 + pc * iinit + pc * icnt * jinit] = d0;
        }
      }
    }
  } else {
    double h = std::pow(n, -0.5);
    for (R_xlen_t jadd = 0; jadd < jcnt; ++jadd) {
      for (R_xlen_t iadd = 0; iadd < icnt; ++iadd) {
        R_xlen_t poffset = pc * iadd + pc * icnt * jadd;
        R_xlen_t boffset = dci * iadd + dci * icnt * jadd;

        // Tierney stochastic approximation of median to start.
        double x = image(iadd, jadd);
        double xi = param(0 + poffset);
        double fn = param(1 + poffset);
        double d0 = param(2 + poffset);
        double I;
        if (std::fabs(xi <= x) <= h) {
          I = 1;
        } else {
          I = 0;
        }
        double fn1 = (1 / (n + 1)) * (n * fn + I / (2 * h));
        double dn = std::min(1 / fn, d0 * std::pow(n, alpha));
        double Z;
        if (x <= xi) {
          Z = 1;
        } else{
          Z = 0;
        }
        double xi1 = xi - (dn / (n + 1)) * (Z - alpha);
        param[0 + poffset] = xi1;
        param[1 + poffset] = fn1;

        // Then do the CI at the boundaries.
        // Here we keep every value that is lower (or higher) rank.
        double max_low = param(ci - 1 + dci * iadd +  dci* icnt * jadd);
        double min_high = param(ci + dci * iadd +  dci* icnt * jadd);
        if (x < max_low) {
          NumericVector::iterator bufit = buf.begin() + boffset;
          buf[boffset + ci - 1] = x;
          auto next_max_iter = std::max_element(bufit, bufit + ci);
          // Translate to integer b/c you can't write to Rcpp iterator.
          auto next_max_cnt = std::distance(bufit, next_max_iter);
          double item = buf(boffset + ci - 1);
          buf[boffset + ci - 1] = buf(boffset + next_max_cnt);
          buf[boffset + next_max_cnt] = item;

        } else if (x > min_high) {
          NumericVector::iterator bufit = buf.begin() + boffset;
          buf[boffset + ci] = x;
          auto next_min_iter = std::min_element(bufit + ci, bufit + dci);
          // Translate to integer b/c you can't write to Rcpp iterator.
          auto next_min_cnt = std::distance(bufit + ci, next_min_iter);
          double item = buf(boffset + ci);
          buf[boffset + ci] = buf(boffset + ci + next_min_cnt);
          buf[boffset + ci + next_min_cnt] = item;
        }
      }
    }
    n++;
  }

  return List::create(
    Named("alpha") = base["alpha"],
    Named("n") = n,
    Named("total") = base["total"],
    Named("param") = param,
    Named("buf") = buf
  );
}


// Teach me about writing in the correct order.
// [[Rcpp::export]]
IntegerVector sample_matrix(IntegerVector deep) {
  IntegerVector buf_dim = deep.attr("dim");
  assert(buf_dim.size() == 3);
  R_xlen_t n = buf_dim(0);
  R_xlen_t icnt = buf_dim(1);
  R_xlen_t jcnt = buf_dim(2);
  int cnt = 0;
  for (R_xlen_t jidx = 0; jidx < jcnt; jidx++) {
    for (R_xlen_t iidx = 0; iidx < icnt; iidx++) {
      for (R_xlen_t nidx = 0; nidx < n; nidx++) {
        R_xlen_t offset = nidx + n * iidx + n * icnt * jidx;
        deep[offset] = cnt;
        cnt++;
      }
    }
  }
  return(deep);
}


// Teach me how iterators work.
// [[Rcpp::export]]
IntegerVector iter_matrix(IntegerVector deep) {
  IntegerVector buf_dim = deep.attr("dim");
  assert(buf_dim.size() == 3);
  R_xlen_t n = buf_dim(0);
  R_xlen_t icnt = buf_dim(1);
  R_xlen_t jcnt = buf_dim(2);
  int cnt = 0;

  IntegerVector::iterator diter = deep.begin();
  for (R_xlen_t jidx = 0; jidx < jcnt; jidx++) {
    for (R_xlen_t iidx = 0; iidx < icnt; iidx++) {
      IntegerVector::iterator coliter = diter;
      for (R_xlen_t nidx = 0; nidx < n; nidx++) {
        *diter = cnt;
        diter++;
      }
//diter += n;
    }
  }
  return(deep);
}
