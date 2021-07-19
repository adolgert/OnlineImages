// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// quantile_add
List quantile_add(List base, NumericMatrix image);
RcppExport SEXP _OnlineImages_quantile_add(SEXP baseSEXP, SEXP imageSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type base(baseSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type image(imageSEXP);
    rcpp_result_gen = Rcpp::wrap(quantile_add(base, image));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_OnlineImages_quantile_add", (DL_FUNC) &_OnlineImages_quantile_add, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_OnlineImages(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
