// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// check_continuous_points_above_thr
bool check_continuous_points_above_thr(NumericVector v, int i_start, int num, double thr, int n_skip);
RcppExport SEXP _Met4DX_check_continuous_points_above_thr(SEXP vSEXP, SEXP i_startSEXP, SEXP numSEXP, SEXP thrSEXP, SEXP n_skipSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type v(vSEXP);
    Rcpp::traits::input_parameter< int >::type i_start(i_startSEXP);
    Rcpp::traits::input_parameter< int >::type num(numSEXP);
    Rcpp::traits::input_parameter< double >::type thr(thrSEXP);
    Rcpp::traits::input_parameter< int >::type n_skip(n_skipSEXP);
    rcpp_result_gen = Rcpp::wrap(check_continuous_points_above_thr(v, i_start, num, thr, n_skip));
    return rcpp_result_gen;
END_RCPP
}
// get_continuous_points_above_thr_idx
NumericVector get_continuous_points_above_thr_idx(NumericVector v, int i_start, int num, double thr, int n_skip);
RcppExport SEXP _Met4DX_get_continuous_points_above_thr_idx(SEXP vSEXP, SEXP i_startSEXP, SEXP numSEXP, SEXP thrSEXP, SEXP n_skipSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type v(vSEXP);
    Rcpp::traits::input_parameter< int >::type i_start(i_startSEXP);
    Rcpp::traits::input_parameter< int >::type num(numSEXP);
    Rcpp::traits::input_parameter< double >::type thr(thrSEXP);
    Rcpp::traits::input_parameter< int >::type n_skip(n_skipSEXP);
    rcpp_result_gen = Rcpp::wrap(get_continuous_points_above_thr_idx(v, i_start, num, thr, n_skip));
    return rcpp_result_gen;
END_RCPP
}
// find_greater_equal_than
NumericVector find_greater_equal_than(NumericVector data, NumericVector values);
RcppExport SEXP _Met4DX_find_greater_equal_than(SEXP dataSEXP, SEXP valuesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type data(dataSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type values(valuesSEXP);
    rcpp_result_gen = Rcpp::wrap(find_greater_equal_than(data, values));
    return rcpp_result_gen;
END_RCPP
}
// find_dense_min
NumericVector find_dense_min(NumericVector values, int i_max);
RcppExport SEXP _Met4DX_find_dense_min(SEXP valuesSEXP, SEXP i_maxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type values(valuesSEXP);
    Rcpp::traits::input_parameter< int >::type i_max(i_maxSEXP);
    rcpp_result_gen = Rcpp::wrap(find_dense_min(values, i_max));
    return rcpp_result_gen;
END_RCPP
}
// rect_unique
NumericVector rect_unique(NumericMatrix data, NumericVector order, NumericVector diff);
RcppExport SEXP _Met4DX_rect_unique(SEXP dataSEXP, SEXP orderSEXP, SEXP diffSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type data(dataSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type order(orderSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type diff(diffSEXP);
    rcpp_result_gen = Rcpp::wrap(rect_unique(data, order, diff));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_Met4DX_check_continuous_points_above_thr", (DL_FUNC) &_Met4DX_check_continuous_points_above_thr, 5},
    {"_Met4DX_get_continuous_points_above_thr_idx", (DL_FUNC) &_Met4DX_get_continuous_points_above_thr_idx, 5},
    {"_Met4DX_find_greater_equal_than", (DL_FUNC) &_Met4DX_find_greater_equal_than, 2},
    {"_Met4DX_find_dense_min", (DL_FUNC) &_Met4DX_find_dense_min, 2},
    {"_Met4DX_rect_unique", (DL_FUNC) &_Met4DX_rect_unique, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_Met4DX(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
