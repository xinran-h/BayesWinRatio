// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// update_theta
List update_theta(int N_iter, NumericMatrix dd, int n, arma::vec m0, arma::mat L0, arma::mat S0, double v0, double time_max);
RcppExport SEXP _BayesWinRatio_update_theta(SEXP N_iterSEXP, SEXP ddSEXP, SEXP nSEXP, SEXP m0SEXP, SEXP L0SEXP, SEXP S0SEXP, SEXP v0SEXP, SEXP time_maxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type N_iter(N_iterSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type dd(ddSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type m0(m0SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type L0(L0SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type S0(S0SEXP);
    Rcpp::traits::input_parameter< double >::type v0(v0SEXP);
    Rcpp::traits::input_parameter< double >::type time_max(time_maxSEXP);
    rcpp_result_gen = Rcpp::wrap(update_theta(N_iter, dd, n, m0, L0, S0, v0, time_max));
    return rcpp_result_gen;
END_RCPP
}
// update_theta_univariate
List update_theta_univariate(int N_iter, NumericMatrix dd, int n, double L0, double m0, double v0, double S0, double time_max);
RcppExport SEXP _BayesWinRatio_update_theta_univariate(SEXP N_iterSEXP, SEXP ddSEXP, SEXP nSEXP, SEXP L0SEXP, SEXP m0SEXP, SEXP v0SEXP, SEXP S0SEXP, SEXP time_maxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type N_iter(N_iterSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type dd(ddSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< double >::type L0(L0SEXP);
    Rcpp::traits::input_parameter< double >::type m0(m0SEXP);
    Rcpp::traits::input_parameter< double >::type v0(v0SEXP);
    Rcpp::traits::input_parameter< double >::type S0(S0SEXP);
    Rcpp::traits::input_parameter< double >::type time_max(time_maxSEXP);
    rcpp_result_gen = Rcpp::wrap(update_theta_univariate(N_iter, dd, n, L0, m0, v0, S0, time_max));
    return rcpp_result_gen;
END_RCPP
}
// compare
List compare(int M_iter, int n_current_ctrl, int n_current_trt, const arma::cube& postData, arma::field<arma::vec> time_trt, arma::field<arma::vec> time_ctrl, arma::field<arma::vec> surv_trt, arma::field<arma::vec> surv_ctrl);
RcppExport SEXP _BayesWinRatio_compare(SEXP M_iterSEXP, SEXP n_current_ctrlSEXP, SEXP n_current_trtSEXP, SEXP postDataSEXP, SEXP time_trtSEXP, SEXP time_ctrlSEXP, SEXP surv_trtSEXP, SEXP surv_ctrlSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type M_iter(M_iterSEXP);
    Rcpp::traits::input_parameter< int >::type n_current_ctrl(n_current_ctrlSEXP);
    Rcpp::traits::input_parameter< int >::type n_current_trt(n_current_trtSEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type postData(postDataSEXP);
    Rcpp::traits::input_parameter< arma::field<arma::vec> >::type time_trt(time_trtSEXP);
    Rcpp::traits::input_parameter< arma::field<arma::vec> >::type time_ctrl(time_ctrlSEXP);
    Rcpp::traits::input_parameter< arma::field<arma::vec> >::type surv_trt(surv_trtSEXP);
    Rcpp::traits::input_parameter< arma::field<arma::vec> >::type surv_ctrl(surv_ctrlSEXP);
    rcpp_result_gen = Rcpp::wrap(compare(M_iter, n_current_ctrl, n_current_trt, postData, time_trt, time_ctrl, surv_trt, surv_ctrl));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_BayesWinRatio_update_theta", (DL_FUNC) &_BayesWinRatio_update_theta, 8},
    {"_BayesWinRatio_update_theta_univariate", (DL_FUNC) &_BayesWinRatio_update_theta_univariate, 8},
    {"_BayesWinRatio_compare", (DL_FUNC) &_BayesWinRatio_compare, 8},
    {NULL, NULL, 0}
};

RcppExport void R_init_BayesWinRatio(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
