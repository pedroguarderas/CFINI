// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// cf_diff_solv_euls
List cf_diff_solv_euls(const Eigen::MatrixXd& alpha, const Eigen::VectorXd& u0, const Eigen::VectorXd& u1, const Eigen::VectorXd& u2, const Eigen::VectorXd& t, const Eigen::VectorXd& x, const bool is_initial);
RcppExport SEXP _CFINI_cf_diff_solv_euls(SEXP alphaSEXP, SEXP u0SEXP, SEXP u1SEXP, SEXP u2SEXP, SEXP tSEXP, SEXP xSEXP, SEXP is_initialSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type u0(u0SEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type u1(u1SEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type u2(u2SEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type t(tSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const bool >::type is_initial(is_initialSEXP);
    rcpp_result_gen = Rcpp::wrap(cf_diff_solv_euls(alpha, u0, u1, u2, t, x, is_initial));
    return rcpp_result_gen;
END_RCPP
}
// cf_diff_solv_cns
List cf_diff_solv_cns(const double& theta, const Eigen::MatrixXd& alpha, const Eigen::VectorXd& u0, const Eigen::VectorXd& u1, const Eigen::VectorXd& u2, const Eigen::VectorXd& t, const Eigen::VectorXd& x, const bool is_initial);
RcppExport SEXP _CFINI_cf_diff_solv_cns(SEXP thetaSEXP, SEXP alphaSEXP, SEXP u0SEXP, SEXP u1SEXP, SEXP u2SEXP, SEXP tSEXP, SEXP xSEXP, SEXP is_initialSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type u0(u0SEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type u1(u1SEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type u2(u2SEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type t(tSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const bool >::type is_initial(is_initialSEXP);
    rcpp_result_gen = Rcpp::wrap(cf_diff_solv_cns(theta, alpha, u0, u1, u2, t, x, is_initial));
    return rcpp_result_gen;
END_RCPP
}
// cf_black_scholes_solv_cns
List cf_black_scholes_solv_cns(const double& sigma, const double& rate, const double& theta, const Eigen::VectorXd& u0, const Eigen::VectorXd& u1, const Eigen::VectorXd& u2, const Eigen::VectorXd& t, const Eigen::VectorXd& x);
RcppExport SEXP _CFINI_cf_black_scholes_solv_cns(SEXP sigmaSEXP, SEXP rateSEXP, SEXP thetaSEXP, SEXP u0SEXP, SEXP u1SEXP, SEXP u2SEXP, SEXP tSEXP, SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double& >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< const double& >::type rate(rateSEXP);
    Rcpp::traits::input_parameter< const double& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type u0(u0SEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type u1(u1SEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type u2(u2SEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type t(tSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(cf_black_scholes_solv_cns(sigma, rate, theta, u0, u1, u2, t, x));
    return rcpp_result_gen;
END_RCPP
}
// cf_uniform_grid
Eigen::VectorXd cf_uniform_grid(const double& a, const double& b, const double& n);
RcppExport SEXP _CFINI_cf_uniform_grid(SEXP aSEXP, SEXP bSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double& >::type a(aSEXP);
    Rcpp::traits::input_parameter< const double& >::type b(bSEXP);
    Rcpp::traits::input_parameter< const double& >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(cf_uniform_grid(a, b, n));
    return rcpp_result_gen;
END_RCPP
}
// cf_adapt_grid
Eigen::VectorXd cf_adapt_grid(const double& l, const double& a, const double& b, const double& n, const double& E);
RcppExport SEXP _CFINI_cf_adapt_grid(SEXP lSEXP, SEXP aSEXP, SEXP bSEXP, SEXP nSEXP, SEXP ESEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double& >::type l(lSEXP);
    Rcpp::traits::input_parameter< const double& >::type a(aSEXP);
    Rcpp::traits::input_parameter< const double& >::type b(bSEXP);
    Rcpp::traits::input_parameter< const double& >::type n(nSEXP);
    Rcpp::traits::input_parameter< const double& >::type E(ESEXP);
    rcpp_result_gen = Rcpp::wrap(cf_adapt_grid(l, a, b, n, E));
    return rcpp_result_gen;
END_RCPP
}
// cf_edo_solv_precor
List cf_edo_solv_precor(const Eigen::VectorXd& t, const Eigen::VectorXd& v0, Function f, const int& m, const double& err);
RcppExport SEXP _CFINI_cf_edo_solv_precor(SEXP tSEXP, SEXP v0SEXP, SEXP fSEXP, SEXP mSEXP, SEXP errSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type t(tSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type v0(v0SEXP);
    Rcpp::traits::input_parameter< Function >::type f(fSEXP);
    Rcpp::traits::input_parameter< const int& >::type m(mSEXP);
    Rcpp::traits::input_parameter< const double& >::type err(errSEXP);
    rcpp_result_gen = Rcpp::wrap(cf_edo_solv_precor(t, v0, f, m, err));
    return rcpp_result_gen;
END_RCPP
}
// cf_tri_diag_solv
void cf_tri_diag_solv(Eigen::VectorXd& a, Eigen::VectorXd& b, Eigen::VectorXd& c, Eigen::VectorXd& d);
RcppExport SEXP _CFINI_cf_tri_diag_solv(SEXP aSEXP, SEXP bSEXP, SEXP cSEXP, SEXP dSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::VectorXd& >::type a(aSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd& >::type b(bSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd& >::type c(cSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd& >::type d(dSEXP);
    cf_tri_diag_solv(a, b, c, d);
    return R_NilValue;
END_RCPP
}
// cf_psor_solv
List cf_psor_solv(const Eigen::VectorXd& u0, const Eigen::MatrixXd& A, const Eigen::VectorXd& b, const Eigen::VectorXd& c, const double& w, const int& n, const double& e);
RcppExport SEXP _CFINI_cf_psor_solv(SEXP u0SEXP, SEXP ASEXP, SEXP bSEXP, SEXP cSEXP, SEXP wSEXP, SEXP nSEXP, SEXP eSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type u0(u0SEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type A(ASEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type b(bSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type c(cSEXP);
    Rcpp::traits::input_parameter< const double& >::type w(wSEXP);
    Rcpp::traits::input_parameter< const int& >::type n(nSEXP);
    Rcpp::traits::input_parameter< const double& >::type e(eSEXP);
    rcpp_result_gen = Rcpp::wrap(cf_psor_solv(u0, A, b, c, w, n, e));
    return rcpp_result_gen;
END_RCPP
}
// cf_wiener
Eigen::MatrixXd cf_wiener(const int& d, const Eigen::VectorXd& t);
RcppExport SEXP _CFINI_cf_wiener(SEXP dSEXP, SEXP tSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int& >::type d(dSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type t(tSEXP);
    rcpp_result_gen = Rcpp::wrap(cf_wiener(d, t));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_CFINI_cf_diff_solv_euls", (DL_FUNC) &_CFINI_cf_diff_solv_euls, 7},
    {"_CFINI_cf_diff_solv_cns", (DL_FUNC) &_CFINI_cf_diff_solv_cns, 8},
    {"_CFINI_cf_black_scholes_solv_cns", (DL_FUNC) &_CFINI_cf_black_scholes_solv_cns, 8},
    {"_CFINI_cf_uniform_grid", (DL_FUNC) &_CFINI_cf_uniform_grid, 3},
    {"_CFINI_cf_adapt_grid", (DL_FUNC) &_CFINI_cf_adapt_grid, 5},
    {"_CFINI_cf_edo_solv_precor", (DL_FUNC) &_CFINI_cf_edo_solv_precor, 5},
    {"_CFINI_cf_tri_diag_solv", (DL_FUNC) &_CFINI_cf_tri_diag_solv, 4},
    {"_CFINI_cf_psor_solv", (DL_FUNC) &_CFINI_cf_psor_solv, 7},
    {"_CFINI_cf_wiener", (DL_FUNC) &_CFINI_cf_wiener, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_CFINI(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
