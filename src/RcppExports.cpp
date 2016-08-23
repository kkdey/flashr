// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// threadcheck
int threadcheck(int i);
RcppExport SEXP flashr_threadcheck(SEXP iSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< int >::type i(iSEXP);
    __result = Rcpp::wrap(threadcheck(i));
    return __result;
END_RCPP
}
// vec_to_mat
SEXP vec_to_mat(SEXP x_);
RcppExport SEXP flashr_vec_to_mat(SEXP x_SEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< SEXP >::type x_(x_SEXP);
    __result = Rcpp::wrap(vec_to_mat(x_));
    return __result;
END_RCPP
}
// SquaredVec
arma::vec SquaredVec(arma::vec Efac);
RcppExport SEXP flashr_SquaredVec(SEXP EfacSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< arma::vec >::type Efac(EfacSEXP);
    __result = Rcpp::wrap(SquaredVec(Efac));
    return __result;
END_RCPP
}
// twopi
double twopi();
RcppExport SEXP flashr_twopi() {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    __result = Rcpp::wrap(twopi());
    return __result;
END_RCPP
}
// ProdVec
arma::vec ProdVec(arma::vec E, arma::vec F);
RcppExport SEXP flashr_ProdVec(SEXP ESEXP, SEXP FSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< arma::vec >::type E(ESEXP);
    Rcpp::traits::input_parameter< arma::vec >::type F(FSEXP);
    __result = Rcpp::wrap(ProdVec(E, F));
    return __result;
END_RCPP
}
// sumFac
double sumFac(Rcpp::NumericVector E);
RcppExport SEXP flashr_sumFac(SEXP ESEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type E(ESEXP);
    __result = Rcpp::wrap(sumFac(E));
    return __result;
END_RCPP
}
// sebeta_const
double sebeta_const(double sigmae2, arma::vec E);
RcppExport SEXP flashr_sebeta_const(SEXP sigmae2SEXP, SEXP ESEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< double >::type sigmae2(sigmae2SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type E(ESEXP);
    __result = Rcpp::wrap(sebeta_const(sigmae2, E));
    return __result;
END_RCPP
}
// ashBuild3
Rcpp::List ashBuild3(NumericMatrix Y, NumericVector Ef, NumericMatrix sigma2_mat);
RcppExport SEXP flashr_ashBuild3(SEXP YSEXP, SEXP EfSEXP, SEXP sigma2_matSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericMatrix >::type Y(YSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Ef(EfSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type sigma2_mat(sigma2_matSEXP);
    __result = Rcpp::wrap(ashBuild3(Y, Ef, sigma2_mat));
    return __result;
END_RCPP
}
// ashBuild1
Rcpp::List ashBuild1(Rcpp::NumericMatrix Y, Rcpp::NumericVector Ef, double sigma2);
RcppExport SEXP flashr_ashBuild1(SEXP YSEXP, SEXP EfSEXP, SEXP sigma2SEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type Y(YSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type Ef(EfSEXP);
    Rcpp::traits::input_parameter< double >::type sigma2(sigma2SEXP);
    __result = Rcpp::wrap(ashBuild1(Y, Ef, sigma2));
    return __result;
END_RCPP
}
// ashBuild2
Rcpp::List ashBuild2(Rcpp::NumericMatrix Y, Rcpp::NumericVector Ef, Rcpp::NumericVector sigma2_vec, bool isRow);
RcppExport SEXP flashr_ashBuild2(SEXP YSEXP, SEXP EfSEXP, SEXP sigma2_vecSEXP, SEXP isRowSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type Y(YSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type Ef(EfSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type sigma2_vec(sigma2_vecSEXP);
    Rcpp::traits::input_parameter< bool >::type isRow(isRowSEXP);
    __result = Rcpp::wrap(ashBuild2(Y, Ef, sigma2_vec, isRow));
    return __result;
END_RCPP
}
// SigmaV_Est
arma::mat SigmaV_Est(arma::mat Y, arma::vec El, arma::vec Ef);
RcppExport SEXP flashr_SigmaV_Est(SEXP YSEXP, SEXP ElSEXP, SEXP EfSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< arma::mat >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type El(ElSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type Ef(EfSEXP);
    __result = Rcpp::wrap(SigmaV_Est(Y, El, Ef));
    return __result;
END_RCPP
}
// SigmaV_Est_2
arma::mat SigmaV_Est_2(arma::mat Y, arma::vec El, arma::vec Ef, arma::vec El_listed, arma::vec Ef_listed);
RcppExport SEXP flashr_SigmaV_Est_2(SEXP YSEXP, SEXP ElSEXP, SEXP EfSEXP, SEXP El_listedSEXP, SEXP Ef_listedSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< arma::mat >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type El(ElSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type Ef(EfSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type El_listed(El_listedSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type Ef_listed(Ef_listedSEXP);
    __result = Rcpp::wrap(SigmaV_Est_2(Y, El, Ef, El_listed, Ef_listed));
    return __result;
END_RCPP
}
// calcLik2
double calcLik2(NumericMatrix sigma2v, NumericMatrix sigma2_true, int n_threads);
RcppExport SEXP flashr_calcLik2(SEXP sigma2vSEXP, SEXP sigma2_trueSEXP, SEXP n_threadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericMatrix >::type sigma2v(sigma2vSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type sigma2_true(sigma2_trueSEXP);
    Rcpp::traits::input_parameter< int >::type n_threads(n_threadsSEXP);
    __result = Rcpp::wrap(calcLik2(sigma2v, sigma2_true, n_threads));
    return __result;
END_RCPP
}
// calcLik0
double calcLik0(NumericMatrix sigma2v, double sigma2_true);
RcppExport SEXP flashr_calcLik0(SEXP sigma2vSEXP, SEXP sigma2_trueSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericMatrix >::type sigma2v(sigma2vSEXP);
    Rcpp::traits::input_parameter< double >::type sigma2_true(sigma2_trueSEXP);
    __result = Rcpp::wrap(calcLik0(sigma2v, sigma2_true));
    return __result;
END_RCPP
}
// rcpparma_hello_world
arma::mat rcpparma_hello_world();
RcppExport SEXP flashr_rcpparma_hello_world() {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    __result = Rcpp::wrap(rcpparma_hello_world());
    return __result;
END_RCPP
}
// rcpparma_outerproduct
arma::mat rcpparma_outerproduct(const arma::colvec& x);
RcppExport SEXP flashr_rcpparma_outerproduct(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const arma::colvec& >::type x(xSEXP);
    __result = Rcpp::wrap(rcpparma_outerproduct(x));
    return __result;
END_RCPP
}
// rcpparma_innerproduct
double rcpparma_innerproduct(const arma::colvec& x);
RcppExport SEXP flashr_rcpparma_innerproduct(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const arma::colvec& >::type x(xSEXP);
    __result = Rcpp::wrap(rcpparma_innerproduct(x));
    return __result;
END_RCPP
}
// rcpparma_bothproducts
Rcpp::List rcpparma_bothproducts(const arma::colvec& x);
RcppExport SEXP flashr_rcpparma_bothproducts(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const arma::colvec& >::type x(xSEXP);
    __result = Rcpp::wrap(rcpparma_bothproducts(x));
    return __result;
END_RCPP
}
