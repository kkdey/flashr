#include <RcppArmadillo.h>
#include <cmath>
#include <thread>
#include <typeinfo>
#include <omp.h>
#include <chrono>
#include <vector>
#include <stdlib.h>

// [[Rcpp::plugins(openmp)]]

// For Rcpp+OpenMP parallelization, check: https://wbnicholson.wordpress.com/2014/07/10/parallelization-in-rcpp-via-openmp/
// For Armadillo Support, check: http://arma.sourceforge.net/docs.html

// LDFLAGS:  -L/usr/local/opt/openblas/lib
// CPPFLAGS: -I/usr/local/opt/openblas/include

using namespace Rcpp;
using namespace arma;
using namespace std;

int f(int i) {
  std::this_thread::sleep_for (std::chrono::seconds(1));
  return i;
}

// [[Rcpp::depends(RcppArmadillo)]]


// [[Rcpp::export]]
int threadcheck(int i)
{
  int n_threads = i;
  std::cout << n_threads << " threads ..." << std::endl;
  std::vector<int> M(12);
#pragma omp parallel for num_threads(n_threads)
  for (int i=0; i<12; i++)
    M[i] = f(i);
  return 0;
}

// [[Rcpp::export]]
SEXP vec_to_mat(SEXP x_) {
  std::vector<double> x = as< std::vector<double> >(x_);
  /* ... do something with x ... */
  NumericVector output = wrap(x);
  output.attr("dim") = Dimension(x.size(), 1);
  return output;
}

// [[Rcpp::export]]
arma::vec SquaredVec(arma::vec Efac) {
  arma::vec Efac2 = Efac % Efac;
  return Efac2;
}

// [[Rcpp::export]]
double twopi() { return M_PI; }

// [[Rcpp::export]]
arma::vec ProdVec(arma::vec E, arma::vec F) {
  arma::vec W = E % F;
  return W;
}

// [[Rcpp::export]]
double sumFac(Rcpp::NumericVector E){
  double acc = 0;
  for(int i =0; i < E.size(); i++){
    acc+=E[i];
  }
  return acc;
}

// [[Rcpp::export]]
double sebeta_const(double sigmae2, arma::vec E){
  double sum = sumFac(Rcpp::wrap(SquaredVec(E)));
  double tmp = sqrt (sigmae2/sum);
  return tmp;
}

// [[Rcpp::export]]
Rcpp::List ashBuild3(NumericMatrix Y,
                     NumericVector Ef,
                     NumericMatrix sigma2_mat){
  arma::vec Ef_arma = as<arma::vec>(Ef);
  arma::vec Ef2_arma = Ef_arma % Ef_arma;
  arma::mat Y_arma = as<arma::mat>(Y);
  arma::mat sigma2_mat_arma = as<arma::mat>(sigma2_mat);

  arma::vec tmp = (1.0/(sigma2_mat_arma))*Ef2_arma;
  arma::vec sebeta = sqrt(1.0/tmp);
  arma::vec Z = (Y_arma % (1.0/(sigma2_mat_arma)))*(Ef_arma);
  arma::vec betahat = Z/tmp;
  return List::create(_["betahat"]= betahat,
                      _["sebetahat"]=sebeta);
}

// [[Rcpp::export]]
Rcpp::List ashBuild1(Rcpp::NumericMatrix Y, Rcpp::NumericVector Ef, double sigma2){
  NumericMatrix Efmat = vec_to_mat(Ef);
  // arma::vec Ef2 = SquaredVec(Ef);
  NumericMatrix Efmat2 = vec_to_mat(Rcpp::wrap(SquaredVec(Ef)));
  arma::mat Efmat_arma = as<arma::mat>(Efmat);
  arma::mat Ymat = as<arma::mat>(Y);
  int K_Y = Y.ncol() ;
  int N_Ef = Efmat.nrow();
  if (N_Ef != K_Y){
    stop("Inadmissible dimensions of the matrices");
  }

  arma::mat tmp = Efmat_arma.t() * Ymat.t();
  double s = sumFac(Efmat2);
  arma::mat tmp2 = tmp/s;
  double sebeta = sebeta_const(sigma2, as<arma::vec>(Ef));

  return List::create(_["betahat"]= tmp2,
                      _["sebetahat"]=sebeta);
}

// [[Rcpp::export]]
Rcpp::List ashBuild2(Rcpp::NumericMatrix Y, Rcpp::NumericVector Ef, Rcpp::NumericVector sigma2_vec, bool isRow){
  if (isRow==TRUE){
    arma::vec Ef_arma = as<arma::vec>(Ef);
    arma::mat Y_arma = as<arma::mat>(Y);
    arma::vec sigma2_vec_arma = as<arma::vec>(sigma2_vec);

    arma::vec Ef2_arma = Ef_arma % Ef_arma;
    arma::vec invsigma2_vec_arma = 1.0/sigma2_vec_arma;
    arma::vec tmp = invsigma2_vec_arma % Ef2_arma;
    double tmp_sum = arma::sum(tmp);
    double sebeta = sqrt(1.0/tmp_sum);
    arma::vec tmp2 = invsigma2_vec_arma % Ef_arma;
    arma::mat betahat = (Y_arma * tmp2)/ tmp_sum;
    return List::create(_["betahat"]= betahat,
                        _["sebetahat"]=sebeta);
  }else{
    arma::vec sigma2_vec_arma = as<arma::vec>(sigma2_vec);
    double tmp_sum = arma::sum (SquaredVec(Ef));
    arma::vec sebeta = sqrt((sigma2_vec_arma)/tmp_sum);
    arma::mat Efmat_arma = as<arma::mat>(vec_to_mat(Ef));
    arma::mat Ymat = as<arma::mat>(Y);
    arma::mat betahat = (Efmat_arma.t() * Ymat.t())/tmp_sum;
    return List::create(_["betahat"]= betahat,
                        _["sebetahat"]=sebeta);
  }
}


// [[Rcpp::export]]

arma::mat SigmaV_Est (arma::mat Y, arma::vec El, arma::vec Ef)
{
  arma::mat Efmat = arma::mat(Ef);
  arma::mat Elmat = arma::mat(El);
  arma::mat Ef2mat = arma::mat(SquaredVec(Ef));
  arma::mat El2mat = arma::mat(SquaredVec(El));
  arma::mat res = (Y%Y) - 2*(Y % (Elmat * Efmat.t())) + (El2mat*Ef2mat.t());
  return (res);
}

// [[Rcpp::export]]

arma::mat SigmaV_Est_2 (arma::mat Y, arma::vec El, arma::vec Ef,
                        arma::vec El_listed, arma::vec Ef_listed)
{
  arma::mat Efmat = arma::mat(Ef);
  arma::mat Elmat = arma::mat(El);
  arma::mat Elmat_listed = arma::mat(El_listed);
  arma::mat Efmat_listed = arma::mat(Ef_listed);

  arma::mat Ef2mat = arma::mat(SquaredVec(Ef));
  arma::mat El2mat = arma::mat(SquaredVec(El));
  arma::mat El2mat_listed = arma::mat(SquaredVec(El_listed));
  arma::mat Ef2mat_listed = arma::mat(SquaredVec(Ef_listed));

  arma::mat Yhat = Y + Elmat_listed * Efmat_listed.t();
  arma::mat fl_norm1 = (Elmat*Efmat.t()) + (Elmat_listed * Efmat_listed.t());
  arma::mat sq_fl_norm1 = fl_norm1 % fl_norm1;
  //  arma::mat fl_norm2 = El2mat * Ef2mat.t() + El2mat_listed * Ef2mat_listed.t();
  arma::mat res = (sq_fl_norm1) + (Yhat%Yhat) - 2*(Yhat % ((Elmat * Efmat.t()) + (Elmat_listed * Efmat_listed.t())));
  return (res);
}


double  calcLik1 (NumericMatrix sigma2v,
                  NumericVector sigma2_true,
                  int n_threads)
{
  int N = sigma2v.nrow();
  int P = sigma2v.ncol();
  NumericVector res(P);
#pragma omp parallel for num_threads(n_threads)
  for( int j=0; j < P; j++ ) {
    res[j] = (-N/2)*((mean(sigma2v.column(j))/sigma2_true[j]) + log(2*M_PI*(sigma2_true[j])));
  }
  return(sum(res));
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

double  calcLik2 (NumericMatrix sigma2v,
                  NumericMatrix sigma2_true,
                  int n_threads)
{
  //  int N = sigma2v.nrow();
  int P = sigma2v.ncol();
  NumericVector outvec(P);
#pragma omp parallel for num_threads(n_threads)
  for(int j=0; j<P; j++){
    outvec[j] = sum(log(2*M_PI*(sigma2_true.column(j)))+ ((sigma2v.column(j))/(sigma2_true.column(j))));
  }
  double res = (-0.5)*(sum(outvec));
  return(res);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double  calcLik0 (NumericMatrix sigma2v, double sigma2_true){
  int N = sigma2v.nrow();
  int P = sigma2v.ncol();
  double mean_globe = (sum(sigma2v)/(N*P));
  double out = -((N*P)/2) * (log(2*M_PI*sigma2_true) +
                 (mean_globe)/sigma2_true + log((N*P)/2) - Rf_digamma((N*P)/2));
  return(out);
}

