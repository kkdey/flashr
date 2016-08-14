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
// [[Rcpp::depends(RcppEigen)]]

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


/***R
system.time(threadcheck(4))
*/

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

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

/*

sigmae2_v <- matrix(rchisq(10000, 4), nrow = 100);
sigmae2_true <- matrix(rchisq(10000, 1), nrow = 100);

dim(sigmae2_true)
dim(sigmae2_v)
res1 <- calcLik2(sigmae2_v, sigmae2_true, n_threads=1)
N <- dim(sigmae2_v)[1];
P <- dim(sigmae2_v)[2];
res <- C_likelihood(N,P,sigmae2_v,sigmae2_true)
*/