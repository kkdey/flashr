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
// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]


arma::vec SquaredVec(arma::vec Efac) {
  arma::vec Efac2 = Efac % Efac;
  return Efac2;
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

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppEigen)]]

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

/***R
Y <- matrix(rnorm(1000000), nrow = 200);
Ef <- as.vector(rnorm(5000,3,1))
  El <- as.vector(rnorm(200,2,1));
system.time(for(m in 1:10){
  res11 <- SigmaV_Est(Y, El, Ef)
})
  
  system.time(for(m in 1:10){
    El2 <- El^2; Ef2 <- Ef^2;
    res12 <- Y^2 - 2*Y*(El %*% t(Ef)) + (El2 %*% t(Ef2))
  })
sum(abs(res11 - res12))

Y <- matrix(rnorm(1000000), nrow = 200);
Ef <- as.vector(rnorm(5000,3,1))
El <- as.vector(rnorm(200,2,1));
Ef_listed <- as.vector(rnorm(5000, 10,1));
El_listed <- as.vector(rnorm(200, 1,1));
  
system.time(for(m in 1:10){
  res13 <- SigmaV_Est_2(Y, El, Ef, El_listed, Ef_listed)
})

system.time(for(m in 1:10){
  El2 <- El^2; Ef2 <- Ef^2;
  El2_listed <- El_listed^2; Ef2_listed <- Ef_listed^2;
  Yhat = Y + El_listed %*% t(Ef_listed)
  # the residual matrix should be 
# this is for E(l_1f_1+l_2f_2+...+l_kf_k)^2
    fl_norm = (El%*%t(Ef) + El_listed%*%t(Ef_listed))^2 - 
      (El^2 %*% t(Ef^2) + (El_listed)^2 %*% t(Ef_listed)^2) +
      (El2 %*% t(Ef2) + El2_listed %*% t(Ef2_listed)) ### ???????????
    
    sigmae2_v = Yhat^2 - 2* Yhat * (El%*%t(Ef) + El_listed%*%t(Ef_listed)) + fl_norm
    res14 <- sigmae2_v
    })

sum(abs(res13 - res14))

*/


