#include <RcppArmadillo.h>
#include <cmath>
#include <thread>
#include <typeinfo>
#include <omp.h>
#include <chrono>
#include <vector>
#include <stdlib.h> 

// For Rcpp+OpenMP parallelization, check: https://wbnicholson.wordpress.com/2014/07/10/parallelization-in-rcpp-via-openmp/
// For Armadillo Support, check: http://arma.sourceforge.net/docs.html

// LDFLAGS:  -L/usr/local/opt/openblas/lib
// CPPFLAGS: -I/usr/local/opt/openblas/include

using namespace Rcpp;
using namespace arma;
using namespace std;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]

arma::vec ProdVec(arma::rowvec v1, arma::rowvec v2){
  if(v1.n_elem != v2.n_elem){
    stop("vector lengths do not match");
  }

  arma::vec out = zeros<vec>(v1.n_elem);
  #pragma omp parallel for num_threads(100)
  for(int i=0; i < v1.n_elem; i++){
       out[i]=v1[i]*v2[i];
  }
  return out;
}

// [[Rcpp::export]]

arma::mat ProdMat(arma::mat m1, arma::mat m2){
  arma::mat out = zeros<mat>(m1.n_rows, m1.n_cols);
  #pragma omp parallel for num_threads(100)
  for(int i=0; i < m1.n_rows; i++){
    for(int j=0; j < m1.n_cols; j++){
      out(i,j)= m1(i,j)*m2(i,j);
    }
  }
  return out;
}
  

arma::vec SquaredVec(arma::vec Efac) {
  arma::vec Efac2 = Efac % Efac;
  return Efac2;
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]

Rcpp::List ashBuild3(arma::mat Y, arma::vec Ef, arma::mat sigma2_mat){
  arma::vec Ef2 = SquaredVec(Ef);
  arma::mat inv_sigma2_mat = (1.0/sigma2_mat);
  arma::vec tmp = inv_sigma2_mat*Ef2;
  arma::vec sebeta = sqrt(1.0/tmp);
  arma::vec betahat = ((Y % inv_sigma2_mat)*(Ef))/tmp;
 // arma::vec betahat = Z/tmp;
  return List::create(_["betahat"]= betahat,
                      _["sebetahat"]=sebeta);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]

Rcpp::List ashBuild2(arma::mat Y, arma::vec Ef, arma::vec sigma2_vec, bool isRow){
  if (isRow==TRUE){
    arma::vec Ef2 = SquaredVec(Ef);
    arma::vec invsigma2_vec = 1.0/sigma2_vec;
    arma::vec tmp = invsigma2_vec % Ef2;
    double tmp_sum = arma::sum(tmp);
    double sebeta = sqrt(1.0/tmp_sum);
    arma::vec tmp2 = tmp % Ef;
    arma::mat betahat = (Y * tmp2)/ tmp_sum;
    return List::create(_["betahat"]= betahat,
                        _["sebetahat"]=sebeta);
  }else{
    double tmp_sum = arma::sum (SquaredVec(Ef));
    arma::vec sebeta = sqrt((sigma2_vec)/tmp_sum);
    arma::mat betahat = (Ef.t() * Y.t())/tmp_sum;
    return List::create(_["betahat"]= betahat,
                        _["sebetahat"]=sebeta);
  }
}



/*** R

sigma2_mat <- matrix(rchisq(100000,4), nrow=20);
Y <- matrix(rnorm(100000), nrow = 20);
Ef <- as.vector(rnorm(5000,3,1))

system.time(res <- ProdVec(Y[1,], sigma2_mat[1,]))
  
  system.time(for(m in 1:1){
    res7 <- ashBuild3(Y, Ef, sigma2_mat)
  })
  
  system.time(for(m in 1:1){
    Ef2 <- Ef^2;
    sum_Ef2 = (1/sigma2_mat) %*% Ef2
      sum_Ef2 = as.vector(sum_Ef2)
      sebeta = sqrt(1/(sum_Ef2))
      betahat = as.vector( (Y/sigma2_mat) %*% Ef ) / (sum_Ef2)
      betahat=as.vector(betahat)
      res8 <- list("betahat"=betahat, "sebeta"=sebeta)
  })
  
  
 # sigma2_mat <- matrix(rchisq(100000,4), nrow=20);
  Y <- matrix(rnorm(1000000), nrow = 20000);
  Ef <- as.vector(rnorm(50,3,1))
  sigma2_vec <- as.vector(rchisq(50,4));
  
  system.time(for(m in 1:1000){
    res3 <- ashBuild2(Y, Ef, sigma2_vec, isRow=TRUE)
  })
    
  system.time(for (m in 1:1000){
    Ef2 <- Ef^2
    sum_Ef2 = (1/sigma2_vec) * Ef2
    sum_Ef2 = sum(sum_Ef2)
    sebeta = sqrt(1/(sum_Ef2))
    betahat = as.vector( Y %*% (Ef/sigma2_vec) ) / (sum_Ef2)
    betahat=as.vector(betahat)
    res4 <- list("betahat"=betahat, "sebeta"=sebeta)
  })
    
  
  
*/
