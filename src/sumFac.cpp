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


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppEigen)]]

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

Rcpp::NumericVector inv_vec(Rcpp::NumericVector E){
  NumericVector invE = 1.0/E;
  return invE;
}

// [[Rcpp::export]]
double sebeta_const(double sigmae2, arma::vec E){
  double sum = sumFac(Rcpp::wrap(SquaredVec(E)));
  double tmp = sqrt (sigmae2/sum);
  return tmp;
}


// [[Rcpp::export]]
SEXP vec_to_mat(SEXP x_) {
  std::vector<double> x = as< std::vector<double> >(x_);
  /* ... do something with x ... */
  NumericVector output = wrap(x);
  output.attr("dim") = Dimension(x.size(), 1);
  return output;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]
Rcpp::List ashBuildpool(Rcpp::NumericMatrix Y, Rcpp::NumericVector Ef, arma::mat sigma2, bool isRow){
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
  
  int N_sig = sigma2.n_rows;
  int P_sig = sigma2.n_cols;
  
  
  if (N_sig==1 && P_sig==1){
    if(isRow==TRUE || isRow==FALSE){    
    double sigma2_float = as_scalar(sigma2);
    arma::mat tmp = Efmat_arma.t() * Ymat.t();
    double s = sumFac(Efmat2);
    arma::mat betahat = tmp/s;
    double sebeta = sebeta_const(sigma2_float, as<arma::vec>(Ef));
    return List::create(_["betahat"]= betahat,
                        _["sebetahat"]=sebeta);
    }
  }
  
  if (N_sig > 1 && P_sig == 1){
    
    if(N_sig !=N_Ef){
      stop("the length of sigma vector does not match with loadings vector");
    }
    
    if(isRow==TRUE){
    arma::vec sigma2_vec = vectorise(sigma2);
    arma::vec Ef2 = SquaredVec(Ef);
    arma::vec invsigma2_vec = 1.0/sigma2_vec;
    arma::vec tmp = ProdVec(invsigma2_vec, Ef2);
    double tmp_sum = arma::sum(tmp);
    double sebeta = sqrt(1.0/tmp_sum);
    arma::vec tmp2 = ProdVec(invsigma2_vec, Ef);
    arma::vec betahat = (as<arma::mat>(Y) * tmp2)/ tmp_sum;
    return List::create(_["betahat"]= betahat,
                        _["sebetahat"]=sebeta);
    }
    if(isRow==FALSE){
      arma::vec sigma2_vec = vectorise(sigma2);
      double tmp_sum = arma::sum (SquaredVec(Ef));
      arma::vec sebeta = sqrt((sigma2_vec)/tmp_sum);
      arma::mat Efmat_arma = as<arma::mat>(vec_to_mat(Ef));
      arma::mat Ymat = as<arma::mat>(Y);
      arma::vec betahat = (Efmat_arma.t() * Ymat.t())/tmp_sum;
      return List::create(_["betahat"]= betahat,
                          _["sebetahat"]=sebeta);
    }
  }
  
  if (N_sig > 1 && P_sig > 1){
    if(isRow==TRUE || isRow==FALSE){  
    arma::vec Ef2 = SquaredVec(Ef);
    arma::vec tmp = (1.0/(sigma2))*Ef2;
    arma::vec sebeta = sqrt(1.0/tmp);
    arma::vec Z = ((as<arma::mat>(Y)) % (1.0/(sigma2)))*(as<arma::vec>(Ef));
    arma::vec betahat = Z/tmp;
    return List::create(_["betahat"]= betahat,
                        _["sebetahat"]=sebeta);
    }}
  return 0; 
}



// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppEigen)]]

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
  
//  arma::mat Ytrans = arma::trans(Ymat);
//  arma::mat Eftrans = arma::trans(Efmat);
  arma::mat tmp = Efmat_arma.t() * Ymat.t();
  double s = sumFac(Efmat2);
  arma::mat tmp2 = tmp/s;
  double sebeta = sebeta_const(sigma2, as<arma::vec>(Ef));
//  arma::vec betahat_vec = as<arma::vec>(tmp2);
//  double sebeta =0.5;

  return List::create(_["betahat"]= tmp2,
                      _["sebetahat"]=sebeta);
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]



Rcpp::List ashBuild2(Rcpp::NumericMatrix Y, Rcpp::NumericVector Ef, Rcpp::NumericVector sigma2_vec, bool isRow){
    if (isRow==TRUE){
//  NumericMatrix Efmat = vec_to_mat(Ef);
  arma::vec Ef2 = SquaredVec(Ef);
//  NumericMatrix Efmat2 = vec_to_mat(Rcpp::wrap(Ef2));
  arma::vec invsigma2_vec = 1.0/sigma2_vec;
  arma::vec tmp = ProdVec(invsigma2_vec, Ef2);
  double tmp_sum = arma::sum(tmp);
  double sebeta = sqrt(1.0/tmp_sum);
  arma::vec tmp2 = ProdVec(invsigma2_vec, Ef);
//  arma::mat Ymat = as<arma::mat>(Y);
//  arma::mat tmp2_mat = as<arma::mat>(vec_to_mat(Rcpp::wrap(tmp2)));
//  arma::mat prodtemp = Ymat * tmp2_mat;
//  arma::mat betahat = prodtemp/tmp_sum;
  arma::mat betahat = (as<arma::mat>(Y) * tmp2)/ tmp_sum;
  return List::create(_["betahat"]= betahat,
                      _["sebetahat"]=sebeta);
    }else{
      double tmp_sum = arma::sum (SquaredVec(Ef));
      arma::vec sebeta = sqrt(as<arma::vec>(sigma2_vec)/tmp_sum);
      arma::mat Efmat_arma = as<arma::mat>(vec_to_mat(Ef));
      arma::mat Ymat = as<arma::mat>(Y);
      arma::mat betahat = (Efmat_arma.t() * Ymat.t())/tmp_sum;
      return List::create(_["betahat"]= betahat,
                          _["sebetahat"]=sebeta);
    }
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]

Rcpp::List ashBuild3(arma::mat Y, arma::vec Ef, arma::mat sigma2_mat){
  arma::vec Ef2 = SquaredVec(Ef);
  arma::vec tmp = (1.0/(sigma2_mat))*Ef2;
  arma::vec sebeta = sqrt(1.0/tmp);
  arma::vec Z = (Y % (1.0/(sigma2_mat)))*(Ef);
  arma::vec betahat = Z/tmp;
  return List::create(_["betahat"]= betahat,
                      _["sebetahat"]=sebeta);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]

Rcpp::List initval_postprocess1(arma::mat Y, Rcpp::NumericVector El, bool nonnegative){
  if(nonnegative==TRUE){
    for(int i=0; i< El.size(); i++){
      El[i] = std::max(El[i],0.0);
    }
  }
  
  arma::mat Elmat = as<arma::mat>(vec_to_mat(El));
  arma::vec El2 = SquaredVec(as<arma::vec>(El));
  arma::mat Elt = Elmat.t();
  arma::mat Ef = Elt * Y;
  arma::vec Ef2 = SquaredVec(arma::vectorise(Ef));
  return List::create(_["loading"]= El,
                      _["loading2"]=El2,
                      _["factor"]=arma::vectorise(Ef),
                      _["factor2"]=Ef2
  );
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]


Rcpp::List initval_postprocess_fixfactor(arma::mat Y, arma::vec fixF, bool nonnegative)
{
  arma::vec Ef = fixF;
  arma::vec Ef2 = SquaredVec(fixF);
  arma::vec El = (Y*Ef)/ (arma::sum(Ef2));
  if(nonnegative==TRUE){
    #pragma omp parallel for
    for(int i=0; i< El.size(); i++){
      El[i] = std::max(El[i],0.0);
    }
  }
  arma::vec El2 = SquaredVec(El);
  return List::create(_["loading"]= El,
                      _["loading2"]=El2,
                      _["factor"]=arma::vectorise(Ef),
                      _["factor2"]=Ef2
  );
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]

Rcpp::NumericMatrix SigmaV_Est (arma::mat Y,
                                arma::vec El,
                                arma::vec Ef)
{
  arma::mat Efmat = arma::mat(Ef);
  arma::mat Elmat = arma::mat(El);
  arma::mat Ef2mat = arma::mat(SquaredVec(Ef));
  arma::mat El2mat = arma::mat(SquaredVec(El));
  arma::mat res = (Y%Y) - 2*(Y % (Elmat * Efmat.t())) + (El2mat*Ef2mat.t());
  return Rcpp::wrap(res);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]

Rcpp::NumericMatrix SigmaV_Est_2 (arma::mat Y,
                                arma::vec El,
                                arma::vec Ef,
                                arma::vec El_listed,
                                arma::vec Ef_listed)
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
  arma::mat fl_norm = (Elmat*Efmat.t()) + (Elmat_listed * Ef2mat_listed.t());
  arma::mat res = fl_norm + (Y%Y) - 2*(Y % ((Elmat * Efmat.t()) + (Elmat_listed * Efmat_listed.t())));
  return Rcpp::wrap(res);
}

// [[Rcpp::export]]

double  calcLik0 (arma::mat sigma2v, double sigma2_true)
{
  int N = sigma2v.n_rows;
  int P = sigma2v.n_cols;
  double mean_globe = (arma::accu(sigma2v)/(N*P));
//  double mean2 = 0;
//  #pragma omp parallel for num_threads(4)
//  for(int j=0; j < sigma2v.n_cols; j++){
//    mean2 = mean2 + mean(sigma2v.col(j));
//  }
//  mean2 = mean2/(sigma2v.n_cols);
  double out = -((N*P)/2) * (log(2*M_PI*sigma2_true) + 
                   (mean_globe)/sigma2_true + log((N*P)/2) - Rf_digamma((N*P)/2));
  return(out);
}

// [[Rcpp::export]]

double mean_mat_rcpp (arma::mat sigma)
{
  double res = mean(mean(sigma));
  return(res);
}




// [[Rcpp::export]]

Rcpp::NumericVector rowmeans( Rcpp::NumericMatrix& X ) {
  
  int nRows = X.nrow();
  NumericVector out = no_init(nRows);
  #pragma omp parallel for num_threads(4)
  for( int i=0; i < nRows; i++ ) {
    NumericMatrix::Row tmp = X(i, _);
    out[i] = mean( tmp );
  }
  
  return out;
  
}

// [[Rcpp::export]]


Rcpp::NumericVector colmeans( Rcpp::NumericMatrix X ) {
  
  int nCols = X.ncol();
  NumericVector out = no_init(nCols);
  #pragma omp parallel for num_threads(4)
  for( int j=0; j < nCols; j++ ) {
    NumericMatrix::Column tmp = X(_, j);
    out[j] = mean( tmp );
  }
  
  return out;
}


/*double  calcLik1 (arma::mat sigma2v, arma::vec sigma2_true)
{
  int N = sigma2v.n_rows;
  int P = sigma2v.n_cols;
  arma::vec res = zeros<vec>(P);
//  arma::vec meancol = arma::vectorise(arma::mean(sigma2v, 0));
//  arma::vec res = (-N/2)*((meancol/sigma2_true) + log(2*M_PI*(sigma2_true)));
  #pragma omp parallel for num_threads(10)
  for( int j=0; j < P; j++ ) {
    res[j] = (-N/2)*((mean(sigma2v.col(j))/sigma2_true[j]) + log(2*M_PI*(sigma2_true[j])));
  }
  return(arma::accu(res));
} */

/* // [[Rcpp::export]]

double  calcLik2 (arma::mat sigma2v, arma::mat sigma2_true)
{
  arma::mat res1 = log(2*M_PI*(sigma2_true)) + (sigma2v % (1.0/sigma2_true));
  double out = (-0.5)*arma::accu(res1);
  return(out);
} */


/*double  calcLik22 (arma::mat sigma2v, arma::mat sigma2_true)
{
  int N = sigma2v.n_rows;
  int P = sigma2v.n_cols;
  arma::mat outmat(N,P, fill::zeros);
  #pragma omp parallel for num_threads(4)
    for(int j=0; j<P; j++){
      outmat.col(j) = log(2*M_PI*(sigma2_true.col(j)))+ ((sigma2v.col(j))/(sigma2_true.col(j)));
    }
  double res = (-0.5)*(arma::accu(outmat));
  return(res);
} 
 */

/*
system.time(threadcheck(4))

sumFac(as.vector(c(34,56, 23)))
sebeta_const(0.5, as.vector(c(34,56, 23)))
SquaredVec(as.vector(c(23,45,12)))

m <- c(1, 2, 3)
vec_to_mat(m)

Y <- matrix(rnorm(10000), nrow = 200);
e <- as.vector(rnorm(50,3,1))
sigma2 <- 0.5;

system.time(for(m in 1:1000){
res1 <- ashBuild1(Y, e, sigma2)
})

system.time(for(m in 1:1000){
betahat <- (t(e) %*% t(Y)) / (sum(e^2))
sum_e2 = sum(e^2)
sebeta = sqrt( sigma2/(sum_e2) )
res2 <- list("betahat"=betahat, "sebeta"=sebeta)
})

sigma2_vec <- as.vector(rchisq(50,4));

system.time(for(m in 1:1000){
  res3 <- ashBuild2(Y, e, sigma2_vec, isRow=TRUE)
})

Ef <-e;

system.time(for (m in 1:1000){
Ef2 <- Ef^2
sum_Ef2 = (1/sigma2_vec) * Ef2
sum_Ef2 = sum(sum_Ef2)
sebeta = sqrt(1/(sum_Ef2))
betahat = as.vector( Y %*% (Ef/sigma2_vec) ) / (sum_Ef2)
betahat=as.vector(betahat)
res4 <- list("betahat"=betahat, "sebeta"=sebeta)
})

system.time(for(m in 1:1000){
  res5 <- ashBuild2(Y, e, sigma2_vec, isRow=FALSE)
})

system.time(for (m in 1:1000){
  Ef2 <- Ef^2
  sum_Ef2 = sum(Ef2)
  sebeta = sqrt( sigma2_vec/(sum_Ef2) )
  betahat = (t(Ef) %*% t(Y)) / (sum_Ef2)
  betahat=as.vector(betahat)
  res6 <- list("betahat"=betahat, "sebeta"=sebeta)
})

  

ProdVec(c(3,6), c(5,7))
inv_vec(e)


sigma2_mat <- matrix(rchisq(10000,4), nrow=200);
Y <- matrix(rnorm(10000), nrow = 200);
Ef <- as.vector(rnorm(50,3,1))

system.time(for(m in 1:1000){
  res7 <- ashBuild3(Y, Ef, sigma2_mat)
})

system.time(for(m in 1:1000){
  Ef2 <- Ef^2;
  sum_Ef2 = (1/sigma2_mat) %*% Ef2
  sum_Ef2 = as.vector(sum_Ef2)
  sebeta = sqrt(1/(sum_Ef2))
  betahat = as.vector( (Y/sigma2_mat) %*% Ef ) / (sum_Ef2)
  betahat=as.vector(betahat)
  res8 <- list("betahat"=betahat, "sebeta"=sebeta)
})

nonnegative =TRUE;

El <- as.vector(rnorm(200,2,1));
system.time(for(m in 1:1000){
  res8 <- initval_postprocess1(Y,El, nonnegative);
})

system.time(for(m in 1:1000){
  if(nonnegative){
    El = abs(El)
  }
  El2 = El^2
  Ef = as.vector(t(El)%*%Y)
  Ef2 = Ef^2
  res9 <- list("loading"=El, "loading2"=El2, "factor"=Ef, "factor2"=Ef2)
})

nonnegative <- TRUE
fixF <- as.vector(rnorm(50,0,1));
system.time(for(m in 1:1000){
  res10 <- initval_postprocess_fixfactor(Y,fixF, nonnegative);
})

Y <- matrix(rnorm(10000), nrow = 200);
Ef <- as.vector(rnorm(50,3,1))
El <- as.vector(rnorm(200,2,1));
system.time(for(m in 1:1000){
  res11 <- SigmaV_Est(Y, El, Ef)
})

system.time(for(m in 1:1000){
  El2 <- El^2; Ef2 <- Ef^2;
  res12 <- Y^2 - 2*Y*(El %*% t(Ef)) + (El2 %*% t(Ef2))
})
  
sigma2v <- matrix(rchisq(10000,4), nrow=200);
sigma2_true <- 0.5;

system.time(for(m in 1:1000){
res13 <- calcLik0(sigma2v, sigma2_true)
})

system.time(for(m in 1:1000){
  N <- dim(sigma2v)[1]; P <- dim(sigma2v)[2];
  res14 <-  -(N*P)/2 * ( log(2*pi*sigma2_true) + (mean(sigma2v))/(sigma2_true) + log(N*P/2) - digamma(N*P/2))
})

rowmeans(sigma2v)
colmeans(sigma2v)

sigma2v <- matrix(rchisq(1e+8,4), nrow=200);
sigma2_true <- as.vector(rchisq(500000,3))
system.time(for(m in 1:1){
  res15 <- calcLik1(sigma2v, sigma2_true)
})
  
#system.time(for(m in 1:1000){
#  N <- dim(sigma2v)[1]; P <- dim(sigma2v)[2];
#  tmp <- colMeans(sigma2v)
#  res16 <-  -(N/2) * sum( log(2*pi*sigma2_true) + (tmp)/(sigma2_true) )
#})
  
  
sigma2v <- matrix(rchisq(10000,4), nrow=200);
sigma2_true <- matrix(rchisq(10000,0.5), nrow=200);
system.time(for(m in 1:1000){
  res17 <- calcLik2(sigma2v, sigma2_true)
})

system.time(for(m in 1:1000){
  res17b <- calcLik22(sigma2v, sigma2_true)
})

system.time(for(m in 1:1000){
  N <- dim(sigma2v)[1]; P <- dim(sigma2v)[2];
  res18 <-  -(1/2) * sum( log(2*pi*sigma2_true) + (sigma2v)/(sigma2_true) )
})

*/