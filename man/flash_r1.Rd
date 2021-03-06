% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/flash_r1.R
\name{flash_r1}
\alias{flash_r1}
\title{Factor Loading Adaptive Shrinkage (VEM version)}
\usage{
flash_r1(Y, tol = 1e-06, maxiter_r1 = 500)
}
\arguments{
\item{Y}{is the data matrix (N by P)}

\item{tol}{which is the stop criterion for the convergence, default is 1e-5}

\item{maxiter_r1}{number of iteration, default is 500. for the backfitting case, the number of tau should be 5 or 10.}
}
\value{
list of factor, loading and variance of noise matrix
 \itemize{
  \item{\code{l}} {is a N vector for loadings}
  \item{\code{f}} {is a P vector for factors}
  \item{\code{sigmae2}} {is mean of sigma square which is estimation for the noise variance}
 }
}
\description{
flash provide rank one matrix decomposition
}
\details{
flash_r1 privide rank one matrix decomposition with variational EM algorithm.
}
\examples{
N = 100
P = 200
Y = matrix(rnorm(N*P,0,1),ncol=P)
g = flash_r1(Y)

}

