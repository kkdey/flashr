% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tkron.R
\name{diag_mle}
\alias{diag_mle}
\title{MLE of mean zero Kronecker variance model.}
\usage{
diag_mle(R, itermax = 100, tol = 10^-3, homo_modes = NULL)
}
\arguments{
\item{R}{An array of numerics.}

\item{itermax}{A positive integer. The maximium number of
iterations to perform.}

\item{tol}{A positive numeric. The stopping criterion.}

\item{homo_modes}{A vector of integers. If \code{var_type =
"kronecker"} then \code{homo_modes} indicates which modes are
assumed to be homoscedastic.}
}
\value{
\code{esig_list} A list of vectors of numerics. The MLEs of
    the precisions, not the variances.
}
\description{
MLE of mean zero Kronecker variance model.
}
\author{
David Gerard
}

