% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/flash.R
\name{initial_value}
\alias{initial_value}
\title{inital value for flash}
\usage{
initial_value(Y, nonnegative = FALSE, factor_value = NA,
  fix_factor = FALSE, fl_list = list())
}
\arguments{
\item{Y}{the data matrix}

\item{nonnegative}{if the facotor and loading are nonnegative or not.
TRUE for nonnegative
FALSE for no constraint}

\item{factor_value}{is the factor value if the factor is fixed}

\item{fix_factor}{whether the factor is fixed or not
TRUE for fix_factor
FALSE for non-constraint}
}
\value{
list of factor, loading and variance of noise matrix
 \itemize{
  \item{\code{El}} {is a N vector for mean of loadings}
  \item{\code{El2}} {is a N vector for second moment of loadings}
  \item{\code{Ef}} {is a N vector for mean of factors}
  \item{\code{Ef2}} {is a N vector for second moment of factors}
  \item{\code{sigmae2_v}}{is a N by P matrix for residual square}
  \item{\code{sigmae2_true}}{is a N by P matrix for estimated value for the variance structure}
 }
}
\description{
description inital value for flash
}
\keyword{internal}

