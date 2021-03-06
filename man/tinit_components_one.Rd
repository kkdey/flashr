% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tflash.R
\name{tinit_components_one}
\alias{tinit_components_one}
\title{Obtain initial estimates of each mode's components.}
\usage{
tinit_components_one(Y)
}
\arguments{
\item{Y}{An array of numerics.}
}
\value{
\code{ex_list} A list of vectors of the starting expected
    values of each component.


    \code{ex2_vec} A vector of starting values of \eqn{E[x'x]}
}
\description{
The inital values of each component are taken to be proportional to
the first singular vector of each mode-k matricization of the data
array. The constant is placed entirely on the first mode.
}
\author{
David Gerard
}

