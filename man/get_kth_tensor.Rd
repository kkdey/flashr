% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tbackfit.R
\name{get_kth_tensor}
\alias{get_kth_tensor}
\title{Form the rank-1 tensor from the kth components.}
\usage{
get_kth_tensor(factor_list, k)
}
\arguments{
\item{factor_list}{A list of matrices of numerics, each with the
same number of columns. The kth column of the nth element of
\code{factor_list} is the kth component of the nth mode of mean
tensor.}

\item{k}{A positive integer. The component to form.}
}
\description{
Form the rank-1 tensor from the kth components.
}
\author{
David Gerard
}

