% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/flash.R
\name{one_step_update}
\alias{one_step_update}
\title{title one step update in flash iteration using ash}
\usage{
one_step_update(Y, El, El2, Ef, Ef2, N, P, sigmae2_v, sigmae2_true, sigmae2,
  nonnegative = FALSE, partype = "constant", objtype = "margin_lik",
  fix_factor = FALSE, ash_para = list(), fl_list = list())
}
\arguments{
\item{Y}{the data matrix}

\item{El}{mean for the loadings}

\item{El2}{second moment for the loadings}

\item{Ef}{mean for the factors}

\item{Ef2}{second moment for the factors}

\item{N}{dimension of Y}

\item{P}{dimension of Y}

\item{sigmae2_v}{residual square}

\item{sigmae2_true}{the (true) known variance structure
Here, sigmae2 is the estimated variance structure in each step
sigmae2_true is the truth we know, some times sigmae2 is noisy version of sigmae2_true}

\item{sigmae2}{the estimation of the variance structure}

\item{nonnegative}{if the facotor and loading are nonnegative or not.
TRUE for nonnegative
FALSE for no constraint}

\item{partype}{parameter type for the variance,
"constant" for constant variance,
"var_col" for nonconstant variance for column,
"known" for the kown variance,
"Bayes_var" for Bayes version of the nonconstant variance for row and column
"loganova" is anova estiamtion for the log residual square}

\item{objtype}{objective function type,
"margin_lik" for conditional likelihood,
"lowerbound_lik" for full objective function}

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
  \item{\code{obj_val}}{the value of objectice function}
 }
}
\description{
description one step update in flash iteration using ash
}
\keyword{internal}

