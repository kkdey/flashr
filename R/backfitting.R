#' backfitting to correct K-factor model
#'
#' Correct the factor and loading matrix estimates using backfitting algorithm
#'
#' @param Y is the data matrix (N by P)
#' @param Lest is estimate for l to correct
#' @param Fest is estimate for f to correct
#' @param maxiter_bf maximum number of iterations
#' @param maxiter_r1 maximum number of iterations in flash rank one.
#' @param r1_type is the flag to choose the constant variance for columns or not, defulat is constant
#'
#' @details Repeatedly applies rank 1 algorithm to Y-L[,-i]F[,-i]'
#'
#' @export backfitting
#'
#' @importFrom ashr ash
#'
#' @return list of factor, loading and variance of noise matrix
#'  \itemize{
#'   \item{\code{l}} {is a N by K matrix for loadings}
#'   \item{\code{f}} {is a P by K matrix for factors}
#'  }
#' @examples
#' N = 100
#' P = 200
#' Y = matrix(rnorm(N*P,0,1),ncol=P)
#' g = greedy(Y,K = 10)
#' gb = backfitting(Y,g$l,g$f,maxiter_bf=100,maxiter_r1 = 5)
#'
backfitting = function(Y,Lest,Fest,maxiter_bf=100,maxiter_r1 = 500,r1_type ="constant"){
  # backfitting with initial values
  epsilon = 1
  tau = 1
  if(is.vector(Lest)){Lest = matrix(Lest,ncol=1)}
  if(is.vector(Fest)){Fest = matrix(Fest,ncol=1)}
  if(nrow(Lest)!=nrow(Y)){stop("Lest of wrong dimension for Y")}
  if(nrow(Fest)!=ncol(Y)){stop("Fest of wrong dimension for Y")}

  while(epsilon>1e-5 & tau < maxiter_bf){
    tau = tau + 1
    K = dim(Lest)[2]
    if(K==0){break} #tests for case where all factors disappear!
    # this one can be put out of the while loop
    residual = Y - Lest %*% t(Fest)
    preRMSfl = sqrt(mean((Lest %*% t(Fest))^2))
    for(i in 1:K){
      residual = residual + Lest[,i] %*% t(Fest[,i])
      if(r1_type == "nonconstant"){
        r_flash = flash_r1c(residual,maxiter_r1 = maxiter_r1)
      }else{
        r_flash = flash_r1(residual,maxiter_r1 = maxiter_r1)
      }
      Lest[,i] = r_flash$l
      Fest[,i] = r_flash$f
      residual = residual - Lest[,i] %*% t(Fest[,i])
    }
    # remove the zero in the l and f
    zeros = is_zero_factor(Lest) || is_zero_factor(Fest)
    Lest = Lest[,!zeros,drop=FALSE]
    Fest = Fest[,!zeros,drop=FALSE]
    RMSfl = sqrt(mean((Lest %*% t(Fest))^2))
    epsilon = abs(preRMSfl - RMSfl)
  }
  return(list(l = Lest, f = Fest))
}

# returns whether a vector is all 0
allzeros = function(x){return (all(x==0))}

#returns vector of TRUE and FALSE for which factors (columns) of l are all 0
is_zero_factor=function(l){
  apply(l,2,allzeros)
}
