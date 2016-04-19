#' backfitting to correct K-factor model
#'
#' Correct the factor and loading matrix estimates using backfitting algorithm
#'
#' @param Y is the data matrix (N by P)
#' @param intial_list is the list from intial_list algorithm
#' @param Fest is estimate for f to correct
#' @param maxiter_bf maximum number of iterations
#' @param flash_parais the list for flash parameters setting up
#' @param  gvalue is for output of the backfit,
#' "eigen" means just provide the sum of square
#' "lik" mean provide the lowerbound
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
#'   \item{\code{sigmae2}} {is variance for the error, scalar for constant case and vector for nonconstant case}
#'  }
#' @examples
#' N = 100
#' P = 200
#' Y = matrix(rnorm(N*P,0,1),ncol=P)
#' g = intial_list(Y,K = 10)
#' gb = backfitting(Y,g$l,g$f,maxiter_bf=100,maxiter_r1 = 5)
#'
backfitting = function(Y,intial_list, maxiter_bf=100,
                       flash_para = list(), gvalue = c("lik","eigen")){
  epsilon = 1
  tau = 1
  N = dim(Y)[1]
  P = dim(Y)[2]
  # match the input parameter
  gvalue = match.arg(gvalue, c("lik","eigen"))
  # set the default value for flash
  flash_default = list(tol=1e-5, maxiter_r1 = 500,
                       partype = "constant", 
                       sigmae2_true = NA, 
                       factor_value = NA,fix_factor = FALSE,
                       nonnegative = FALSE,
                       objtype = "margin_lik",
                       ash_para = list(),
                       fl_list = list())
  if(gvalue == "lik"){
    flash_default$objtype = "lowerbound_lik"
  }
  # this is never used but just a initail value
  flash_default$Y = Y
  flash_para = modifyList(flash_default,flash_para)
  #initial check
  if(is.vector(intial_list$l)){intial_list$l = matrix(intial_list$l,ncol=1)}
  if(is.vector(intial_list$f)){intial_list$f = matrix(intial_list$f,ncol=1)}
  if(nrow(intial_list$l)!=nrow(Y)){stop("L of wrong dimension for Y")}
  if(nrow(intial_list$f)!=ncol(Y)){stop("F of wrong dimension for Y")}
  # now we need to specify the fl_list for each
  # at first we put all the result from intial_list into the Lest Fest
  Lest = intial_list$l
  Fest = intial_list$f
  L2est = intial_list$l2
  F2est = intial_list$f2
  priorpost_vec = intial_list$priorpost_vec
  # print((sum(priorpost_vec) + intial_list$clik_vec[length(intial_list$clik_vec)]))
  # print(priorpost_vec)
  # track the obj value
  track_obj = c((sum(priorpost_vec) + intial_list$clik_vec[length(intial_list$clik_vec)]))
  # in the begining we have all the factors storing in the fl_list
  while(epsilon>1e-5 & tau < maxiter_bf){
    tau = tau + 1
    K = dim(Lest)[2]
    if(K==0){break} #tests for case where all factors disappear!
    # this one can be put out of the while loop
    preRMSfl = sqrt(mean((Lest %*% t(Fest))^2))
    sigmae2_out = rep(0,K)
    for(k in 1:K){
      residual = Y - Lest[,-k] %*% t(Fest[,-k])
      flash_para$Y = residual
      flash_para$fl_list$El = Lest[,-k]
      flash_para$fl_list$Ef = Fest[,-k]
      flash_para$fl_list$Ef2 = F2est[,-k]
      flash_para$fl_list$El2 = L2est[,-k]
      # run the rank one flash
      r_flash = do.call(flash,flash_para)
      Lest[,k] = r_flash$l
      Fest[,k] = r_flash$f
      L2est[,k] = r_flash$l2
      F2est[,k] = r_flash$f2
      sigmae2_out[k] = r_flash$sigmae2
      priorpost = r_flash$obj - r_flash$c_lik_val
      priorpost_vec[k] = priorpost
      clik = r_flash$c_lik_val
      obj_lik = clik + sum(priorpost_vec)
      track_obj = c(track_obj,obj_lik)
      # print(obj_lik)
      # print(priorpost_vec)
      # print(sqrt(mean((Y - Lest %*% t(Fest) -E)^2)) / sqrt(mean((Y - E)^2)))
    }
    #sigmae2_out = sigmae2_out/K
    # print(sigmae2_out)
    # remove the zero in the l and f
    # zeros = is_zero_factor(Lest) || is_zero_factor(Fest)  #this is wrong
    zeros = is_zero_factor(Lest)
    Lest = Lest[,!zeros,drop=FALSE]
    Fest = Fest[,!zeros,drop=FALSE]
    L2est = L2est[,!zeros,drop=FALSE]
    F2est = F2est[,!zeros,drop=FALSE]
    priorpost_vec = priorpost_vec[!zeros,drop = FALSE]
    RMSfl = sqrt(mean((Lest %*% t(Fest))^2))
    epsilon = abs(preRMSfl - RMSfl)
  }
  return(list(l = Lest, f = Fest , sigmae2 = sigmae2_out,track_obj = track_obj))
}

# returns whether a vector is all 0
allzeros = function(x){return (all(x==0))}

#returns vector of TRUE and FALSE for which factors (columns) of l are all 0
is_zero_factor=function(l){
  apply(l,2,allzeros)
}
