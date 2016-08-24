#' backfitting to correct K-factor model
#'
#' Correct the factor and loading matrix estimates using backfitting algorithm
#'
#' @param Y is the data matrix (N by P)
#' @param initial_list is the list from initial_list algorithm
#' @param Fest is estimate for f to correct
#' @param maxiter_bf maximum number of iterations
#' @param flash_parais the list for flash parameters setting up
#' @param initial_list initial list for the starting value which contains
#' l,f,l2,f2
#' priorpost_vec
#' clik_vec
#' @param  gvalue is for output of the backfit,
#' "eigen" means just provide the sum of square
#' "lik" mean provide the lowerbound
#'
#' @details Repeatedly applies rank 1 algorithm to Y-L[,-i]F[,-i]'
#'
#' @export flash.backfitting
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
#' g = initial_list(Y,K = 10)
#' gb = flash.backfitting(Y,g$l,g$f,maxiter_bf=100,maxiter_r1 = 5)
#'
flash.backfitting = function(Y,initial_list, maxiter_bf=100,
                       flash_para = list(), gvalue = c("lik","eigen"),
                       parallel = FALSE){
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
  flash_para = modifyList(flash_default,flash_para, keep.null = TRUE)
  #initial check
  if(is.vector(initial_list$l)){initial_list$l = matrix(initial_list$l,ncol=1)}
  if(is.vector(initial_list$f)){initial_list$f = matrix(initial_list$f,ncol=1)}
  if(nrow(initial_list$l)!=nrow(Y)){stop("L of wrong dimension for Y")}
  if(nrow(initial_list$f)!=ncol(Y)){stop("F of wrong dimension for Y")}
  # now we need to specify the fl_list for each
  # at first we put all the result from initial_list into the Lest Fest
  Lest = initial_list$l
  Fest = initial_list$f
  L2est = initial_list$l2
  F2est = initial_list$f2
  priorpost_vec = initial_list$priorpost_vec
  # print((sum(priorpost_vec) + initial_list$clik_vec[length(initial_list$clik_vec)]))
  # print(priorpost_vec)
  # track the obj value
  obj_lik = (sum(priorpost_vec) + initial_list$clik_vec[length(initial_list$clik_vec)])
  track_obj = c(obj_lik)
  # in the begining we have all the factors storing in the fl_list
  while(epsilon>1e-5 & tau < maxiter_bf){
    tau = tau + 1
    K = dim(Lest)[2]
    if(K==0){break} #tests for case where all factors disappear!
    # this one can be put out of the while loop
    # preRMSfl = sqrt(mean((Lest %*% t(Fest))^2))
    pre_obj_lik = obj_lik
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
      priorpost = r_flash$obj_val - r_flash$c_lik_val
      priorpost_vec[k] = priorpost
      clik = r_flash$c_lik_val
      obj_lik = clik + sum(priorpost_vec)
      track_obj = c(track_obj,obj_lik)
    }
    #sigmae2_out = sigmae2_out/K
    print(sigmae2_out)
    # remove the zero in the l and f
    # zeros = is_zero_factor(Lest) || is_zero_factor(Fest)  #this is wrong
    zeros = is_zero_factor(Lest)
    Lest = Lest[,!zeros,drop=FALSE]
    Fest = Fest[,!zeros,drop=FALSE]
    L2est = L2est[,!zeros,drop=FALSE]
    F2est = F2est[,!zeros,drop=FALSE]
    priorpost_vec = priorpost_vec[!zeros,drop = FALSE]
    # RMSfl = sqrt(mean((Lest %*% t(Fest))^2))
    epsilon = abs(obj_lik - pre_obj_lik)
  }
  return(list(l = Lest, f = Fest , sigmae2 = sigmae2_out,track_obj = track_obj))
}

# returns whether a vector is all 0
allzeros = function(x){return (all(x==0))}

#returns vector of TRUE and FALSE for which factors (columns) of l are all 0
is_zero_factor=function(l){
  apply(l,2,allzeros)
}
