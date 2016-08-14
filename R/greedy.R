#' Factor Loading Adaptive Shrinkage (K factors version using greedy algorithm)
#'
#' flash provide rank K matrix decomposition with greedy algorithm
#'
#' @param Y is the data matrix (N by P)
#' @param K is the max number of factor user want. the output will provide the actual number of factor
#'   automaticaly.
#' @param flash_para is the list for input parameters for flash
#' @param gvalue is the output style of greedy algorithm, 
#' "lik" provide the lowerbound of loglikelihood
#' "eigen" provide the sudo eigenvalue for each factor.
#'
#' @details greedy algorithm on the residual matrix to get a rank one matrix decomposition
#'
#' @export greedy
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
#' g = greedy(Y,10)
#'

greedy = function(Y,K,flash_para = list(),
                  plugin = FALSE,
                  gvalue = c("lik","eigen")){
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
  #initial the residual for the greedy algorithm
  residual = Y
  flash_default$Y = residual
  # update the input parameters for flash
  flash_para = modifyList(flash_default,flash_para,keep.null = TRUE)
  # this is the first rank one decomposition
  r_flash = do.call(flash,flash_para)
  # get the parameters and the gvalue which is for tracking the factors
  l_temp = r_flash$l
  f_temp = r_flash$f
  # keep the second moment in case
  l2_temp = r_flash$l2
  f2_temp = r_flash$f2
  priorpost = r_flash$obj_val - r_flash$c_lik_val
  clik = r_flash$c_lik_val
  priorpost_vec = c(priorpost)
  clik_vec = c(clik)
  # test whether it is zero
  if(sum(l_temp^2)==0 | sum(f_temp^2)==0){
    L_out = rep(0,N)
    F_out = rep(0,P)
    return(list(l = L_out,f = F_out))
  }else{
    L_out = l_temp
    F_out = f_temp
    # keep second moment
    L2_out = l2_temp
    F2_out = f2_temp
    # restore the mean and second moment into fl_list
    flash_para$fl_list$El = L_out
    flash_para$fl_list$Ef = F_out
    if(plugin == TRUE){
      flash_para$fl_list$El2 = L_out^2
      flash_para$fl_list$Ef2 = F_out^2
    }else{
      flash_para$fl_list$El2 = L2_out
      flash_para$fl_list$Ef2 = F2_out
    }
    #get the new residual
    residual = residual - l_temp %*% t(f_temp)
    #itreation for the rank K-1
    for(k in 2:K){
      cat("We are at iteration", k, "\n");
      # use the residual as the input of the next flash rank one
      flash_para$Y = residual
      #rank one for residual
      r_flash = do.call(flash,flash_para)
      # get the parameters and the gvalue which is for tracking the factors
      l_temp = r_flash$l
      f_temp = r_flash$f
      # try to output the sigmae2
     # print(r_flash$sigmae2)
      # get the new residual
      residual = residual - l_temp %*% t(f_temp)
      #check if we should stop at this iteration
      if(sum(l_temp^2)==0 | sum(f_temp^2)==0){
        break
      }else{
        # keep the second moment in case
        l2_temp = r_flash$l2
        f2_temp = r_flash$f2
        priorpost = r_flash$obj_val - r_flash$c_lik_val
        clik = r_flash$c_lik_val
        priorpost_vec = c(priorpost_vec, priorpost)
        clik_vec = c(clik_vec,clik)
        # if not go to next step and restore the l and f
        L_out = cbind(L_out,l_temp)
        F_out = cbind(F_out,f_temp)
        # for second moment 
        L2_out = cbind(L2_out,l2_temp)
        F2_out = cbind(F2_out,f2_temp)
        # restore the mean and second moment into fl_list
        flash_para$fl_list$El = L_out
        flash_para$fl_list$Ef = F_out
        if(plugin == TRUE){
          flash_para$fl_list$El2 = L_out^2
          flash_para$fl_list$Ef2 = F_out^2
        }else{
          flash_para$fl_list$El2 = L2_out
          flash_para$fl_list$Ef2 = F2_out
        }
      }
      # for loop can control the number of the factors if needed
    }
    # delete the column names for l f l2 f2
    colnames(L_out) = NULL
    colnames(F_out) = NULL
    colnames(L2_out) = NULL
    colnames(F2_out) = NULL
    return(list(l = L_out,f = F_out,
                l2 = L2_out, f2 = F2_out, 
                priorpost_vec = priorpost_vec, clik_vec =clik_vec))
  }
  # no return here, since return before
}

# eigen plot
#' title eigen plot
#'
#' description eigen plot
#'
#' @return return eigen value
#' @param g  greedy object which is a list
#' @param ssY sum square of Y
#' @keywords internal
#' 
eigen_plot = function(g,ssY){
  K = dim(g$l)[2]
  eigen_val = rep(0,K)
  for(k in 1:K){
    eigen_val[k] = sum( (g$l[,k] %*% t(g$f[,k]))^2 ) / ssY
  }
  return(eigen_val)
}
