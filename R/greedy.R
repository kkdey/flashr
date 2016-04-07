#' Factor Loading Adaptive Shrinkage (K factors version using greedy algorithm)
#'
#' flash provide rank K matrix decomposition with greedy algorithm
#'
#' @param Y is the data matrix (N by P)
#' @param K is the max number of factor user want. the output will provide the actual number of factor
#'   automaticaly.
#' @param r1_type is the flag to choose the constant variance for columns or not, defulat is constant
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

greedy = function(Y,K,r1_type = "constant"){
  N = dim(Y)[1]
  P = dim(Y)[2]
  #initial the residual for the greedy algorithm
  residual = Y
  # do rank one decomposition
  if(r1_type == "nonconstant"){
    r_flash = flash_r1c(residual)
  }else{
    r_flash = flash_r1(residual)
  }
  l_temp = r_flash$l
  f_temp = r_flash$f
  # test whether it is zero
  if(sum(l_temp^2)==0 | sum(f_temp^2)==0){
    L_out = 0
    F_out = 0
    return(list(l = L_out,f = F_out))
  }else{
    L_out = l_temp
    F_out = f_temp
    #get the new residual
    residual = residual - l_temp %*% t(f_temp)
    #itreation for the rank K-1
    for(k in 2:K){
      #rank one for residual
      if(r1_type == "nonconstant"){
        r_flash = flash_r1c(residual)
      }else{
        r_flash = flash_r1(residual)
      }
      l_temp = r_flash$l
      f_temp = r_flash$f
      # get the new residual
      residual = residual - l_temp %*% t(f_temp)
      #check if we should stop at this iteration
      if(sum(l_temp^2)==0 | sum(f_temp^2)==0){
        break
      }else{
        # if not go to next step and restore the l and f
        L_out = cbind(L_out,l_temp)
        F_out = cbind(F_out,f_temp)
      }
      # for loop can control the number of the factors if needed
    }
    return(list(l = L_out,f = F_out))
  }
  # no return here, since return before
}
