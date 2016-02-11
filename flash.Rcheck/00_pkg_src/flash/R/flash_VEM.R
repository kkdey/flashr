
# El is expectation of l, and El2 is the second moment of l, sigmae2 is estimation of sigmae^2
ATM_f = function(Y,El,El2,sigmae2){
  sum_El2 = sum(El2)
  sebeta = sqrt( sigmae2/(sum_El2) )
  betahat = (t(El) %*% Y) / (sum_El2)
  # betahat=(sum(l^2))^(-1)*(t(l)%*%Y)
  betahat=as.vector(betahat)
  ATM = ash(betahat, sebeta, method="fdr", mixcompdist="normal")
  Ef = ATM$PosteriorMean
  SDf = ATM$PosteriorSD
  Ef2 = SDf^2 + Ef^2
  return(list(Ef = Ef, Ef2 = Ef2))
}

ATM_l = function(Y,Ef,Ef2,sigmae2){
  sum_Ef2 = sum(Ef2)
  sebeta = sqrt(sigmae2/(sum_Ef2))
  betahat = (t(Ef) %*% t(Y)) / (sum_Ef2)
  # betahat=(sum(f^2))^(-1)*(t(f)%*%t(Y))
  betahat=as.vector(betahat)
  ATM = ash(betahat, sebeta, method="fdr", mixcompdist="normal")
  El = ATM$PosteriorMean
  SDl = ATM$PosteriorSD
  El2 = SDl^2 + El^2
  return(list(El = El, El2 = El2))
}


# Fval = function()







#' Multivariate Adaptive Shrinkage (original version)
#'
#' mvash provide adaptive shrinkage estimation of effect size for correlated data
#'
#' @param Y zzz
#' @param tol zzz
#' @param numtau zzz
#'
#' @details mvash provide shrinkage estimation of effect size for correlated data with
#' mixed normal prior of which the grids (standard deviation) are predefined. The proportion
#' for each component in the prior can be initialized as arbitury value. In updating the parameters
#' and hyper parameters, mvash uses Variantional Bayes procedure.
#'
#' @export flash_VEM
#'
#' @return list of posterior mean, standard deviation, proportion for each mixed component
#' and the estimation of the hyper patameters pi_hat in the prior
#'  \itemize{
#'   \item{\code{mu}} {is a P by K matrix of posterior mean for each gene on each component}
#'   \item{\code{s}} {is a P by K matrix of posterior standard deviation for each gene on each component}
#'   \item{\code{alpha}} {is a P by K matrix of posterior proportion for each gene on each component}
#'   \item{\code{pihat}} {is a estimation for the hyper parameters pi in the prior which the proportion
#'     for each fixed grids.}
#'  }
#' @examples
#' NULL

# set the number of iteration as numtau
flash_VEM = function(Y,tol=1e-6,numtau = 500){
  #dealing with missing value
  Y[is.na(Y)] = 0
  # get initial value for l and f and sigmae
  El = svd(Y)$u[,1]
  El2 = El^2
  Ef = as.vector(t(El)%*%Y)
  Ef2 = Ef^2

  #start iteration
  sigmae2_v = mean( Y^2 - 2*Y*(El %*% t(Ef)) + (El2 %*% t(Ef2)) )

  par_f = ATM_f(Y,El,El2,sigmae2_v)
  Ef = par_f$Ef
  Ef2 = par_f$Ef2

  sigmae2_v = mean( Y^2 - 2*Y*(El %*% t(Ef)) + (El2 %*% t(Ef2)) )
  #sigmae2
  par_l = ATM_l(Y,Ef,Ef2,sigmae2_v)
  El = par_l$El
  El2 = par_l$El2

  epsilon = 1
  tau = 1
  while(epsilon >= tol & tau < numtau ){
    tau = tau + 1
    pre_sigmae2 = sigmae2_v

    sigmae2_v = mean( Y^2 - 2*Y*(El %*% t(Ef)) + (El2 %*% t(Ef2)) )

    par_f = ATM_f(Y,El,El2,sigmae2_v)
    Ef = par_f$Ef
    Ef2 = par_f$Ef2
    if(sum(Ef^2)==0){
      l = 0
      f = 0
      break
    }
    sigmae2_v = mean( Y^2 - 2*Y*(El %*% t(Ef)) + (El2 %*% t(Ef2)) )
    #sigmae2
    par_l = ATM_l(Y,Ef,Ef2,sigmae2_v)
    El = par_l$El
    El2 = par_l$El2
    if(sum(El^2)==0){
      l = 0
      f = 0
      break
    }
    epsilon = abs(pre_sigmae2 - sigmae2_v )
  }
  return(list(l = El, f = Ef, sigmae2 = sigmae2_v))
}
