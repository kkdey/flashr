#' title ash type model for f
#'
#' description use ash type model to maxmization
#'
#' @return Ef is the mean of f Ef2 is mean of f^2
#'
#' @keywords internal

# El is expectation of l, and El2 is the second moment of l, sigmae2 is estimation of sigmae^2 which is a P-vector
ATM_f_c = function(Y,El,El2,sigmae2){
  # sigmae2 is a P-vector
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

#' title ash type model for l
#'
#' description use ash type model to maxmization
#'
#' @return El is the mean of l,  El2 is mean of l^2
#'
#' @keywords internal
#'
ATM_l_c = function(Y,Ef,Ef2,sigmae2){
  sum_Ef2 = (1/sigmae2) * Ef2
  sum_Ef2 = sum(sum_Ef2)
  sebeta = sqrt(1/(sum_Ef2))
  betahat = as.vector( Y %*% (Ef/sigmae2) ) / (sum_Ef2)
  betahat=as.vector(betahat)
  ATM = ash(betahat, sebeta, method="fdr", mixcompdist="normal")
  El = ATM$PosteriorMean
  SDl = ATM$PosteriorSD
  El2 = SDl^2 + El^2
  return(list(El = El, El2 = El2))
}


# Fval = function()







#' Factor Loading Adaptive Shrinkage (non-constant on columns)
#'
#' flash provide rank one matrix decomposition
#'
#' @param Y is the data matrix (N by P)
#' @param tol which is the stop criterion for the convergence, default is 1e-5
#' @param numtau number of iteration, default is 500. for the backfitting case, the number of tau should be 5 or 10.
#'
#' @details flash_r1 privide rank one matrix decomposition with variational EM algorithm.
#'
#' @export flash_r1c
#'
#' @importFrom ashr ash
#'
#' @return list of factor, loading and variance of noise matrix
#'  \itemize{
#'   \item{\code{l}} {is a N vector for loadings}
#'   \item{\code{f}} {is a P vector for factors}
#'   \item{\code{sigmae2}} {is mean of sigma square which is estimation for the noise variance}
#'  }
#' @examples
#' sim_hd = function(N, P, SF, SL, signal, a = rchisq(N,3),b = rchisq(P,1),mu = 0){
#' E = matrix(rep(0,N*P),nrow=N)
#' sig2_true = matrix(rep(0,N*P),nrow=N)
#' for(i in 1:N){
#'   for(j in 1:P){
#'     sig2_true[i,j] = mu + a[i] + b[j]
#'     E[i,j] = rnorm(1,0,sqrt(mu + a[i] + b[j]))
#'   }
#' }
#'
#' K=1
#' lstart = rnorm(N, 0, signal)
#'
#' fstart = rnorm(P, 0, signal)
#'
#' index = sample(seq(1:N),(N*SL))
#' lstart[index] = 0
#' index = sample(seq(1:P),(P*SF))
#' fstart[index] = 0
#'
#' Y = lstart %*% t(fstart) + E
#'
#' return(list(Y = Y, L_true = lstart, F_true = fstart, Error = E,sig2_true = sig2_true))
#' }
#' N = 200
#' P = 500
#' SF = 0.5
#' SL = 0.5
#' signal = 1
#' data = sim_hd(N, P, SF, SL, signal, a = rep(1,N),b = rchisq(P,2),mu = 0)
#' sigmae2_true = data$sig2_true
#' Y = data$Y
#' E = data$Error
#' gc = flash_r1c(Y,tol = 1e-04,maxiter_r1 = 300)
#' f = gc$f
#' l = gc$l
#' sqrt(mean(((Y - l %*% t(f))-E)^2))/sqrt(mean((Y - E)^2))
#' # check the colmean
#' col_mean = gc$sigmae2
#' col_true = as.vector(sigmae2_true[1,])
#' plot(col_true,col_mean)
#' abline(a = 0,b=1,col = "red")
#' gc = flash_r1(Y,tol = 1e-04,maxiter_r1 = 300)
#' f = gc$f
#' l = gc$l
#' sqrt(mean(((Y - l %*% t(f))-E)^2))/sqrt(mean((Y - E)^2))


# set the number of iteration as numtau
flash_r1c = function(Y, tol=1e-6, maxiter_r1 = 500){
  #dealing with missing value
  Y[is.na(Y)] = 0
  # get initial value for l and f and sigmae
  El = svd(Y)$u[,1]
  El2 = El^2
  Ef = as.vector(t(El)%*%Y)
  Ef2 = Ef^2

  #start iteration
  sigmae2_v = colMeans( Y^2 - 2*Y*(El %*% t(Ef)) + (El2 %*% t(Ef2)) )

  par_f = ATM_f_c(Y,El,El2,sigmae2_v)
  Ef = par_f$Ef
  Ef2 = par_f$Ef2

  sigmae2_v = colMeans( Y^2 - 2*Y*(El %*% t(Ef)) + (El2 %*% t(Ef2)) )
  #sigmae2
  par_l = ATM_l_c(Y,Ef,Ef2,sigmae2_v)
  El = par_l$El
  El2 = par_l$El2

  epsilon = 1
  tau = 1
  while(epsilon >= tol & tau < maxiter_r1 ){
    tau = tau + 1
    pre_sigmae2 = sigmae2_v

    sigmae2_v = colMeans( Y^2 - 2*Y*(El %*% t(Ef)) + (El2 %*% t(Ef2)) )

    par_f = ATM_f_c(Y,El,El2,sigmae2_v)
    Ef = par_f$Ef
    Ef2 = par_f$Ef2
    if(sum(Ef^2)==0){
      El = rep(0,length(El))
      Ef = rep(0,length(Ef))
      break
    }
    sigmae2_v = colMeans( Y^2 - 2*Y*(El %*% t(Ef)) + (El2 %*% t(Ef2)) )
    #sigmae2
    par_l = ATM_l_c(Y,Ef,Ef2,sigmae2_v)
    El = par_l$El
    El2 = par_l$El2
    if(sum(El^2)==0){
      El = rep(0,length(El))
      Ef = rep(0,length(Ef))
      break
    }
    epsilon = sum(abs(pre_sigmae2 - sigmae2_v ))
  }
  return(list(l = El, f = Ef, sigmae2 = sigmae2_v))
}
