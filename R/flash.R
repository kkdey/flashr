# ATM function
#' title ash type model for l
#'
#' description use ash type model to maxmization
#'
#' @return list of factor, loading and variance of noise matrix
#'  \itemize{
#'   \item{\code{El}} {is a N vector for mean of loadings}
#'   \item{\code{El2}} {is a N vector for second moment of loadings}
#'  }
#' @param Y is data matrix
#' @param Ef is mean for the factor
#' @param Ef2 is second moment for the factor
#' @param sigmae2 is variance structure for noise matrix
#' @param col_var is just for column variance, column is for column, row is for row
#' @param nonnegative is flag for whether output nonnegative value
#' @param output is format of output, "mean" is for mean value, "matrix" is for flash data matrix in ash
#'
#' @keywords internal
#'

ATM_r1 = function(Y, Ef, Ef2, 
                  sigmae2, col_var = "row", 
                  nonnegative = FALSE, output = "mean",
                  partype = "constant"){
  if(is.matrix(sigmae2)){
    # this is for all variance are different
    sum_Ef2 = (1/sigmae2) %*% Ef2
    sum_Ef2 = as.vector(sum_Ef2)
    sebeta = sqrt(1/(sum_Ef2))
    betahat = as.vector( (Y/sigmae2) %*% Ef ) / (sum_Ef2)
    betahat=as.vector(betahat)
  } else if(is.vector(sigmae2) & length(sigmae2) == length(Ef)){
    # this is for the non-constant variance in column
    if(col_var == "row"){
      sum_Ef2 = (1/sigmae2) * Ef2
      sum_Ef2 = sum(sum_Ef2)
      sebeta = sqrt(1/(sum_Ef2))
      betahat = as.vector( Y %*% (Ef/sigmae2) ) / (sum_Ef2)
      betahat=as.vector(betahat)
    }else {
      # for column
      # here I have alread change the dimension by taking transpose to Y
      sum_Ef2 = sum(Ef2)
      sebeta = sqrt( sigmae2/(sum_Ef2) )
      betahat = (t(Ef) %*% t(Y)) / (sum_Ef2)
      betahat=as.vector(betahat)
    }
  }else{
    # for constant case
    sum_Ef2 = sum(Ef2)
    sebeta = sqrt(sigmae2/(sum_Ef2))
    betahat = (t(Ef) %*% t(Y)) / (sum_Ef2)
    betahat=as.vector(betahat)
  }
  # ATM update
  # decide the sign for output
  if(nonnegative){
    mixdist = "+uniform"
  }else {
    mixdist = "normal"
  }
  # decide the output to decide the convergence criterion
  if(output == "matrix"){
    ATM = ash(betahat, sebeta, method="fdr", mixcompdist=mixdist,outputlevel=4)
    Ef = ATM$PosteriorMean
    SDf = ATM$PosteriorSD
    Ef2 = SDf^2 + Ef^2
    mat_post = ATM$flash.data
    fit_g = ATM$fitted.g
    return(list(Ef = Ef,
                Ef2 = Ef2,
                mat = mat_post,
                g = fit_g))
  } else {
    ATM = ash(betahat, sebeta, method="fdr", mixcompdist=mixdist)
    Ef = ATM$PosteriorMean
    SDf = ATM$PosteriorSD
    Ef2 = SDf^2 + Ef^2
    return(list(Ef = Ef, Ef2 = Ef2))
  }
  
}

#' title prior and posterior part in objective function
#'
#' description prior and posterior part in objective function
#'
#' @return PrioPost the value for the proir and posterior in objectice function
#' @param mat matrix of flash.data in ash output which is for posterior
#' @param fit_g in fitted.g in ash output which is for prior
#' @keywords internal
#'
# two parts for the likelihood
Fval = function(mat,fit_g){
  prior_pi = fit_g$pi
  nonzeroindex = which(prior_pi!=0)
  prior_var = (fit_g$sd[nonzeroindex])^2
  prior_pi = prior_pi[nonzeroindex]
  
  mat_postmean = mat$comp_postmean[nonzeroindex,]
  mat_postmean2 = mat$comp_postmean2[nonzeroindex,]
  mat_postprob = mat$comp_postprob[nonzeroindex,]
  mat_postvar = mat_postmean2 - mat_postmean^2
  dimension = dim(mat_postmean)
  K = dimension[1]
  N = dimension[2]
  mat_priorprob = matrix(rep(prior_pi,N),ncol = N)
  mat_priorvar = matrix(rep(prior_var,N),ncol = N)
  probodds = log(mat_priorprob/mat_postprob)
  #log(mat_priorprob) - log(mat_postprob)
  varodds = log(mat_priorvar/mat_postvar)
  varodds[1,] = 0
  likodds = mat_postmean2/mat_priorvar
  likodds[1,] = 0
  # here mat_postprob might be equal to zero, so 0 * inf should be equal to zero
  PrioPost = mat_postprob*(probodds - (1/2)*varodds - (1/2)*likodds )
  # PrioPost[which(mat_postprob ==0 )] = 0
  PrioPost[which(mat_postprob< 1e-100 )] = 0
  return(list(PrioPost = sum(PrioPost) ))
}

#' title conditional likelihood
#'
#' description conditional likelihood in objective function
#'
#' @return c_lik for the onditional likelihood in objectice function
#' @param N  dimension of residual matrix 
#' @param P  dimension of residual matrix
#' @param sigmae2_v  residual matrix
#' @param sigmae2_true  true variance structure (we use the estimated one to replace that if the truth is unknown)
#' @keywords internal
#'
# this version we need to know the truth of sigmae2, we can use sigmae2_v as the truth
C_likelihood = function(N,P,sigmae2_v,sigmae2_true){
  if(is.matrix(sigmae2_true)){
    c_lik = -(1/2) * sum( log(2*pi*sigmae2_true) + (sigmae2_v)/(sigmae2_true) )
  }else if(is.vector(sigmae2_true)){
    # change the format to fit the conditional likelihood
    sigmae2_v = colMeans(sigmae2_v)
    c_lik = -(N/2) * sum( log(2*pi*sigmae2_true) + (sigmae2_v)/(sigmae2_true) )
  } else {
    # change the format to fit the conditional likelihood and accelerate the computation.
    sigmae2_v = mean(sigmae2_v)
    c_lik = -(N*P)/2 * ( log(2*pi*sigmae2_true) + (sigmae2_v)/(sigmae2_true) )
  }
  return(list(c_lik = c_lik))
}

#' title objective function in VEM
#'
#' description  objective function
#'
#' @return obj_val value of the objectice function
#' @param N  dimension of residual matrix 
#' @param P  dimension of residual matrix
#' @param sigmae2_v  residual matrix
#' @param sigmae2_true  true variance structure (we use the estimated one to replace that if the truth is unknown)
#' @param par_l ash output for l
#' @param par_f ash output for f
#' @param objtype  objective function type, 
#' "margin_lik" for conditional likelihood,
#' "lowerbound_lik" for full objective function
#' @keywords internal
#'

# objective function
obj = function(N,P,sigmae2_v,sigmae2_true,par_f,par_l,objtype = "margin_lik"){
  if(is.list(sigmae2_true)){
    sigmae2_true = sigmae2_true$sig2_l %*% t(sigmae2_true$sig2_f)
  }
  if(objtype=="lowerbound_lik"){
    priopost_f = Fval(par_f$mat, par_f$g)$PrioPost
    priopost_l = Fval(par_l$mat, par_l$g)$PrioPost
    c_lik = C_likelihood(N,P,sigmae2_v,sigmae2_true)$c_lik
    obj_val = c_lik + priopost_l + priopost_f
  } else {
    obj_val = C_likelihood(N,P,sigmae2_v,sigmae2_true)$c_lik
  }
  return(obj_val)
}

#' title rescale the lambda_l and lambda_f for identifiablity
#'
#' description rescale the lambda_l and lambda_f for identifiablity
#'
#' @return a list of sig2_l and sig2_f for the rescaled variance of the kronecker product
#' @param sig2_l variance of the kronecker product
#' @param sig2_l variance of the kronecker product
#' @keywords internal
#' 
rescale_sigmae2_true = function(sig2_l,sig2_f){
  norm_l = sqrt(sum(sig2_l^2))
  norm_f = sqrt(sum(sig2_f^2))
  norm_total = norm_l * norm_f
  sig2_l = sig2_l / norm_l
  sig2_f = sig2_f / norm_f
  sig2_l = sig2_l * sqrt(norm_total)
  sig2_f = sig2_f * sqrt(norm_total)
  return(list(sig2_l = sig2_l,sig2_f = sig2_f))
}  

#' title initial value for Bayes variance structure estimation for kronecker productor
#'
#' description initial value for Bayes variance structure estimation for kronecker productor
#'
#' @return sig2_out estimated variance structure
#' @param sigmae2_v  residual matrix
#' @param sigmae2_true  true variance structure (we use the estimated one to replace that if the truth is unknown)
#' @keywords internal
#' 
inital_Bayes_var = function(sigmae2_v){
  N = dim(sigmae2_v)[1]
  P = dim(sigmae2_v)[2]
  # estimate the initial value of lambda_l and lambda_f
  sig2_l_pre = rowMeans(sigmae2_v)
  sig2_f_pre = colMeans( sigmae2_v / matrix(rep(sig2_l_pre,P), ncol = P) )
  sig2_pre_list = rescale_sigmae2_true(sig2_l_pre,sig2_f_pre)
  sig2_l_pre = sig2_pre_list$sig2_l
  sig2_f_pre = sig2_pre_list$sig2_f
  # start the iteration
  maxiter = 100
  inital_tol = 1e-3
  tau = 0
  epsilon = 1
  while(epsilon > inital_tol & tau <= maxiter){
    tau = tau + 1
    sig2_l = rowMeans( sigmae2_v / matrix(rep(sig2_f_pre,each = N),ncol = P) )
    sig2_f = colMeans( sigmae2_v / matrix(rep(sig2_l,P), ncol = P) )
    sig2_list = rescale_sigmae2_true(sig2_l,sig2_f)
    sig2_l = sig2_list$sig2_l
    sig2_f = sig2_list$sig2_f
    epsilon = sqrt(mean((sig2_f - sig2_f_pre)^2 )) + sqrt(mean((sig2_l - sig2_l_pre)^2))
    sig2_l_pre = sig2_l
    sig2_f_pre = sig2_f
  }
  return(list(sig2_l = sig2_l,sig2_f = sig2_f))
}

#' title Bayes variance structure estimation for kronecker productor
#'
#' description prior and posterior part in objective function
#'
#' @return sig2_out estimated variance structure
#' @param sigmae2_v  residual matrix
#' @param sigmae2_true  true variance structure (we use the estimated one to replace that if the truth is unknown)
#' and it is a list of two vectors in this case.
#' @keywords internal
#'
Bayes_var = function(sigmae2_v,sigmae2_true){
  N = dim(sigmae2_v)[1]
  P = dim(sigmae2_v)[2]
  if( is.na(sigmae2_true) || !is.list(sigmae2_true) ){
    # we don't know the truth this is in the first iteration
    sigmae2_true = inital_Bayes_var(sigmae2_v)
  }
  sig2_l_pre = sigmae2_true$sig2_l
  sig2_f_pre = sigmae2_true$sig2_f
  # this has already beed rescaled
  # here we use alpha_l = alpha_f = beta_l = beta_f = 0
  sig2_l = rowMeans( sigmae2_v / matrix(rep(sig2_f_pre,each = N),ncol = P) )
  sig2_f = colMeans( sigmae2_v / matrix(rep(sig2_l,P), ncol = P) )
  #rescaled the variance
  sig2_list = rescale_sigmae2_true(sig2_l,sig2_f)
  sig2_l = sig2_list$sig2_l
  sig2_f = sig2_list$sig2_f
  # sig2_out = matrix(rep(sig2_l,P),ncol = P) * matrix(rep(sig2_f,each = N),ncol = P)
  return(list(sig2_l = sig2_l,sig2_f = sig2_f))
}
  
#' title module for estiamtion of the variance structure
#'
#' description estiamtion of the variance structure
#'
#' @return sigmae2 estimated variance structure
#' @param partype parameter type for the variance, 
#' "constant" for constant variance, 
#' "var_col" for nonconstant variance for column, 
#' "known" for the kown variance,
#' "Bayes_var" for Bayes version of the nonconstant variance for row and column
#' "loganova" is anova estiamtion for the log residual square
#' @param sigmae2_v  residual matrix
#' @param sigmae2_true  true variance structure (we use the estimated one to replace that if the truth is unknown)
#' @keywords internal
#'
# sigma estimation function
sigma_est = function(sigmae2_v,sigmae2_true,partype = "constant"){
  if(partype == "var_col"){
    sigmae2 = colMeans(sigmae2_v)
  } else if(partype == "loganova"){
    sigmae2 = log_anova(sigmae2_v)
  } else if (partype == "Bayes_var"){
    # here sigmae2_true is a list
    sigmae2 = Bayes_var(sigmae2_v,sigmae2_true)
  }else if (partype == "known"){
    sigmae2 = sigmae2_true
  } else {
    # this is for constant case
    sigmae2 = mean(sigmae2_v)
  }
  return(sigmae2)
}
  
#' title one step update in flash iteration using ash
#'
#' description one step update in flash iteration using ash
#' @return list of factor, loading and variance of noise matrix
#'  \itemize{
#'   \item{\code{El}} {is a N vector for mean of loadings}
#'   \item{\code{El2}} {is a N vector for second moment of loadings}
#'   \item{\code{Ef}} {is a N vector for mean of factors}
#'   \item{\code{Ef2}} {is a N vector for second moment of factors}
#'   \item{\code{sigmae2_v}}{is a N by P matrix for residual square}
#'   \item{\code{sigmae2_true}}{is a N by P matrix for estimated value for the variance structure}
#'   \item{\code{obj_val}}{the value of objectice function}
#'  }
#' @param Y the data matrix
#' @param N dimension of Y
#' @param P dimension of Y
#' @param El mean for the loadings
#' @param El2 second moment for the loadings
#' @param Ef mean for the factors
#' @param Ef2 second moment for the factors
#' @param sigmae2_v residual square
#' @param sigmae2_true estimated value for the variance structure
#' @param nonnegative if the facotor and loading are nonnegative or not. 
#' TRUE for nonnegative
#' FALSE for no constraint
#' @param partype parameter type for the variance, 
#' "constant" for constant variance, 
#' "var_col" for nonconstant variance for column, 
#' "known" for the kown variance,
#' "Bayes_var" for Bayes version of the nonconstant variance for row and column
#' "loganova" is anova estiamtion for the log residual square
#' @param objtype  objective function type, 
#' "margin_lik" for conditional likelihood,
#' "lowerbound_lik" for full objective function
#' @param fix_factor whether the factor is fixed or not
#' TRUE for fix_factor
#' FALSE for non-constraint 
#' @keywords internal
#'
# one step update function 
one_step_update = function(Y, El, El2, Ef, Ef2,
                           N, P,
                           sigmae2_v, sigmae2_true,
                           nonnegative = FALSE,
                           partype = "constant",
                           objtype = "margin_lik",
                           fix_factor = FALSE){
  # if fix_factor is True, please choose objtype = "margin_lik"
  output = ifelse(objtype == "lowerbound_lik", "matrix", "mean")
  # deal with the missing value
  na_index_Y = is.na(Y)
  is_missing = any(na_index_Y)  # keep the missing result
  if(is_missing){
    Y[na_index_Y] = (El %*% t(Ef))[na_index_Y]
  }
  if(fix_factor){
    # actually we need do nothing here
    output = "mean"
    objtype = "margin_lik"
  }else{
    sigmae2 = sigma_est(sigmae2_v,sigmae2_true,partype)
    # for the list case in kronecker product
    if(is.list(sigmae2)){
      sigmae2 = sigmae2$sig2_l %*% t(sigmae2$sig2_f)
    }
    # Y = lf^T + E and ATM is for l given f, for f given l we need y^T = fl^T + E^T
    # sigmae2_input = ifelse(is.matrix(sigmae2),t(sigmae2),sigmae2) is wrong here
    if(is.matrix(sigmae2)){
      sigmae2_input = t(sigmae2)
    }else{
      sigmae2_input = sigmae2
    }
    par_f = ATM_r1(t(Y), El, El2, sigmae2_input,
                   col_var = "column", nonnegative,
                   output,partype)
    Ef = par_f$Ef
    Ef2 = par_f$Ef2
    # if the Ef is zeros ,just return zeros
    if(sum(Ef^2)==0){
      El = rep(0,length(El))
      Ef = rep(0,length(Ef))
      sigmae2_v = Y^2
      sigmae2_true = sigma_est(sigmae2_v,sigmae2_true,partype)
      obj_val = obj(N, P, sigmae2_v, sigmae2_true, par_f=NA, par_l=NA, objtype="margin_lik")
      return(list(El = El, El2 = El^2,
                  Ef = Ef, Ef2 = Ef^2,
                  sigmae2_v = sigmae2_v,
                  sigmae2_true = sigmae2_true,
                  obj_val = obj_val))
    }
    # update the Y and the Y^2 in this if else block
    if(is_missing){
      ElEf = (El %*% t(Ef))
      El2Ef2 = (El2 %*% t(Ef2))
      Y[na_index_Y] = ElEf[na_index_Y]
      Y2 = Y^2
      if(is.matrix(sigmae2)){
        sigmae2_impute = sigmae2[na_index_Y]
      }else if(is.vector(sigmae2) & length(sigmae2) == length(Ef)){
        sigmae2_impute = (matrix(rep(sigmae2,each = N),ncol = P))[na_index_Y]
      }else{
        sigmae2_impute = sigmae2
      }
      Y2[na_index_Y] = El2Ef2[na_index_Y] + sigmae2_impute
      sigmae2_v =  Y2 - 2*Y*ElEf + El2Ef2
    }else{
      sigmae2_v =  Y^2 - 2*Y*(El %*% t(Ef)) + (El2 %*% t(Ef2))
    }
  }
  
  sigmae2 = sigma_est(sigmae2_v,sigmae2_true,partype)
  if(is.list(sigmae2)){
    # kronecker product
    sigmae2 = sigmae2$sig2_l %*% t(sigmae2$sig2_f)
  }
  par_l = ATM_r1(Y, Ef, Ef2, 
                 sigmae2, col_var = "row",
                 nonnegative, output,partype)
  El = par_l$Ef
  El2 = par_l$Ef2
  # if El is zeros just return
  if(sum(El^2)==0){
    El = rep(0,length(El))
    Ef = rep(0,length(Ef))
    sigmae2_v = Y^2
    sigmae2_true = sigma_est(sigmae2_v,sigmae2_true,partype)
    obj_val = obj(N, P, sigmae2_v, sigmae2_true, par_f=NA, par_l=NA, objtype="margin_lik")
    return(list(El = El, El2 = El^2,
                Ef = Ef, Ef2 = Ef^2,
                sigmae2_v = sigmae2_v,
                sigmae2_true = sigmae2_true,
                obj_val = obj_val))
  }
  if(is_missing){
    ElEf = (El %*% t(Ef))
    El2Ef2 = (El2 %*% t(Ef2))
    Y[na_index_Y] = ElEf[na_index_Y]
    Y2 = Y^2
    if(is.matrix(sigmae2)){
      sigmae2_impute = sigmae2[na_index_Y]
    }else if(is.vector(sigmae2) & length(sigmae2) == length(Ef)){
      sigmae2_impute = (matrix(rep(sigmae2,each = N),ncol = P))[na_index_Y]
    }else{
      sigmae2_impute = sigmae2
    }
    Y2[na_index_Y] = El2Ef2[na_index_Y] + sigmae2_impute
    sigmae2_v =  Y2 - 2*Y*ElEf + El2Ef2
  }else{
    sigmae2_v =  Y^2 - 2*Y*(El %*% t(Ef)) + (El2 %*% t(Ef2))
  }
  #use the estiamtion as the truth
  sigmae2_true = sigma_est(sigmae2_v,sigmae2_true,partype)
  obj_val = obj(N, P, sigmae2_v, sigmae2_true, par_f, par_l, objtype)
  
  return(list(El = El, El2 = El2,
              Ef = Ef, Ef2 = Ef2,
              sigmae2_v = sigmae2_v,
              sigmae2_true = sigmae2_true,
              obj_val = obj_val))
}

#' inital value for flash
#' 
#' description inital value for flash
#' 
#' @return list of factor, loading and variance of noise matrix
#'  \itemize{
#'   \item{\code{El}} {is a N vector for mean of loadings}
#'   \item{\code{El2}} {is a N vector for second moment of loadings}
#'   \item{\code{Ef}} {is a N vector for mean of factors}
#'   \item{\code{Ef2}} {is a N vector for second moment of factors}
#'   \item{\code{sigmae2_v}}{is a N by P matrix for residual square}
#'   \item{\code{sigmae2_true}}{is a N by P matrix for estimated value for the variance structure}
#'  }
#' @param Y the data matrix
#' @param nonnegative if the facotor and loading are nonnegative or not. 
#' TRUE for nonnegative
#' FALSE for no constraint
#' @param fix_factor whether the factor is fixed or not
#' TRUE for fix_factor
#' FALSE for non-constraint 
#' @param factor_value is the factor value if the factor is fixed
#' @keywords internal
#'
initial_value = function(Y, nonnegative = FALSE,
                         factor_value = NA,fix_factor = FALSE){
  # use the total mean as the estimated missing value
  na_index_Y = is.na(Y)
  is_missing = any(na_index_Y)
  if(is_missing){
    Y[na_index_Y] = mean(Y, na.rm = TRUE)
  }
  # the flash with plug-in missing value
  if(fix_factor){
    Ef = factor_value
    Ef2 = Ef^2
    El = as.vector( (Y %*% Ef) / (sum(Ef2)) )
    El2 = El^2
  }else{
    # El = svd(Y)$u[,1]
    El = as.vector(irlba::irlba(Y,nv = 0,nu = 1)$u)
    # the nonnegative value need positive inital value
    if(nonnegative){
      El = abs(El)
    }
    El2 = El^2
    Ef = as.vector(t(El)%*%Y)
    Ef2 = Ef^2
  }
  # residual matrix initialization
  # here we don't have any information about the sigmae2_true
  if(is_missing){
    ElEf = (El %*% t(Ef))
    El2Ef2 = (El2 %*% t(Ef2))
    Y[na_index_Y] = ElEf[na_index_Y]
    Y2 = Y^2
    Y2[na_index_Y] = El2Ef2[na_index_Y]
    sigmae2_v =  Y2 - 2*Y*ElEf + El2Ef2
  }else{
    sigmae2_v =  Y^2 - 2*Y*(El %*% t(Ef)) + (El2 %*% t(Ef2))
  }
  return(list(El = El, El2 = El^2,
              Ef = Ef, Ef2 = Ef^2,
              sigmae2_v = sigmae2_v))
}

#' FLASH
#'
#' factor loading adaptive shrinkage rank one version
#' @return list of factor, loading and variance of noise matrix
#'  \itemize{
#'   \item{\code{El}} {is a N vector for mean of loadings}
#'   \item{\code{Ef}} {is a N vector for mean of factors}
#'   \item{\code{sigmae2}}{is a N by P matrix for estimated value for the variance structure}
#'  }
#' @param Y the data matrix
#' @param tol is for the tolerence for convergence in iterations and ash
#' @param maciter_r1 is maximum of the iteration times for rank one case
#' @param sigmae2_true true value for the variance structure
#' @param nonnegative if the facotor and loading are nonnegative or not. 
#' TRUE for nonnegative
#' FALSE for no constraint
#' @param partype parameter type for the variance, 
#' "constant" for constant variance, 
#' "var_col" for nonconstant variance for column, 
#' "known" for the kown variance,
#' "Bayes_var" for Bayes version of the nonconstant variance for row and column
#' "loganova" is anova estiamtion for the log residual square
#' @param objtype  objective function type, 
#' "margin_lik" for conditional likelihood,
#' "lowerbound_lik" for full objective function
#' @param fix_factor whether the factor is fixed or not
#' TRUE for fix_factor
#' FALSE for non-constraint 
#' @param factor_value is the factor value if the factor is fixed
#'
#' @details flash privide rank one matrix decomposition with variational EM algorithm.
#'
#' @export flash
#'
#' @importFrom ashr ash
#'
flash = function(Y, tol=1e-5, maxiter_r1 = 500,
                 partype = "constant", sigmae2_true = NA, 
                 factor_value = NA,fix_factor = FALSE,
                 nonnegative = FALSE, objtype = "margin_lik" ){
  N = dim(Y)[1]
  P = dim(Y)[2]
  
  # to get the inital values
  g_inital = initial_value(Y, nonnegative, factor_value, fix_factor)
  El = g_inital$El
  Ef = g_inital$Ef
  El2 = g_inital$El2
  Ef2 = g_inital$Ef2
  sigmae2_v = g_inital$sigmae2_v
  
  # start iteration 
  g_update = one_step_update(Y, El, El2, Ef, Ef2,
                             N, P,
                             sigmae2_v, sigmae2_true,
                             nonnegative ,
                             partype ,
                             objtype ,
                             fix_factor)
  # parameters updates
  El = g_update$El
  El2 = g_update$El2
  Ef = g_update$Ef
  Ef2 = g_update$Ef2
  sigmae2_v = g_update$sigmae2_v
  sigmae2_true = g_update$sigmae2_true
  obj_val = g_update$obj_val
  
  epsilon = 1
  tau = 1
  while(epsilon >= tol & tau < maxiter_r1){
    tau = tau + 1
    pre_obj = obj_val
    
    g_update = one_step_update(Y, El, El2, Ef, Ef2,
                               N, P,
                               sigmae2_v, sigmae2_true,
                               nonnegative ,
                               partype ,
                               objtype ,
                               fix_factor)
    # parameters updates
    El = g_update$El
    El2 = g_update$El2
    Ef = g_update$Ef
    Ef2 = g_update$Ef2
    sigmae2_v = g_update$sigmae2_v
    sigmae2_true = g_update$sigmae2_true
    obj_val = g_update$obj_val
    
    if(sum(El^2)==0 || sum(Ef^2)==0){
      El = rep(0,length(El))
      Ef = rep(0,length(Ef))
      break
    }
    epsilon = abs(pre_obj - obj_val)
    print(obj_val)
  }
  sigmae2 = sigma_est(sigmae2_v,sigmae2_true,partype)
  return(list(l = El, f = Ef, sigmae2 = sigmae2))
}