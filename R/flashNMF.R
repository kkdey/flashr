sim_data = function(N, P, SF, SL, signal,
                    a = rchisq(N,3), b = rchisq(P,1),
                    mu = 0, K = 6, positive = FALSE){

  E = matrix(rep(0,N*P),nrow=N)
  sig2_true = matrix(rep(0,N*P),nrow=N)
  for(i in 1:N){
    for(j in 1:P){
      sig2_true[i,j] = mu + a[i] + b[j]
      E[i,j] = rnorm(1,0,sqrt(mu + a[i] + b[j]))
    }
  }

  Y = E
  L_true = array(0, dim = c(N,K))
  F_true = array(0, dim = c(P,K))

  for(k in 1:K){
    lstart = rnorm(N, 0, signal)
    fstart = rnorm(P, 0, signal)

    index = sample(seq(1:N),(N*SL))
    lstart[index] = 0
    index = sample(seq(1:P),(P*SF))
    fstart[index] = 0

    if(positive == TRUE){
      lstart = abs(lstart)
      fstart = abs(fstart)
    }

    L_true[,k] = lstart
    F_true[,k] = fstart

    Y = Y + lstart %*% t(fstart)
  }

  return(list(Y = Y, L_true = L_true, F_true = F_true, Error = E,sig2_true = sig2_true))
}
#data = sim_hd(N, P, SF, SL, signal, a = rchisq(N,3),b = rchisq(P,1),mu = 0)

# test the flash_hd
#data = sim(N, P, SF = SF, SL = SL, signal = signal)
#Y = data$Y
#ghd = flash_hd(Y,partype = "constant")
#gvem = flash(Y,tol=1e-6,numtau = 500)


flash_nmf_r1 = function(Y, tol=1e-3, maxiter_r1 = 50,
                     partype = c("constant","known","Bayes_var","var_col","noisy"),
                     sigmae2_true = NA,
                     factor_value = NA,fix_factor = FALSE,
                     nonnegative = FALSE,
                     objtype = c("margin_lik","lowerbound_lik"),
                     ash_para = list(),
                     fl_list=list(),
                     initial_nmf){
  partype = match.arg(partype, c("constant","known","Bayes_var","var_col","noisy"))
  objtype = match.arg(objtype, c("margin_lik","lowerbound_lik"))
  # to get the intial value.
  N = dim(Y)[1]
  P = dim(Y)[2]
  nmf_L = initial_nmf$L
  nmf_F = initial_nmf$F
  # initialization
  El = nmf_L
  Ef = nmf_F
  El2 = El^2
  Ef2 = Ef^2
  sigmae2_v = sigmae2_true

  g_update = one_step_update(Y, El, El2, Ef, Ef2,
                             N, P,
                             sigmae2_v, sigmae2_true,
                             sigmae2 = NA,
                             nonnegative ,
                             partype ,
                             objtype ,
                             fix_factor,
                             ash_para,
                             fl_list)
  # parameters updates
  El = g_update$El
  El2 = g_update$El2
  Ef = g_update$Ef
  Ef2 = g_update$Ef2
  sigmae2_v = g_update$sigmae2_v
  sigmae2 = g_update$sigmae2
  obj_val = g_update$obj_val


  obj_val_track = c(obj_val)

  # we should also return when the first run get all zeros
  if(sum(El^2)==0 || sum(Ef^2)==0){
    sigmae2 = sigma_est(sigmae2_v,sigmae2_true,sigmae2 ,partype)
    if(is.list(sigmae2)){
      # print("here using kronecker product")
      sigmae2 = sigmae2$sig2_l %*% t(sigmae2$sig2_f)
    }
    # add one more output for greedy algorithm which not useful here
    c_lik_val = C_likelihood(N,P,sigmae2_v,sigmae2)$c_lik
    # the above value is not useful, but is helpful to get the postprior value
    # since obj_val = c_lik_value + priorpost_l + priorpost_f
    return(list(l = El, f = Ef, l2 = El2, f2 = Ef2,
                sigmae2 = sigmae2,
                obj_val = obj_val,
                c_lik_val = c_lik_val))
  }

  epsilon = 1
  tau = 1
  while(epsilon >= tol & tau < maxiter_r1){
    tau = tau + 1
    pre_obj = obj_val

    g_update = one_step_update(Y, El, El2, Ef, Ef2,
                               N, P,
                               sigmae2_v, sigmae2_true,
                               sigmae2,
                               nonnegative ,
                               partype ,
                               objtype ,
                               fix_factor,
                               ash_para,
                               fl_list)
    # parameters updates
    El = g_update$El
    El2 = g_update$El2
    Ef = g_update$Ef
    Ef2 = g_update$Ef2
    sigmae2_v = g_update$sigmae2_v
    sigmae2 = g_update$sigmae2
    obj_val = g_update$obj_val

    if(sum(El^2)==0 || sum(Ef^2)==0){
      El = rep(0,length(El))
      Ef = rep(0,length(Ef))
      break
    }
    epsilon = abs(pre_obj - obj_val)
    obj_val_track = c(obj_val_track,obj_val)
  }
  sigmae2 = sigma_est(sigmae2_v,sigmae2_true,sigmae2 ,partype)
  if(is.list(sigmae2)){
    # print("here, using kronecker product")
    sigmae2 = sigmae2$sig2_l %*% t(sigmae2$sig2_f)
  }
  # add one more output for greedy algorithm which not useful here
  c_lik_val = C_likelihood(N,P,sigmae2_v,sigmae2)$c_lik
  # the above value is not useful, but is helpful to get the postprior value
  # since obj_val = c_lik_value + priorpost_l + priorpost_f
  return(list(l = El, f = Ef, l2 = El2, f2 = Ef2,
              sigmae2 = sigmae2,
              obj_val = obj_val,
              c_lik_val = c_lik_val,
              obj_val_track = obj_val_track))
}

# just for the case that the variance is known
flash_nmf = function(Y,initial_list, maxiter_bf=100,
                       flash_para = list(), gvalue = c("lik","eigen"),
                       parallel = FALSE){
  epsilon = 1
  tau = 1
  N = dim(Y)[1]
  P = dim(Y)[2]
  # match the input parameter
  gvalue = match.arg(gvalue, c("lik","eigen"))
  # set the default value for flash
  flash_default = list(tol=1e-3, maxiter_r1 = 50,
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
  while(epsilon> flash_para$tol & tau < maxiter_bf){
    tau = tau + 1
    K = dim(Lest)[2]
    if(K==0){break} #tests for case where all factors disappear!
    # this one can be put out of the while loop
    # preRMSfl = sqrt(mean((Lest %*% t(Fest))^2))
    pre_obj_lik = obj_lik
    sigmae2_out = list()
    for(k in 1:K){
      cat("Applying FLASH on factor:", k, "\n");

      residual = Y - Lest[,-k] %*% t(Fest[,-k])
   #   residual = residual * (residual > 0)
      flash_para$Y = residual
      flash_para$initial_nmf$L = Lest[,k]
      flash_para$initial_nmf$F = Fest[,k]
      # flash_para$fl_list$Ef2 = F2est[,-k]
      # flash_para$fl_list$El2 = L2est[,-k]
      # run the rank one flash
      r_flash = do.call(flash_nmf_r1,flash_para)
      Lest[,k] = r_flash$l
      Fest[,k] = r_flash$f
      L2est[,k] = r_flash$l2
      F2est[,k] = r_flash$f2
      sigmae2_out[[k]] = r_flash$sigmae2
      priorpost = r_flash$obj_val - r_flash$c_lik_val
      priorpost_vec[k] = priorpost
      clik = r_flash$c_lik_val
      obj_lik = clik + sum(priorpost_vec)
      track_obj = c(track_obj,obj_lik)
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

