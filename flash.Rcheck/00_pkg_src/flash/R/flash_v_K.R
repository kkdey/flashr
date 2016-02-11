flash_v_K = function(Y,K){
  N = dim(Y)[1]
  P = dim(Y)[2]
  #initial the residual for the greedy algorithm
  residual = Y
  # do rank one decomposition
  r_flash = flash_VEM(residual)
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
    for(k in 1:K){
      #rank one for residual
      r_flash = flash_VEM(residual)
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
