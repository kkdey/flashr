
# tautol is the number of iterations here
backfitting = function(Y,Lest,Fest,tautol=100,numtau = 500){
  # backfitting with initial values
  epsilon = 1
  tau = 1
  while(epsilon>1e-5 & tau < tautol){
    tau = tau + 1
    # if the l or f is a vector
    if(is.vector(Lest) || is.vector(Fest)){
      residual = Y - Lest %*% t(Fest)
      preRMSfl = sqrt(mean((Lest %*% t(Fest))^2))
      residual = residual + Lest %*% t(Fest)
      r_flash = flash_VEM(residual,numtau = numtau)
      Lest = r_flash$l
      Fest = r_flash$f
      residual = residual - Lest %*% t(Fest)
    }else{
      K = dim(Lest)[2]
      # this one can be put out of the while loop
      residual = Y - Lest %*% t(Fest)
      preRMSfl = sqrt(mean((Lest %*% t(Fest))^2))
      for(i in 1:K){
        residual = residual + Lest[,i] %*% t(Fest[,i])
        r_flash = flash_VEM(residual,numtau = numtau)
        Lest[,i] = r_flash$l
        Fest[,i] = r_flash$f
        residual = residual - Lest[,i] %*% t(Fest[,i])
      }
      # remove the zero in the l and f
      while(i <= dim(Lest)[2] ){
        if(sum((Lest[,i])^2)==0 || sum((Fest[,i])^2)==0){
          Lest = Lest[,-i]
          Fest = Fest[,-i]
        }
        numfactor = ifelse(is.vector(Lest),1,dim(Lest)[2])
        if(numfactor == 1){
          break
        }
        i = i+1
      }
    }

    RMSfl = sqrt(mean((Lest %*% t(Fest))^2))
    epsilon = abs(preRMSfl - RMSfl)
  }
  return(list(Lest = Lest,Fest = Fest))
}
