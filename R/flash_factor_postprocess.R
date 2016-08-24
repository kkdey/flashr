
#' post-processing the loadings and factors from flash output
#'
#' @return sparsity and PVE and big genes
#' @export

flash_factor_postprocess = function(loadings, factors, data){
  PVE <- array(0, (dim(factors)[2]-1));
  sparsity <- array(0, dim(factors)[2])
  big_genes_2 <- array(0, dim(factors)[2])
  big_genes_5 <- array(0, dim(factors)[2])
  prop_positive <- array(0, dim(factors)[2])
  prop_negative <- array(0, dim(factors)[2])
  for(num in 1:dim(factors)[2]){
    max_gene <- max(factors[,num]);
    abs_fac <- abs(factors[,num]);
    abs_fac_max <- max(abs_fac);
    zero_indices <- which(((abs_fac)/abs_fac_max) < 1e-04);
    factors[zero_indices,num]=0;
    abs_fac[zero_indices]=0;
    big_genes_2[num] <- length(which(abs_fac > (0.2*abs_fac_max)))/(length(abs_fac));
    big_genes_5[num] <- length(which(abs_fac > (0.5*abs_fac_max)))/(length(abs_fac));
    prop_positive[num] <- length(which(factors[,num] > 0))/length(factors[,num]);
    prop_negative[num] <- length(which(factors[,num] < 0))/length(factors[,num]);
    sparsity[num] <- (length(which(abs_fac==0)))/(length(abs_fac))
    if(num==1){
      temp <- loadings[,num]%*%t(factors[,num]);
      res <- data - temp;
    }else{
      temp <- loadings[,num]%*%t(factors[,num]);
      PVE[num-1] <- (sum(temp^2))/(sum(res^2));
    }
  }

  post <- list("PVE"=PVE,
               "sparsity_prop"=sparsity,
               "prop_positive_features"=prop_positive,
               "prop_negative_features"=prop_negative,
               "big_genes_2"=big_genes_2,
               "big_genes_5"=big_genes_5)
  return(post)
}
