

############  Deng et al TF - FLASH-tpx application  ######################

source("../R/flash.R")
source("../R/flashNMF.R")

read.data1 = function() {
  x = tempfile()
  download.file('https://cdn.rawgit.com/kkdey/singleCellRNASeqMouseDeng2014/master/data/Deng2014MouseEsc.rda', destfile=x, quiet=TRUE)
  z = get(load((x)))
  return(z)
}

Deng2014MouseESC <- read.data1()
deng.counts <- Biobase::exprs(Deng2014MouseESC)
deng.meta_data <- Biobase::pData(Deng2014MouseESC)
deng.gene_names <- rownames(deng.counts)

gene_names <- get(load(file="../../count-clustering/project/external_data/Deng_Data/TF_gene_names.rda"))

matched_indices <- match(gene_names, deng.gene_names)

deng_counts_TF <- t(deng.counts[matched_indices,])
deng_counts_TF <- deng_counts_TF[,-which(colSums(deng_counts_TF)==0)]

gnmf = NMF::nmf(deng_counts_TF,rank = 6)

projected <- gnmf@fit@W %*% gnmf@fit@H
data_temp <- pmin(deng_counts_TF, projected)
gnmf2 <- NMF::nmf(data_temp,rank = 6)

resid <- deng_counts_TF - gnmf2@fit@W%*% gnmf2@fit@H

initial_nmf_list = list(l = gnmf2@fit@W, f = t(gnmf2@fit@H),
                        l2 = (gnmf2@fit@W)^2, f2 = (t(gnmf2@fit@H))^2,
                        priorpost_vec = rep(1,6),clik_vec = rep(1,6))


gflash = flash_nmf(deng_counts_TF,initial_list = initial_nmf_list,
                   maxiter_bf=2,
                   flash_para = list(partype = "known",
                                     sigmae2_true = deng_counts_TF+0.001,
                                     nonnegative = TRUE),
                   gvalue = "eigen",parallel = FALSE)


flash_para = list(partype = "known",
                  sigmae2_true = deng_counts_TF+0.001,
                  nonnegative = TRUE)
gvalue = "eigen"
parallel = FALSE
maxiter_bf=2
initial_list = initial_nmf_list
Y = deng_counts_TF






topic_clus <- maptpx::topics(deng_counts_TF, K=6, tol=0.1)

lib_size <- rowSums(deng_counts_TF)

prob_adjusted <- topic_clus$omega%*%t(topic_clus$theta)

omega <- topic_clus$omega
theta <- topic_clus$theta

loadings <- matrix(0, dim(omega)[1], dim(omega)[2])

for(n in 1:dim(omega)[1]){
  loadings[n,] <- omega[n,]*lib_size[n];
}

factors <- theta

initial_nmf_list = list(l = loadings, f = factors,
                        l2 = loadings^2, f2 = factors^2,
                        priorpost_vec = rep(1,6),clik_vec = rep(1,6))

Y <- deng_counts_TF;
initial_list <- initial_nmf_list
res <- Y - initial_list$l%*%t(initial_list$f);

gflash = flash_nmf(deng_counts_TF,initial_list = initial_nmf_list,
                   maxiter_bf=2,
                   flash_para = list(partype = "constant", nonnegative = TRUE),
                   gvalue = "eigen",parallel = FALSE)




initial_nmf_list = list(l = gnmf@fit@W, f = t(gnmf@fit@H),
                        l2 = (gnmf@fit@W)^2, f2 = t((gnmf@fit@H))^2,
                        priorpost_vec = rep(1,6),clik_vec = rep(1,6))

Y <- deng_counts_TF;
initial_list <- initial_nmf_list
maxiter_bf <- 2
flash_para = list(partype = "constant", nonnegative = TRUE)
gvalue <- "eigen"
parallel <- FALSE

res <- Y - gnmf@fit@W[,-1] %*% gnmf@fit@H[-1,]

gflash = flash_nmf(deng_counts_TF,initial_list = initial_nmf_list,
                   maxiter_bf=2,
                   flash_para = list(partype = "constant", nonnegative = TRUE),
                   gvalue = "eigen",parallel = FALSE)

