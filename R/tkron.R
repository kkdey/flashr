
#' Tensor FLASH assuming a diagonal Kronecker structured covariance model.
#'
#' @inheritParams tflash
#'
#' @export
#'
#' @author David Gerard
tflash_kron <- function(Y, tol = 10^-5, itermax = 100, alpha = 0, beta = 0,
                        mixcompdist = "normal", nullweight = 10, print_update = FALSE,
                        start = c("first_sv", "random"), known_modes = NULL,
                        known_factors = NULL, homo_modes = NULL) {
    p <- dim(Y)
    n <- length(p)

    if (is.null(known_modes)) {
        unknown_modes <- 1:n
    } else {
        unknown_modes <- (1:n)[-known_modes]
    }

    start <- match.arg(start, c("first_sv", "random"))

    which_na <- is.na(Y)
    if (all(!which_na)) {
        which_na <- NULL
    }

    if (is.null(homo_modes)) {
        hetero_modes <- 1:n
    } else {
        hetero_modes <- (1:n)[-homo_modes]
    }

    init_return <- tinit_kron_components(Y = Y, which_na = which_na, start = start,
                                         known_factors = known_factors,
                                         known_modes = known_modes,
                                         homo_modes = homo_modes)
    ex_list <- init_return$ex_list # list of expected value of components.
    ex2_list <- init_return$ex2_list # list of expected value of x^2
    esig_list <- init_return$esig_list
    
    post_rate <- list()
    post_shape <- list()
    for(mode_index in 1:n) {
        post_shape[[mode_index]] <- rep(prod(p[-mode_index])) / 2 + alpha
        post_rate[[mode_index]] <- post_shape[[mode_index]] / esig_list[[mode_index]]
    }


    prob_zero <- list()
    pi0vec <- rep(NA, length = n)
    
    iter_index <- 1
    err <- tol + 1
    not_all_zero <- TRUE

    while (iter_index <= itermax & err > tol) {
        esig_list_old <- esig_list

        for (mode_index in 1:n) {
            ## update mean
            if (mode_index %in% unknown_modes) {
                if(sum(abs(ex_list[[mode_index]])) < 10^-6) {
                    ex_list <- lapply(ex_list, FUN = function(x) { rep(0, length = length(x)) })
                    not_all_zero <- FALSE
                    break
                }
                
                t_out <- tupdate_kron_modek(Y = Y, ex_list = ex_list, ex2_list = ex2_list,
                                            esig_list = esig_list, k = mode_index,
                                            mixcompdist = mixcompdist, which_na = which_na,
                                            nullweight = nullweight)
                ex_list <- t_out$ex_list
                ex2_list <- t_out$ex2_list
                prob_zero[[mode_index]] <- t_out$prob_zero
                pi0vec[mode_index] <- t_out$pi0
                
                if(sum(abs(ex_list[[mode_index]])) < 10^-6) {
                    ex_list <- lapply(ex_list, FUN = function(x) { rep(0, length = length(x)) })
                    not_all_zero <- FALSE
                    break
                }
            }
            
            ## update variance
            if (mode_index %in% hetero_modes) {
                post_rate[[mode_index]] <- tupdate_kron_sig(Y = Y, ex_list = ex_list,
                                                            ex2_list = ex2_list,
                                                            esig_list = esig_list,
                                                            k = mode_index, beta = beta,
                                                            which_na = which_na)
                esig_list[[mode_index]] <- post_shape[[mode_index]] / post_rate[[mode_index]]
            }
        }
        
        ## if(not_all_zero) {
        ##     max_ex <- sapply(ex_list, function(x) { max(abs(x)) })
        ##     if (max(outer(max_ex, max_ex, "/")) > 10^4)
        ##     {
        ##         ex_list <- rescale_factors(ex_list)
        ##         ex2_list <- lapply(ex_list, function(x) { x ^ 2 })
        ##     }
        ## }
        ## esig_list <- rescale_factors(esig_list)
        
        iter_index <- iter_index + 1

        ## calculate stopping criterion
        err <- 0
        for (sig_mode_index in hetero_modes) {
            old_scaled <- esig_list_old[[sig_mode_index]] /
                sqrt(sum(esig_list_old[[sig_mode_index]] ^ 2))
            new_scaled <- esig_list[[sig_mode_index]] /
                sqrt(sum(esig_list[[sig_mode_index]] ^ 2))
            err <- err + sum(abs(old_scaled - new_scaled))
        }

        if (print_update & iter_index %% 5 == 0) {
            cat("Iteration:", iter_index, "\n")
            cat("Stop Crit:", err, "\n\n")
        }
    }
    return(list(post_mean = ex_list, sigma_est = esig_list, prob_zero = prob_zero,
                pi0vec = pi0vec,
                num_iter = iter_index))
}

#' Update the variational density for the kth mode when assuming a
#' diagonal Kronecker structured covariance matrix.
#'
#' @inheritParams tupdate_modek
#' @param esig_list A list of vectors of positive numerics. The
#'     expected variances.
#' @param ex2_list A list of of vectors of positive numerics. The
#'     variational second moments of the components.
#'
#' @author David Gerard
#'
#' @export
#'
tupdate_kron_modek <- function(Y, ex_list, ex2_list, esig_list, k,
                               mixcompdist = "normal",
                               which_na = NULL, nullweight = 10) {
    p <- dim(Y)

    if (!is.null(which_na)) {
        Y[which_na] <- form_outer(ex_list)[which_na]
    }
    
    a <- prod(sapply(list_prod(ex2_list, esig_list), sum)[-k])

    left_mult <- list_prod(ex_list, esig_list)
    left_mult[[k]] <- diag(p[k])
    b <- as.numeric(tensr::atrans(Y, lapply(left_mult, t)))

    sebetahat <- 1 / sqrt(esig_list[[k]] * a)
    betahat <- b / a

    ATM = ashr::ash(betahat = betahat, sebetahat = sebetahat,
                    method = "fdr", mixcompdist = mixcompdist, nullweight = nullweight)

    post_mean <- ATM$PosteriorMean
    post_sd <- ATM$PosteriorSD

    ex_list[[k]] <- post_mean
    ex2_list[[k]] <- post_mean ^ 2 + post_sd ^ 2

    prob_zero <- ATM$ZeroProb
    pi0 <- ATM$fitted.g$pi[1]

    return(list(ex_list = ex_list, ex2_list = ex2_list, prob_zero = prob_zero,
                pi0 = pi0))
}

#' Update the kth mode diagional Kronecker structured covariance
#' matrix.
#'
#'
#' @inheritParams tupdate_kron_modek
#' @param beta A non-negative vector of numerics. Prior rate
#'        parameters for the mode-k variances.
#'
#' @return \code{post_rate} The variational posterior rate parameter.
#'
#'
tupdate_kron_sig <- function(Y, ex_list, ex2_list, esig_list, k, beta = 0, which_na = NULL) {
    p <- dim(Y)
    n <- length(p)

    current_mean <- form_outer(ex_list)
    current_var <- 1 / form_outer(esig_list)

    A <- Y
    if (!is.null(which_na)) {
        A[which_na] <- sqrt(current_mean[which_na] ^ 2 + current_var[which_na])
    }
    
    for (a_index in (1:n)[-k]) {
        AM <- sqrt(esig_list[[a_index]]) * tensr::mat(A, a_index)
        AMA <- array(AM, dim = c(p[a_index], p[-a_index]))
        A <- aperm(AMA, match(1:n, c(a_index, (1:n)[-a_index])))
    }
    a <- rowSums(tensr::mat(A, k) ^ 2)

    ## Alternative way to do calculation.
    ## tempa <- lapply(lapply(esig_list, sqrt), diag)
    ## tempa[[k]] <- diag(p[k])
    ## rowSums((tensr::mat(tensr::atrans(Y, tempa), k)) ^ 2)

    if (!is.null(which_na)) {
        Y[which_na] <- current_mean[which_na]
    }
    
    b_mult <- list_prod(ex_list, esig_list)
    b_mult[[k]] <- diag(p[k])
    b <- as.numeric(tensr::atrans(Y, lapply(b_mult, t)))

    c_mult <- list_prod(ex2_list, esig_list)
    c <- prod(sapply(c_mult[-k], sum))

    post_rate <- (a - 2 * b * ex_list[[k]] + c * ex2_list[[k]]) / 2 + beta

    return(post_rate)
}


#' Obtain initial estimates of each mode's components.
#'
#' The inital values of each component are taken to be proportional to
#' the first singular vector of each mode-k matricization of the data
#' array.
#'
#'
#'
#' @inheritParams tflash
#' @param which_na Either NULL or an array the same dimension as
#'     \code{Y} indicating if an element of \code{Y} is missing
#'     (\code{TRUE}) or observed (\code{FALSE}).
#'
#' @return \code{ex_list} A list of vectors of the starting expected
#'     values of each component.
#'
#'
#'     \code{ex2_list} A vector of starting values of \eqn{E[x^2]}
#'
#'
#' @author David Gerard
#'
tinit_kron_components <- function(Y, which_na = NULL, start = c("first_sv", "random"),
                                  known_factors = NULL, known_modes = NULL,
                                  homo_modes = NULL) {
    p <- dim(Y)
    n <- length(p)

    start <- match.arg(start, c("first_sv", "random"))
    
    if (!is.null(which_na)) {
        Y[which_na] <- mean(Y, na.rm = TRUE)
    }
    
    x <- vector(mode = "list", length = n)

    if (is.null(known_modes)) {
        unknown_modes <- 1:n
    } else {
        unknown_modes <- (1:n)[-known_modes]
        known_f_index <- 1
        for(k in known_modes) {
            x[[k]] <- known_factors[[known_f_index]]
            known_f_index <- known_f_index + 1
        }
    }
    

    
    if (start == "first_sv") {
        for(k in unknown_modes) {
            sv_out <- tryCatch(irlba::irlba(tensr::mat(Y, k), nv = 0, nu = 1),
                               error = function(){"do_full"})
            if (identical(sv_out, "do_full")) {
                sv_out <- svd(tensr::mat(Y, k))$u[,1]
            }
            x[[k]] <-  c(sv_out$u) * sign(c(sv_out$u)[1]) ## for identifiability reasons
        }
    } else if (start == "random") {
        for (k in unknown_modes) {
            x[[k]] <- rnorm(p[k])
            x[[k]] <- x[[k]] / sqrt(sum(x[[k]] ^ 2))
        }
    }
    d1 <- as.numeric(tensr::atrans(Y, lapply(x, t)))
    if (d1 < 0) {
        x[[min(unknown_modes)]] <- x[[min(unknown_modes)]] * -1
        d1 <- abs(d1)
    }

    ## scale the unknown factors
    ex_list <- x
    xmult <- d1 ^ (1 / length(unknown_modes))
    for (mode_index in unknown_modes) {
        ex_list[[mode_index]] <- x[[mode_index]] * xmult
    }

    ex2_list <- lapply(ex_list, FUN = function(x) { x ^ 2 })

    ## huberized initial sigma est
    R <- Y - form_mean(ex_list)
    esig_list <- diag_mle(R, homo_modes = homo_modes)

    ## If want to run a few iterations of updating sig
    ## itermax <- 10
    ## for(iter_index in 1:itermax){
    ##     gamma <- prod(p) / (2 * p)
    ##     for(mode_index in 1:n) {
    ##         delta <- tupdate_kron_sig(Y = Y, ex_list = ex_list, ex2_list = ex2_list,
    ##                                   esig_list = esig_list,
    ##                                   k = mode_index, which_na = which_na)
    ##         esig_list[[mode_index]] <- gamma[mode_index] / delta
    ##         cat(esig_list[[mode_index]][1], "\n")
    ##     }
    ## }

    return(list(ex_list = ex_list, ex2_list = ex2_list, esig_list = esig_list))
}

#' MLE of mean zero Kronecker variance model.
#' 
#' @param R An array of numerics.
#' @param itermax A positive integer. The maximium number of
#'     iterations to perform.
#' @param tol A positive numeric. The stopping criterion.
#' @param homo_modes  A vector of integers. If \code{var_type =
#'     "kronecker"} then \code{homo_modes} indicates which modes are
#'     assumed to be homoscedastic.
#' 
#' 
#' 
#' @return \code{esig_list} A list of vectors of numerics. The MLEs of
#'     the precisions, not the variances.
#' 
#' @author David Gerard
#' 
#' @export
diag_mle <- function(R, itermax = 100, tol = 10^-3, homo_modes = NULL) {
    p <- dim(R)
    n <- length(p)

    if (is.null(homo_modes)) {
        hetero_modes <- 1:n
    } else {
        hetero_modes <- (1:n)[-homo_modes]
    }
    
    ## huberized initial sigma est
    resid2 <- R ^ 2
    esig_list <- list()
    for(mode_index in hetero_modes) {
        z <- apply(resid2, mode_index, mean)
        quants <- quantile(z, c(0.25, 0.75))
        z[z > quants[2]] <- quants[2]
        z[z < quants[1]] <- quants[1]
        esig_list[[mode_index]] <- z
    }

    if (!is.null(homo_modes)) {
        for(mode_index in homo_modes) {
            esig_list[[mode_index]] <- rep(1, length = p[mode_index])
        }
    }
    naive_est <- rescale_factors(esig_list)

    ## Now run MLE algorithm
    err <- tol + 1

    iter_index <- 1
    while(err > tol & iter_index < itermax) {
        esig_list_old <- esig_list
        
        for(mode_index in hetero_modes) {
            esig_list <-
        temp <-         mle_update_modek(R = R, esig_list = esig_list, k = mode_index)
        }
        esig_list <- rescale_factors(esig_list)
        iter_index <- iter_index + 1
        err <- sum(abs(esig_list[[1]] - esig_list_old[[1]]))
    }

    esig_list <- lapply(esig_list, function(x) { 1/x }) # to make it precisions
    return(esig_list)
}


#' In MLE algorithm for mean zero diagonal Kronecker structured
#' variance model, update mode k.
#' 
#' @param R An array of numerics.
#' @param esig_list A list of current values of the variances.
#' @param k The current mode to update
#' 
#' @author David Gerard
#' 
#' 
#' 
#' 
mle_update_modek <- function(R, esig_list, k) {
    p <- dim(R)
    n <- length(p)
    
    A <- R
    for (a_index in (1:n)[-k]) {
        AM <- (1 / sqrt(esig_list[[a_index]])) * tensr::mat(A, a_index)
        AMA <- array(AM, dim = c(p[a_index], p[-a_index]))
        A <- aperm(AMA, match(1:n, c(a_index, (1:n)[-a_index])))
    }
    a <- rowSums(tensr::mat(A, k) ^ 2)

    esig_list[[k]] <- a / prod(p[-k])
    return(esig_list)
}


#' Element-wise multiplication of two lists.
#'
#' @param A A list.
#' @param B A list.
#'
#' @author David Gerard
#'
list_prod <- function(A, B) {
    C <- list()
    for(index in 1:length(A)) {
        C[[index]] <- A[[index]] * B[[index]]
    }
    return(C)
}




## p <- c(10,10,10)
## u <- list()
## u[[1]] <- rnorm(p[1])
## u[[2]] <- rnorm(p[2])
## u[[3]] <- rnorm(p[3])

## Theta <- outer(outer(u[[1]], u[[2]], "*"), u[[3]], "*")
## E <- array(rnorm(prod(p)), dim = p)
## Y <- Theta + E

## yup_out <- tflash_kron(Y = Y, alpha = 100, beta = 100)
## yup_out

## plot(yup_out$post_mean[[1]], u[[1]])
## plot(yup_out$post_mean[[2]], u[[2]])
## plot(yup_out$post_mean[[3]], u[[3]])



## library(ggplot2)
## set.seed(31)
## n <- 100
## p <- 100
## u <- rnorm(n, mean = 10)
## v <- rnorm(p, mean = 10)
## row_cov_half <- diag(sqrt(seq(1, 5, length = n)))
## col_cov_half <- diag(sqrt(seq(1, 5, length = p)))
## E <- row_cov_half %*% matrix(rnorm(n * p), nrow = n) %*% col_cov_half
## Y <- u %*% t(v) + E

## pmiss <- 0.6
## Omega <- sort(sample(1:(n * p), size = round(pmiss * n * p)))
## Y[Omega] <- NA

## tout <- tflash_kron(Y = Y)
## qplot(u, tout$post_mean[[1]], xlab = "Truth", ylab = "Estimate", main = "Loadings")
## qplot(v, tout$post_mean[[2]], xlab = "Truth", ylab = "Estimate", main = "Factors")
## qplot(diag(row_cov_half)^2, 1 / tout$sigma_est[[1]], xlab = "Truth",
##       ylab = "Estimate", main = "Row Cov")
## qplot(diag(col_cov_half)^2, 1 / tout$sigma_est[[2]], xlab = "Truth",
##       ylab = "Estimate", main = "Col Cov")

## Y[Omega] <- 0
## tout <- tflash_kron(Y = Y)
## qplot(u, tout$post_mean[[1]], xlab = "Truth", ylab = "Estimate", main = "Loadings")
## qplot(v, tout$post_mean[[2]], xlab = "Truth", ylab = "Estimate", main = "Factors")
## qplot(diag(row_cov_half)^2, 1 / tout$sigma_est[[1]], xlab = "Truth",
##       ylab = "Estimate", main = "Row Cov")
## qplot(diag(col_cov_half)^2, 1 / tout$sigma_est[[2]], xlab = "Truth",
##       ylab = "Estimate", main = "Col Cov")
