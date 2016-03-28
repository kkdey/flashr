
#' Tensor FLASH assuming a diagonal Kronecker structured covariance model.
#'
#' @inheritParams tflash
#'
#' @export
#'
#' @author David Gerard
tflash_kron <- function(Y, tol = 10^-5, itermax = 100, alpha = 0, beta = 0,
                        mixcompdist = "normal") {
    p <- dim(Y)
    n <- length(p)

    init_return <- tinit_kron_components(Y)
    ex_list <- init_return$ex_list # list of expected value of components.
    ex2_list <- init_return$ex2_list # list of expected value of x^2
    esig_list <- init_return$esig_list

    post_rate <- list()
    post_shape <- list()
    for(mode_index in 1:n) {
        post_shape[[mode_index]] <- rep(prod(p[-mode_index])) / 2 + alpha
        post_rate[[mode_index]] <- post_shape[[mode_index]] / esig_list[[mode_index]]
    }

    iter_index <- 1
    not_all_zero <- TRUE
    while (iter_index <= itermax) {
        for (mode_index in 1:n) {
            if(sum(abs(ex_list[[mode_index]])) < 10^-6) {
                ex_list <- lapply(ex_list, FUN = function(x) { rep(0, length = length(x)) })
                not_all_zero <- FALSE
                break
            }

            t_out <- tupdate_kron_modek(Y = Y, ex_list = ex_list, ex2_list = ex2_list,
                                        esig_list = esig_list, k = mode_index,
                                        mixcompdist = mixcompdist)
            ex_list <- t_out$ex_list
            ex2_list <- t_out$ex2_list

            if(sum(abs(ex_list[[mode_index]])) < 10^-6) {
                ex_list <- lapply(ex_list, FUN = function(x) { rep(0, length = length(x)) })
                not_all_zero <- FALSE
                break
            }


            post_rate[[mode_index]] <- tupdate_kron_sig(Y = Y, ex_list = ex_list,
                                                        ex2_list = ex2_list,
                                                        esig_list = esig_list,
                                                        k = mode_index, beta = beta)
            esig_list[[mode_index]] <- post_shape[[mode_index]] / post_rate[[mode_index]]
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
    }
    return(list(post_mean = ex_list, sigma_est = esig_list, num_iter = iter_index))
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
tupdate_kron_modek <- function(Y, ex_list, ex2_list, esig_list, k, mixcompdist = "normal") {
    p <- dim(Y)
    a <- prod(sapply(list_prod(ex2_list, esig_list), sum)[-k])

    left_mult <- list_prod(ex_list, esig_list)
    left_mult[[k]] <- diag(p[k])
    b <- as.numeric(tensr::atrans(Y, lapply(left_mult, t)))

    sebetahat <- 1 / sqrt(esig_list[[k]] * a)
    betahat <- b / a

    ATM = ashr::ash(betahat = betahat, sebetahat = sebetahat,
                    method = "fdr", mixcompdist = mixcompdist, nullweight = 1)

    post_mean <- ATM$PosteriorMean
    post_sd <- ATM$PosteriorSD

    ex_list[[k]] <- post_mean
    ex2_list[[k]] <- post_mean ^ 2 + post_sd ^ 2

    return(list(ex_list = ex_list, ex2_list = ex2_list))
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
tupdate_kron_sig <- function(Y, ex_list, ex2_list, esig_list, k, beta = 0) {
    p <- dim(Y)
    n <- length(p)

    A <- Y
    for (a_index in (1:n)[-k]) {
        AM <- sqrt(esig_list[[a_index]]) * tensr::mat(A, a_index)
        AMA <- array(AM, dim = c(p[k], p[-k]))
        A <- aperm(AMA, match(1:n, c(a_index, (1:n)[-a_index])))
    }
    a <- rowSums(tensr::mat(A, k) ^ 2)

    ## Alternative way to do calculation.
    ## tempa <- lapply(lapply(esig_list, sqrt), diag)
    ## tempa[[k]] <- diag(p[k])
    ## rowSums((tensr::mat(tensr::atrans(Y, tempa), k)) ^ 2)
    
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
#' @param Y An array of numerics.
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
tinit_kron_components <- function(Y) {
    p <- dim(Y)
    n <- length(p)
    x <- vector(mode = "list", length = n)
    for(k in 1:n) {
        sv_out <- irlba::irlba(tensr::mat(Y, k), nv = 0, nu = 1)
        x[[k]] <-  c(sv_out$u) * sign(c(sv_out$u)[1]) ## for identifiability reasons
    }
    d1 <- abs(as.numeric(tensr::atrans(Y, lapply(x, t))))

    ex_list <- lapply(x, FUN = function(x, xmult) { x * xmult }, xmult = d1 ^ (1 / n))

    ex2_list <- lapply(ex_list, FUN = function(x) { x ^ 2 })

    ## huberized initial sigma est
    R <- Y - form_mean(ex_list)
    esig_list <- diag_mle(R)

    return(list(ex_list = ex_list, ex2_list = ex2_list, esig_list = esig_list))
}

#' MLE of mean zero Kronecker variance model.
#' 
#' @param R An array of numerics.
#' 
#' @return \code{esig_list} A list of vectors of numerics. The MLEs of
#'     the precisions, not the variances.
#' 
#' @author David Gerard
#' 
#' @export
diag_mle <- function(R, itermax = 100, tol = 10^-3) {
    p <- dim(R)
    n <- length(p)
    
    ## huberized initial sigma est
    resid2 <- R ^ 2
    esig_list <- list()
    for(mode_index in 1:n) {
        x <- apply(resid2, mode_index, mean)
        quants <- quantile(x, c(0.25, 0.75))
        x[x > quants[2]] <- quants[2]
        x[x < quants[1]] <- quants[1]
        esig_list[[mode_index]] <- x
    }
    naive_est <- rescale_factors(esig_list)

    ## Now run MLE algorithm
    err <- tol + 1
    iter_index <- 1
    while(err > tol & iter_index < itermax) {
        esig_list_old <- esig_list
        for(mode_index in 1:n) {
            esig_list <- mle_update_modek(R = R, esig_list = esig_list, k = mode_index)
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
        AMA <- array(AM, dim = c(p[k], p[-k]))
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



## set.seed(31)
## n <- 50
## p <- 50
## u <- rnorm(n)
## v <- rnorm(p)
## row_cov_half <- diag(sqrt(seq(1, 5, length = n)))
## col_cov_half <- diag(sqrt(seq(1, 5, length = p)))
## E <- row_cov_half %*% matrix(rnorm(n * p), nrow = n) %*% col_cov_half
## Y <- u %*% t(v) + E

