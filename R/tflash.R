#' Tensor Factor Loading Adaptive SHrinkage (T-FLASH).
#'
#' Fits a rank-1 tensor mean model with a homoscedastic error term. It
#' does this via a variational Bayesian approach with a unimodal prior
#' on the components of the mean tensor.
#'
#'
#'
#' @param Y An array of numerics. The data.
#' @param tol A positive numeric. The stopping criterion for the VEM.
#' @param itermax A positive integer. The maximum number of iterations
#'     to run the VEM
#' @param alpha A non-negative numeric. The prior shape parameter for
#'     the variance. Defaults to zero.
#' @param beta A non-negative numeric. The prior rate parameter for
#'     the variance. Defaults to zero.
#' @param mixcompdist The mixing distribution to assume. Defaults to
#'     normal. Options are those available in the \code{ashr} package.
#' @param var_type A string. What variance model should we assume?
#'     Options are homoscedastic noise (\code{"homoscedastic"}) or
#'     Kronecker structured variance (\code{kronecker}).
#' @param sig_start_itermax A positive integer. The number of
#'     iterations to run in initializing the precision before starting
#'     the VEM. Defaults to 10.
#' @param nullweight A numeric greater than or equal to 1. The penalty
#'     term on the probability of zero.
#' @param print_update A logical. Should we print notifications on how
#'     far along the optimization is?
#' @param start How should we choose the starting values? Either using
#'     the first singular vector along each mode (\code{"first_sv"})
#'     or randomly (\code{"random"}).
#' @param known_factors A list of known factors for the modes
#'     indicated in \code{known_modes}. Defaults to \code{NULL}, where
#'     all factors are assumed to be unknown.
#' @param known_modes A vector of integers. The modes that are
#'     known. Should be the same length as \code{known_factors}.
#' @param homo_modes A vector of integers. If \code{var_type =
#'     "kronecker"} then \code{homo_modes} indicates which modes are
#'     assumed to be homoscedastic.
#'
#'
#' @author David Gerard
#'
#' @export
tflash <- function(Y, var_type = c("homoscedastic", "kronecker"), tol = 10^-5,
                   itermax = 100, alpha = 0, beta = 0, mixcompdist = "normal",
                   sig_start_itermax = 10, nullweight = 10, print_update = FALSE,
                   start = c("first_sv", "random"), known_factors = NULL,
                   known_modes = NULL, homo_modes = NULL) {
    p <- dim(Y)
    n <- length(p)

    var_type <- match.arg(var_type, c("homoscedastic", "kronecker"))

    start = match.arg(start, c("first_sv", "random"))

    ## checks
    if (!is.null(known_modes)) {
        dim_factors <- sapply(known_factors, length)
        if (is.null(known_factors)) {
            stop("known_modes is not NULL but known_factors is NULL")
        } else if (length(known_modes) != length(known_factors)) {
            stop ("known_modes and known_factors must be of same length")
        } else if (all(p[known_modes] != dim_factors)) {
            stop("known_factors not the same dimension as modes of Y")
        }
    }

    if (var_type == "homoscedastic") {
        flash_out <- tflash_homo(Y = Y, tol = tol, itermax = itermax, alpha = alpha,
                                 beta = beta, mixcompdist = mixcompdist,
                                 sig_start_itermax = sig_start_itermax,
                                 nullweight = nullweight, print_update = print_update,
                                 start = start, known_factors = known_factors,
                                 known_modes = known_modes)
    } else if (var_type == "kronecker") {
        flash_out <- tflash_kron(Y = Y, tol = tol, itermax = itermax, alpha = alpha,
                                 beta = beta, mixcompdist = mixcompdist,
                                 nullweight = nullweight, print_update = print_update,
                                 start = start, known_factors = known_factors,
                                 known_modes = known_modes, homo_modes = homo_modes)
    }

    return(flash_out)
}


#' Homoscedastic variance implementation of T-FLASH.
#'
#' @inheritParams tflash
#'
#' @export
#'
#' @author David Gerard
#'
tflash_homo <- function(Y, tol = 10^-5, itermax = 100, alpha = 0, beta = 0,
                        mixcompdist = "normal", sig_start_itermax = 10, nullweight = 10,
                        print_update = FALSE, start = c("first_sv", "random"),
                        known_factors = NULL, known_modes = NULL) {
    p <- dim(Y)
    n <- length(p)


    if (is.null(known_modes)) {
        unknown_modes <- 1:n
    } else {
        unknown_modes <- (1:n)[-known_modes]
    }


    ssY_obs <- sum(Y ^ 2, na.rm = TRUE)

    start <- match.arg(start, c("first_sv", "random"))

    which_na <- is.na(Y)
    if(all(!which_na)) {
        which_na <- NULL
    }

    num_na <- sum(which_na)

    ## Initialize Parameters --------------------------------------------
    init_return <- tinit_components(Y = Y, which_na = which_na, start = start,
                                    known_factors = known_factors,
                                    known_modes = known_modes)
    ex_list <- init_return$ex_list # list of expected value of components.
    ex2_vec <- init_return$ex2_vec # vector of expected value of x'x




    ## posterior shape parameter. Does not change.
    gamma <- prod(p) / 2 + alpha

    ## run through a few iterations of updating ssY and esig to get
    ## initial values of esig
    iter_index <- 1
    sig_err <- tol + 1
    esig <- 1 / var(c(Y), na.rm = TRUE)
    while (iter_index < sig_start_itermax & sig_err > tol) {
        esig_old <- esig
        delta <- tupdate_sig(ssY_obs = ssY_obs, Y = Y, ex_list = ex_list, esig = esig,
                             ex2_vec = ex2_vec, beta = beta, which_na = which_na)
        esig <- gamma / delta

        sig_err <- abs(esig / esig_old - 1)
        iter_index <- iter_index + 1
        ##cat(esig, "\n")
    }

    init_return$esig <- esig

    iter_index <- 1
    err <- tol + 1

    prob_zero <- vector(length = n, mode = "list")
    pi0vec <- rep(NA, length = n)

    while(iter_index < itermax & err > tol) {
        old_sig <- esig
        for(mode_index in unknown_modes) {
            if(sum(abs(ex_list[[mode_index]])) < 10^-6) {
                ex_list <- lapply(ex_list, FUN = function(x) { rep(0, length = length(x)) })
                break
            }

            tupdate_out <- tupdate_modek(Y = Y, ex_list = ex_list, ex2_vec = ex2_vec,
                                         esig = esig, k = mode_index, mixcompdist = mixcompdist,
                                         which_na = which_na, nullweight = nullweight)

            prob_zero[[mode_index]] <- tupdate_out$prob_zero
            pi0vec[mode_index] <- tupdate_out$pi0

            if(sum(abs(tupdate_out$ex_list[[mode_index]])) < 10^-6) {
              ex_list <- lapply(ex_list, FUN = function(x) { rep(0, length = length(x)) })
              break
            }

            if(any(is.na(tupdate_out$ex2_vec))) { stop("na in ex2_vec") }


            ex_list <- tupdate_out$ex_list
            ex2_vec <- tupdate_out$ex2_vec



            delta <- tupdate_sig(ssY_obs = ssY_obs, Y = Y, ex_list = ex_list, esig = esig,
                                 ex2_vec = ex2_vec, beta = beta, which_na = which_na)
            esig <- gamma / delta
            ##cat(esig, "\n")
        }
        iter_index <- iter_index + 1
        err <- abs(old_sig/esig - 1)

        if (print_update & iter_index %% 5 == 0) {
            cat("Iteration:", iter_index, "\n")
            cat("Stop Crit:", err, "\n\n")
        }
    }
    return(list(post_mean = ex_list, sigma_est = esig, prob_zero = prob_zero,
                pi0vec = pi0vec,
                num_iter = iter_index, init_return = init_return))
}





#' Update the mode k variational density.
#'
#'
#' @inheritParams tupdate_sig
#' @param esig A positive numeric. The expected value of the
#'     precision.
#' @param k A positive integer. The current mode to update.
#' @param mixcompdist The mixing distribution to assume. Defaults to
#'     normal. Options are those available in the \code{ashr} package.
#' @param which_na Either NULL (when complete data) or an array of
#'     logicals the same dimension as \code{Y}, indicating if the
#'     observation is missing (\code{TRUE}) or observed
#'     (\code{FALSE}).
#' @param nullweight A numeric greater than or equal to 1. The penalty
#'     term on the probability of a component being zero.
#'
#' @return \code{ex_list} A list of vectors of the starting expected
#'     values of each component.
#'
#'     \code{ex2_vec} A vector of starting values of \eqn{E[x'x]}
#'
#'
#' @export
#'
#' @author David Gerard
#'
tupdate_modek <- function(Y, ex_list, ex2_vec, esig, k, mixcompdist = "normal",
                          which_na = NULL, nullweight = 10) {
    p <- dim(Y)
    n <- length(p)

    if (!is.null(which_na)) {
        Y[which_na] <- form_outer(ex_list)[which_na]
    }

    left_list <- ex_list
    left_list[[k]] <- diag(p[k])
    b <- esig * as.numeric(tensr::atrans(Y, lapply(left_list, t)))
    a <- esig * prod(ex2_vec[-k])

    betahat <- b / a
    sebetahat <- 1 / sqrt(a)

    ATM <- ashr::ash(betahat = betahat, sebetahat = sebetahat,
                    method = "fdr", mixcompdist = mixcompdist, nullweight = nullweight,
                    outputlevel = 4, optmethod = "mixIP")

    if(any(is.na(ATM$PosteriorMean))) {
      ATM = ashr::ash(betahat = betahat, sebetahat = sebetahat,
                      method = "fdr", mixcompdist = mixcompdist, nullweight = nullweight,
                      outputlevel = 4, optmethod = "mixEM")
    }

    post_mean <- ATM$PosteriorMean
    post_sd <- ATM$PosteriorSD
    ## mix_prop <- ATM$fitted.g$pi
    ex_list[[k]] <- post_mean
    ex2_vec[k] <- sum(post_sd ^ 2 + post_mean ^ 2)
    prob_zero <- ATM$ZeroProb
    pi0 <- ATM$fitted.g$pi[1]
    return(list(ex_list = ex_list, ex2_vec = ex2_vec, prob_zero = prob_zero, pi0 = pi0))
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
#' @param which_na Either NULL (when complete data) or an array of
#'     logicals the same dimension as \code{Y}, indicating if the
#'     observation is missing (\code{TRUE}) or observed
#'     (\code{FALSE}).
#'
#' @return \code{ex_list} A list of vectors of the starting expected
#'     values of each component.
#'
#'
#'     \code{ex2_vec} A vector of starting values of \eqn{E[x'x]}
#'
#'
#' @author David Gerard
#'
#' @export
tinit_components <- function(Y, which_na = NULL, start = c("first_sv", "random"),
                             known_factors = NULL, known_modes = NULL) {
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
            x[[k]] <- tryCatch(
            {
                sv_out <- irlba::irlba(tensr::mat(Y, k), nv = 0, nu = 1)
                c(sv_out$u) * sign(c(sv_out$u)[1]) ## for identifiability reasons
            },
            error =
                {
                    svd(tensr::mat(Y, k))$u[, 1]
                }
            )
            
        }
    } else if (start == "random") {
        for (k in unknown_modes) {
            x[[k]] <- rnorm(p[k])
            x[[k]] <- x[[k]] / sqrt(sum(x[[k]] ^ 2))
        }
    }

    
    xscaled <- x
    if(!is.null(known_modes)) {
        fnorm_xknown <- rep(NA, length = length(known_modes))
        km_index <- 1
        for (mode_index in known_modes) {
            fnorm_xknown[km_index] <- sqrt(sum(x[[mode_index]] ^ 2))
            xscaled[[mode_index]] <- x[[mode_index]] / fnorm_xknown[km_index]
            km_index <- km_index + 1
        }
    }
    
    d1 <- as.numeric(tensr::atrans(Y, lapply(xscaled, t)))
    if (d1 < 0) {
        x[[min(unknown_modes)]] <- x[[min(unknown_modes)]] * -1
        d1 <- abs(d1)
    }

    if(!is.null(known_modes)) {
        d1 <- d1 / prod(fnorm_xknown)
    }

    ## scale the unknown factors
    ex_list <- x
    xmult <- d1 ^ (1 / length(unknown_modes))
    for (mode_index in unknown_modes) {
        ex_list[[mode_index]] <- x[[mode_index]] * xmult
    }

    ex2_vec <- sapply(ex_list, FUN = function(x) { sum(x ^ 2) })
    return(list(ex_list = ex_list, ex2_vec = ex2_vec))
}



#' Obtain initial estimates of each mode's components.
#'
#' The inital values of each component are taken to be proportional to
#' the first singular vector of each mode-k matricization of the data
#' array. The constant is placed entirely on the first mode.
#'
#'
#'
#' @param Y An array of numerics.
#'
#' @return \code{ex_list} A list of vectors of the starting expected
#'     values of each component.
#'
#'
#'     \code{ex2_vec} A vector of starting values of \eqn{E[x'x]}
#'
#'
#' @author David Gerard
#'
#' @export
tinit_components_one <- function(Y) {
    p <- dim(Y)
    n <- length(p)
    x <- vector(mode = "list", length = n)
    for(k in 1:n) {
        sv_out <- irlba::irlba(tensr::mat(Y, k), nv = 0, nu = 1)
        x[[k]] <-  c(sv_out$u) * sign(c(sv_out$u)[1]) ## for identifiability reasons
    }
    d1 <- as.numeric(tensr::atrans(Y, lapply(x, t)))

    ex_list <- x
    ex_list[[1]] <- ex_list[[1]] * d1

    ex2_vec <- sapply(ex_list, FUN = function(x) { sum(x ^ 2) })
    return(list(ex_list = ex_list, ex2_vec = ex2_vec))
}

#' Update the gamma distribution parameters of the variational
#' distribution of the precision sig.
#'
#' @param ssY_obs A positive numeric. The sum of squares of the data array
#'     \code{Y}.
#' @param Y An array of numerics. The data array.
#' @param ex_list A list of vectors of numerics. The expected values
#'     of the components of each mode.
#' @param ex2_vec A vector of positive numerics. The expected value of
#'     \eqn{x'x} for each mode.
#' @param beta The prior rate parameter.
#' @param which_na Either NULL (when complete data) or an array of
#'     logicals the same dimension as \code{Y}, indicating if the
#'     observation is missing (\code{TRUE}) or observed
#'     (\code{FALSE}).
#' @param esig A positive numeric. The variational expectation of the
#'     variance.
#'
#'
#' @return \code{delta_new} The updated rate parameter of the
#'     variational distribution of the precision.
#'
#' @export
#'
#' @author David Gerard
#'
tupdate_sig <- function(ssY_obs, Y, ex_list, ex2_vec, beta, esig, which_na = NULL) {

    if (!is.null(which_na)) {
        current_mean <- form_outer(ex_list)
        Y[which_na] <- current_mean[which_na]
        ssY <- ssY_obs + sum(current_mean[which_na] ^ 2 + 1 / esig)
    } else {
        ssY <- ssY_obs
    }

    mid_term <- as.numeric(tensr::atrans(Y, lapply(ex_list, t)))
    delta_new <- (ssY - 2 * mid_term + prod(ex2_vec)) / 2 + beta
    return(delta_new)
}

#' Forms outer product from a list of vectors.
#'
#'
#'
#'
#' @param x A list of vectors of positive numerics.
#'
#' @return An array of numerics that is the outer products of all of
#'     the vectors in \code{x}.
#'
#'
#' @export
#'
#' @author David Gerard
form_outer <- function(x) {
    n <- length(x)
    theta <- x[[1]]
    if (n == 1) {
        return(theta)
    }
    for(mode_index in 2:n) {
        theta <- outer(theta, x[[mode_index]], "*")
    }
    return(theta)
}



## n <- 50
## p <- 50
## Theta <- rnorm(n, mean = 1) %*% t(rnorm(p, mean = 1))
## E <- matrix(rnorm(n * p), nrow = n)
## Y <- Theta + E

## flash(Y, nonnegative = TRUE)

## pi0 <- 0.8
## Omega <- sample(1:(n * p), size = n * p * pi0)
## Y[Omega] <- NA

## tout <- tflash(Y)
## tout

## fout <- flash(Y)

## plot(svd(Theta)$u[, 1], tout$post_mean[[1]])
## plot(svd(Theta)$v[, 1], tout$post_mean[[2]])
## plot(svd(Theta)$u[, 1], fout$l)
## plot(svd(Theta)$v[, 1], fout$f)

## cor(svd(Theta)$u[, 1], tout$post_mean[[1]])
## cor(svd(Theta)$v[, 1], tout$post_mean[[2]])
## cor(svd(Theta)$u[, 1], fout$l)
## cor(svd(Theta)$v[, 1], fout$f)
