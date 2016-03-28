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
#'
#' @author David Gerard
#'
#' @export
tflash <- function(Y, tol = 10^-5, itermax = 100, alpha = 0, beta = 0, mixcompdist = "normal") {
    p <- dim(Y)
    n <- length(p)
    ssY <- sum(Y ^ 2)

    ## Initialize Parameters --------------------------------------------
    init_return <- tinit_components(Y)
    ex_list <- init_return$ex_list # list of expected value of components.
    ex2_vec <- init_return$ex2_vec # vector of expected value of x'x

    gamma <- prod(p) / 2 + alpha # posterior shape parameter. Does not change.
    delta <- tupdate_sig(ssY = ssY, Y = Y, ex_list = ex_list, # posterior rate parameter.
                         ex2_vec = ex2_vec, beta = beta)
    esig <- gamma / delta # expected value of precision.


    iter_index <- 1
    err <- tol + 1
    
    while(iter_index < itermax & err > tol) {
        old_sig <- esig
        for(mode_index in 1:n) {
            if(sum(abs(ex_list[[mode_index]])) < 10^-6) {
                ex_list <- lapply(ex_list, FUN = function(x) { rep(0, length = length(x)) })
                break
            }
            
            tupdate_out <- tupdate_modek(Y = Y, ex_list = ex_list, ex2_vec = ex2_vec,
                                         esig = esig, k = mode_index, mixcompdist = mixcompdist)
            ex_list <- tupdate_out$ex_list
            ex2_vec <- tupdate_out$ex2_vec
            
            if(sum(abs(ex_list[[mode_index]])) < 10^-6) {
                ex_list <- lapply(ex_list, FUN = function(x) { rep(0, length = length(x)) })
                break
            }
            
            delta <- tupdate_sig(ssY = ssY, Y = Y, ex_list = ex_list,
                                 ex2_vec = ex2_vec, beta = beta)
            esig <- gamma / delta
        }
        iter_index <- iter_index + 1
        err <- abs(old_sig/esig - 1)
    }
    return(list(postMean = ex_list, sigma_est = esig))
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
tupdate_modek <- function(Y, ex_list, ex2_vec, esig, k, mixcompdist = "normal") {
    p <- dim(Y)
    n <- length(p)
    
    left_list <- ex_list
    left_list[[k]] <- diag(p[k])
    b <- esig * as.numeric(tensr::atrans(Y, lapply(left_list, t)))
    a <- esig * prod(ex2_vec[-k])

    betahat <- b / a
    sebetahat <- 1 / sqrt(a)
    
    ATM = ashr::ash(betahat = betahat, sebetahat = sebetahat,
                    method = "fdr", mixcompdist = mixcompdist)

    post_mean <- ATM$PosteriorMean
    post_sd <- ATM$PosteriorSD
    ## mix_prop <- ATM$fitted.g$pi
    ex_list[[k]] <- post_mean
    ex2_vec[k] <- sum(post_sd ^ 2 + post_mean ^ 2)

    return(list(ex_list = ex_list, ex2_vec = ex2_vec))
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
#'     \code{ex2_vec} A vector of starting values of \eqn{E[x'x]}
#'
#'
#' @author David Gerard
#'
#' @export
tinit_components <- function(Y) {
    p <- dim(Y)
    n <- length(p)
    x <- vector(mode = "list", length = n)
    for(k in 1:n) {
        sv_out <- irlba::irlba(tensr::mat(Y, k), nv = 0, nu = 1)
        x[[k]] <-  c(sv_out$u) * sign(c(sv_out$u)[1]) ## for identifiability reasons
    }
    d1 <- abs(as.numeric(tensr::atrans(Y, lapply(x, t))))

    ex_list <- lapply(x, FUN = function(x, xmult) { x * xmult }, xmult = d1 ^ (1 / n))

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
#' @param ssY A positive numeric. The sum of squares of the data array
#'     \code{Y}.
#' @param Y An array of numerics. The data array.
#' @param ex_list A list of vectors of numerics. The expected values
#'     of the components of each mode.
#' @param ex2_vec A vector of positive numerics. The expected value of
#'     \eqn{x'x} for each mode.
#' @param beta The prior rate parameter.
#'
#' @return \code{delta_new} The updated rate parameter of the
#'     variational distribution of the precision.
#'
#' @export
#'
#' @author David Gerard
#'
tupdate_sig <- function(ssY, Y, ex_list, ex2_vec, beta) {
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



## n <- 10
## p <- 10
## Theta <- rnorm(n) %*% t(rnorm(p))
## E <- matrix(rnorm(n * p), nrow = n)
## Y <- Theta + E
