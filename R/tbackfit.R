#' Iteratively do T-FLASH on the residuals.
#'
#' This function will estimate a rank-1 tensor using independent
#' unimodal priors on the components via a variational EM
#' algorithm. It will subtract off the posterior mean from the data
#' array, then repeat the estimation procedure on the residuals. It
#' continues to do this until the maximum cp-rank \code{k} is reached
#' or when all of the components are estimated to be 0.
#'
#'
#'
#'
#'
#' @inheritParams tflash
#' @param k The maximum cp-rank of the mean tensor.
#'
#'
#' @seealso \code{\link{tflash}} for fitting the rank-1 mean tensor
#'     model.
#'
#' @export
#'
#' @author David Gerard
#'
tgreedy <- function(Y, k = max(dim(Y)), tol = 10^-5, itermax = 100, alpha = 0, beta = 0,
                    mixcompdist = "normal") {
    p <- dim(Y)
    n <- length(p)
    factor_list <- list()
    sig_vec <- c()

    resids <- Y
    rank_index <- 1
    while (rank_index <= k) {
        tflash_out <- tflash(resids, tol = tol, itermax = itermax, alpha = alpha, beta = beta,
                             mixcompdist = mixcompdist)

        sig_vec <- c(sig_vec, tflash_out$sigma_est)
        
        if (sum(abs(tflash_out$postMean[[1]])) < 10^-6) {
            break
        }
        resids <- resids - form_outer(tflash_out$postMean)
        
        for(mode_index in 1:n) {
            if(rank_index == 1) {
                factor_list[[mode_index]] <- tflash_out$postMean[[mode_index]]
            } else {
                factor_list[[mode_index]] <- cbind(factor_list[[mode_index]],
                                                   tflash_out$postMean[[mode_index]])
            }
        }
        rank_index <- rank_index + 1
    }

    rank_final <- unique(sapply(factor_list, ncol))
    return(list(factor_list = factor_list, sig_vec = sig_vec, rank_final = rank_final))
}

#' Perform backfitting starting at output of \code{tgreedy}.
#'
#' @param Y An array of numerics. The data.
#' @param factor_list A list of matrices with the same number of
#'     columns. These are the starting values for the backfitting
#'     algorithm. The intended starting values are can be obtained
#'     from \code{\link{tgreedy}}.
#' @param maxiter_bf A positive integer. The maximum number of
#'     backfitting steps to perform.
#' @param maxiter_vem A positive integer. The maximum number of steps
#'     in each VEM algorithm to perform at each iteration of the
#'     backfitting algorithm.
#' @param mixcompdist What should the mixing distribution be? See
#'     options from the `ashr` package.
#' @param tol_bf A positive numeric. The stopping criterion for the
#'     backfitting algorithm.
#' @param alpha The prior shape parameter for the variance.
#' @param beta The prior rate parameter for the variance.
#' @param sig_vec A vector of positive numerics. The estimates of
#'     variances that were returned by \code{\link{tgreedy}}.
#'
#' @export
#'
#' @author David Gerard
#'
#'
tbackfitting <- function(Y, factor_list, sig_vec, maxiter_bf = 100, tol_bf = 10^-6,
                         maxiter_vem = 100,
                         mixcompdist = "normal", alpha = 0, beta = 0) {
    p <- dim(Y)
    n <- length(p)
    factor_list <- lapply(factor_list, change_to_mat)

    k <- unique(sapply(factor_list, ncol))
    if(length(k) > 1) {
        stop("matrices in factor_list need to have the same number of columns")
    }

    resids <- Y
    for (mode_index in 1:k) {
        resids <- resids - get_kth_tensor(factor_list = factor_list, k = mode_index)
    }

    iter_bf <- 1
    err_bf <- tol_bf + 1
    while (iter_bf < maxiter_bf & err_bf > tol_bf) {
        sig_old <- sig_vec
        
        for (factor_index in 1:k) {
            resids <- resids + get_kth_tensor(factor_list = factor_list, k = factor_index)
            t_out <- tflash(Y = resids, itermax = maxiter_vem, alpha = alpha, beta = beta,
                            mixcompdist = mixcompdist)
            new_factors <- rescale_factors(t_out$postMean)
            resids <- resids - form_outer(new_factors)
            factor_list <- replace_factors(factor_list = factor_list,
                                           new_factors = new_factors,
                                           k = factor_index)
            sig_vec[factor_index] <- t_out$sigma_est
        }
        err_bf <- sum(abs(sig_old[1:k] - sig_vec[1:k]))
        iter_bf <- iter_bf + 1
    }
    return(list(factor_list = factor_list, sig_vec = sig_vec[1:k]))
}


#' Form tensor from a list of matrices of factors.
#' 
#' @param factor_list A list of matrices of numerics, each with the
#'     same number of columns. The kth column of the nth element of
#'     \code{factor_list} is the kth component of the nth mode of mean
#'     tensor.
#' 
#' @export
#' 
#' @author David Gerard
#' 
form_mean <- function(factor_list) {
    factor_list <- lapply(factor_list, change_to_mat)
    
    k <- unique(sapply(factor_list, ncol))
    if(length(k) > 1) {
        stop("matrices in factor_list need to have the same number of columns")
    }
    
    p <- sapply(factor_list, nrow)
    theta <- array(0, dim = p)

    for (mode_index in 1:k) {
        theta <- theta + get_kth_tensor(factor_list = factor_list, k = mode_index)
    }
    
    return(theta)
}

#' Replace the \code{k}th columns in the matrices in \code{factor_list} with the
#' vectors in \code{new_factors}.
#' 
#' @param factor_list A list of matrices of numerics, each with the
#'     same number of columns. The kth column of the nth element of
#'     \code{factor_list} is the kth component of the nth mode of mean
#'     tensor.
#' @param new_factors A list of vectors. The new factors to place in
#'     the \code{k}th column of \code{factor_list}.
#' @param k A positive integer. The column of the matrices in
#'     \code{factor_list} to replace.
#' 
#' @author David Gerard
#' 
replace_factors <- function(factor_list, new_factors, k) {
    for(mode_index in 1:length(factor_list)) {
        factor_list[[mode_index]][, k] <- new_factors[[mode_index]]
    }
    return(factor_list)
}

#' Rescale factors so some don't get too big.
#'
#' @param x A list of vectors.
#'
#' @export
#'
#' @author David Gerard
rescale_factors <- function(x) {
    mults <- sapply(x, FUN = function(x) { sqrt(sum(x^2)) })
    total_mult <- prod(mults)
    for (mode_index in 1:length(x)) {
        x[[mode_index]] <- x[[mode_index]] / mults[mode_index]
    }
    
    y <- lapply(x, FUN = function(x, d) { x * d }, d = total_mult ^ (1 / length(x)))

    return(y)
}

#' Check to see if a vector and if so, change to a matrix.
#'
#' @param x Either a vector or a matrix.
#'
#' @author David Gerard
change_to_mat <- function(x) {
    if(is.vector(x)) {
        x <- matrix(x, ncol = 1)
    }
    return(x)
}

#' Form the rank-1 tensor from the kth components.
#'
#' @param factor_list A list of matrices of numerics, each with the
#'     same number of columns. The kth column of the nth element of
#'     \code{factor_list} is the kth component of the nth mode of mean
#'     tensor.
#' @param k A positive integer. The component to form.
#'
#' @author David Gerard
#'
get_kth_tensor <- function(factor_list, k) {
    return(form_outer(lapply(factor_list, FUN = get_kth_column, k = k)))
}

#' Extract the kth column from a matrix.
#'
#' @param X A matrix.
#' @param k A positive integer. The column to extract.
#'
#' @author David Gerard
get_kth_column <- function(X, k) {
    return(X[, k])
}

## n <- 50
## p <- 50
## Theta <- rnorm(n) %*% t(rnorm(p)) * 1
## E <- matrix(rnorm(n * p), nrow = n)
## Y <- Theta + E


## p <- c(50, 50, 50)
## u <- list()
## u[[1]] <- rnorm(p[1])
## u[[2]] <- rnorm(p[2])
## u[[3]] <- rnorm(p[3])

## v <- list()
## v[[1]] <- rnorm(p[1])
## v[[2]] <- rnorm(p[2])
## v[[3]] <- rnorm(p[3])

## Theta <- form_outer(u) + form_outer(v)
## E <- array(rnorm(prod(p)), dim = p)
## Y <- Theta + E
## t_out <- tgreedy(Y)
## factor_list <- t_out$factor_list
## sig_vec <- t_out$sig_vec

## b_out <- tbackfitting(Y = Y, factor_list = factor_list, sig_vec = sig_vec)

## t_out

## sum((Theta - form_mean(t_out$factor_list))^2)
## sum((Theta - form_mean(b_out$factor_list))^2)
## sum((Theta - Y)^2)

## plot(t_out$factor_list[[1]][,1], v[[1]])
## plot(t_out$factor_list[[2]][,1], v[[2]])
## plot(t_out$factor_list[[3]][,1], v[[3]])
