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
#' @return \code{factor_list}: A list of matrices of
#'     numerics. \code{factor_list[[i]][, j]} contains the \eqn{j}th
#'     factor of the \eqn{i}th mode.
#'
#'     \code{sigma_est}: If \code{var_type = "homoscedastic"}, then
#'     \code{sigma_est} is a vector of numerics. \code{sigma_est[i]}
#'     is the estimate of the precision during the \eqn{i} iteration
#'     of the greedy algorithm. Only the last one (if that) should
#'     actually be used for any sort of precision estimate.
#'
#'     If \code{var_type = "kronecker"}, then \code{sigma_est} is a
#'     list of matices. \code{sigma_est[[i]][, j]} is the estimate of
#'     the variances for the \eqn{j}th mode during the \eqn{i}th run
#'     of the greedy algorithm. Only the final columns (if those) of
#'     these matrices should actually be used as any sort of precision
#'     estimate.
#'
#'     \code{rank_final} A non-negative integer. The final estimated
#'     cp-rank of the mean.
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
                    mixcompdist = "normal", var_type = c("homoscedastic", "kronecker"),
                    nullweight = 1, print_update = FALSE, known_factors = NULL,
                    known_modes = NULL) {
    p <- dim(Y)
    n <- length(p)
    factor_list <- list()
    prob_zero_list <- list()
    pi0_list <- list()

    var_type <- match.arg(var_type, c("homoscedastic", "kronecker"))

    if (!is.null(known_factors)) {
        known_factors <-lapply(known_factors, change_to_mat)
    }

    ## checks
    if (!is.null(known_modes)) {
        if (is.null(known_factors)) {
            stop("known_modes is not NULL but known_factors is NULL")
        } else if(!all(sapply(known_factors, is.matrix))) {
            stop("all known_factors must be a matrix")
        }
        dim_factors <- sapply(known_factors, nrow)
        num_known_factors <- sapply(known_factors, ncol)
        if (length(known_modes) != length(known_factors)) {
            stop ("known_modes and known_factors must be of same length")
        } else if (any(p[known_modes] != dim_factors)) {
            stop("known_factors not the same dimension as modes of Y")
        }
    }



    sigma_est <- c()

    resids <- Y
    rank_index <- 1
    while (rank_index <= k) {
        if(print_update) {
            cat("Current Rank:", rank_index, "\n\n")
        }

        ## Extract known factors
        if (!is.null(known_factors)) {
            which_now_known <- num_known_factors >= rank_index
            if(any(which_now_known)) {
                current_known_modes <- known_modes[which_now_known]
                current_known_factors <- vector(length = length(current_known_modes),
                                                mode = "list")
                wnk_tracker <- 1
                for (known_mode_index in 1:length(known_modes)) {
                    if (which_now_known[known_mode_index]) {
                        current_known_factors[[wnk_tracker]] <-
                            known_factors[[known_mode_index]][, rank_index]
                        wnk_tracker <- wnk_tracker + 1
                    }
                }
            } else {
                current_known_factors <- NULL
                current_known_modes <- NULL
            }
        } else {
            current_known_factors <- NULL
            current_known_modes <- NULL
        }


        if(var_type == "homoscedastic") {
            tflash_out <- tflash_homo(resids, tol = tol, itermax = itermax, alpha = alpha,
                                      beta = beta, mixcompdist = mixcompdist,
                                      nullweight = nullweight, print_update = print_update,
                                      known_factors = current_known_factors,
                                      known_modes = current_known_modes)

            sigma_est <- c(sigma_est, tflash_out$sigma_est)
        } else if (var_type == "kronecker") {
            tflash_out <- tflash_kron(resids, tol = tol, itermax = itermax, alpha = alpha,
                                      beta = beta, mixcompdist = mixcompdist,
                                      nullweight = nullweight, print_update = print_update,
                                      known_factors = current_known_factors,
                                      known_modes = current_known_modes)
            for(mode_index in 1:n) {
                if(rank_index > 1) {
                    sigma_est[[mode_index]] <- cbind(sigma_est[[mode_index]],
                                                     tflash_out$sigma_est[[mode_index]])
                } else {
                    sigma_est[[mode_index]] <- tflash_out$sigma_est[[mode_index]]
                }
            }
        }

        if (sum(abs(tflash_out$post_mean[[1]])) < 10^-6) {
            break
        }
        resids <- resids - form_outer(tflash_out$post_mean)

        for(mode_index in 1:n) {
            if(rank_index == 1) {
                factor_list[[mode_index]] <- tflash_out$post_mean[[mode_index]]
                if (is.null(tflash_out$prob_zero[[mode_index]])) {
                    prob_zero_list[[mode_index]] <-
                        rep(NA, length = length(tflash_out$post_mean[[mode_index]]))
                    pi0_list[[mode_index]] <- NA
                } else {
                    prob_zero_list[[mode_index]] <- tflash_out$prob_zero[[mode_index]]
                    pi0_list[[mode_index]] <- tflash_out$pi0vec[mode_index]
                }
            } else {
                factor_list[[mode_index]] <- cbind(factor_list[[mode_index]],
                                                   tflash_out$post_mean[[mode_index]])
                if (is.null(tflash_out$prob_zero[[mode_index]])) {
                    prob_zero_list[[mode_index]] <-
                        cbind(prob_zero_list[[mode_index]],
                              rep(NA, length = length(tflash_out$post_mean[[mode_index]])))
                    pi0_list[[mode_index]] <- c(pi0_list[[mode_index]], NA)
                } else {
                    prob_zero_list[[mode_index]] <-
                        cbind(prob_zero_list[[mode_index]],
                              tflash_out$prob_zero[[mode_index]])
                    pi0_list[[mode_index]] <- c(pi0_list[[mode_index]],
                                                tflash_out$pi0vec[mode_index])
                }
            }
        }
        rank_index <- rank_index + 1

    }

    rank_final <- unique(sapply(factor_list, ncol))
    return(list(factor_list = factor_list, prob_zero_list = prob_zero_list,
                pi0_list = pi0_list, sigma_est = sigma_est, rank_final = rank_final))
}

#' Perform backfitting starting at output of \code{tgreedy}.
#'
#' @param factor_list A list of matrices with the same number of
#'     columns. These are the starting values for the backfitting
#'     algorithm. The intended starting values are can be obtained
#'     from \code{\link{tgreedy}}.
#' @param maxiter_bf A positive integer. The maximum number of
#'     backfitting steps to perform.
#' @param maxiter_vem A positive integer. The maximum number of steps
#'     in each VEM algorithm to perform at each iteration of the
#'     backfitting algorithm.
#' @param tol_bf A positive numeric. The stopping criterion for the
#'     backfitting algorithm.
#' @param sigma_est Either a vector of estimated precisions (when
#'     \code{var_type = "homoscedastic"}) or a list of matrices whose
#'     columns are estimated precisions (when \code{var_type =
#'     "kronecker"}).
#' @param tol_r1 A positive numeric. The tolerance for the rank 1 runs
#'     of tflash.
#' 
#' @inheritParams tflash
#'
#' @export
#'
#' @author David Gerard
#'
#'
tbackfitting <- function(Y, factor_list, sigma_est, maxiter_bf = 100, tol_bf = 10^-6,
                         maxiter_vem = 100, var_type = c("homoscedastic", "kronecker"),
                         mixcompdist = "normal", alpha = 0, beta = 0, nullweight = 10,
                         known_factors = NULL, known_modes = NULL, tol_r1 = 10 ^ -3) {
    p <- dim(Y)
    n <- length(p)
    factor_list <- lapply(factor_list, change_to_mat)
    prob_zero_list <- list()
    pi0_list <- list()
    for (mode_index in 1:n) {
        prob_zero_list[[mode_index]] <-
            matrix(NA, nrow = nrow(factor_list[[mode_index]]),
                   ncol = ncol(factor_list[[mode_index]]))
        pi0_list[[mode_index]] <- rep(NA, length = ncol(factor_list[[mode_index]]))
    }

    k <- unique(sapply(factor_list, ncol))
    if(length(k) > 1) {
        stop("matrices in factor_list need to have the same number of columns")
    }

    var_type <- match.arg(var_type, c("homoscedastic", "kronecker"))

    ## checks
    if (!is.null(known_modes)) {
        if (is.null(known_factors)) {
            stop("known_modes is not NULL but known_factors is NULL")
        } else if(!all(sapply(known_factors, is.matrix))) {
            stop("all known_factors must be a matrix")
        }
        dim_factors <- sapply(known_factors, nrow)
        num_known_factors <- sapply(known_factors, ncol)
        if (length(known_modes) != length(known_factors)) {
            stop ("known_modes and known_factors must be of same length")
        } else if (any(p[known_modes] != dim_factors)) {
            stop("known_factors not the same dimension as modes of Y")
        }
    }

    resids <- Y
    for (factor_index in 1:k) {
        resids <- resids - get_kth_tensor(factor_list = factor_list, k = factor_index)
    }

    iter_bf <- 1
    err_bf <- tol_bf + 1
    while (iter_bf < maxiter_bf & err_bf > tol_bf) {
        sig_old <- sigma_est

        for (factor_index in 1:k) {
            resids <- resids + get_kth_tensor(factor_list = factor_list, k = factor_index)


            ## Extract known factors
            if (!is.null(known_factors)) {
                which_now_known <- num_known_factors >= factor_index
                if(any(which_now_known)) {
                    current_known_modes <- known_modes[which_now_known]
                    current_known_factors <- vector(length = length(current_known_modes),
                                                    mode = "list")
                    wnk_tracker <- 1
                    for (known_mode_index in 1:length(known_modes)) {
                        if (which_now_known[known_mode_index]) {
                            current_known_factors[[wnk_tracker]] <-
                                known_factors[[known_mode_index]][, factor_index]
                            wnk_tracker <- wnk_tracker + 1
                        }
                    }
                } else {
                    current_known_factors <- NULL
                    current_known_modes <- NULL
                }
            } else {
                current_known_factors <- NULL
                current_known_modes <- NULL
            }


            if (var_type == "homoscedastic") {
                t_out <- tflash_homo(Y = resids, itermax = maxiter_vem, alpha = alpha, beta = beta,
                                     mixcompdist = mixcompdist, nullweight = nullweight,
                                     known_factors = current_known_factors,
                                     known_modes = current_known_modes,
                                     tol = tol_r1)
                sigma_est[factor_index] <- t_out$sigma_est
            } else if (var_type == "kronecker") {
                t_out <- tflash_kron(Y = resids, itermax = maxiter_vem, alpha = alpha, beta = beta,
                                     mixcompdist = mixcompdist, nullweight = nullweight,
                                     known_factors = current_known_factors,
                                     known_modes = current_known_modes,
                                     tol = tol_r1)
                sigma_est <- replace_factors(factor_list = sigma_est,
                                             new_factors = t_out$sigma_est,
                                             k = factor_index)
            }

            any_zero <- sapply(t_out$post_mean, function(x) { all(abs(x) < 10 ^ -6) })
            if (any(any_zero)) {
                new_factors <- list()
                for(nf_index in 1:n) {
                  new_factors[[nf_index]] <- rep(0, length = p[nf_index])
                }
            } else if (is.null(current_known_factors)) {
                new_factors <- rescale_factors(t_out$post_mean)
            } else {
                new_factors <- t_out$post_mean
            }
            new_probzero <- list()
            for (pz_index in 1:n) {
                if (is.null(t_out$prob_zero[[pz_index]])) {
                    new_probzero[[pz_index]] <-
                        rep(NA, length = length(t_out$post_mean[[pz_index]]))
                } else
                {
                    new_probzero[[pz_index]] <- t_out$prob_zero[[pz_index]]
                }
            }


            resids <- resids - form_outer(new_factors)
            factor_list <- replace_factors(factor_list = factor_list,
                                           new_factors = new_factors,
                                           k = factor_index)
            prob_zero_list <- replace_factors(factor_list = prob_zero_list,
                                              new_factors = new_probzero,
                                              k = factor_index)
            for (p0index in 1:n) {
                pi0_list[[p0index]][factor_index] <- t_out$pi0vec[p0index]
            }


        }

        if (var_type == "homoscedastic") {
            err_bf <- sum(abs(sig_old[1:k] - sigma_est[1:k]))
        } else if (var_type == "kronecker") {
            err_bf <- 0
            for (mode_index in 1:n) {
                err_bf <- sum(abs(sigma_est[[mode_index]] - sig_old[[mode_index]]))
            }
        }
        message(paste("Diff:", err_bf))
        iter_bf <- iter_bf + 1
    }

    return(list(factor_list = factor_list, sigma_est = sigma_est,
                prob_zero_list = prob_zero_list,
                pi0_list = pi0_list))
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


## p <- c(10, 10, 10)
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
## t_out <- tgreedy(Y, var_type = "kronecker")
## factor_list <- t_out$factor_list
## sigma_est <- t_out$sigma_est

## b_out <- tbackfitting(Y = Y, factor_list = factor_list, sigma_est = sigma_est, var_type = "kronecker")

## sum((Theta - form_mean(t_out$factor_list))^2)
## sum((Theta - form_mean(b_out$factor_list))^2)
## sum((Theta - Y)^2)

## plot(t_out$factor_list[[1]][,1], v[[1]])
## plot(t_out$factor_list[[2]][,1], v[[2]])
## plot(t_out$factor_list[[3]][,1], v[[3]])

## plot(t_out$factor_list[[1]][,2], u[[1]])
## plot(t_out$factor_list[[2]][,2], u[[2]])
## plot(t_out$factor_list[[3]][,2], u[[3]])

## plot(b_out$factor_list[[1]][,1], v[[1]])
## plot(b_out$factor_list[[2]][,1], v[[2]])
## plot(b_out$factor_list[[3]][,1], v[[3]])

## plot(b_out$factor_list[[1]][,2], u[[1]])
## plot(b_out$factor_list[[2]][,2], u[[2]])
## plot(b_out$factor_list[[3]][,2], u[[3]])

