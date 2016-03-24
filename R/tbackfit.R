#' Iteratively do tflash on the residuals.
#'
#'
#' @inheritParams tflash
#' @param k The maximum cp-rank of the mean tensor.
#'
#'
#'
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

    return(list(factor_list = factor_list, k = k))
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
## Theta <- form_outer(u)
## E <- array(rnorm(prod(p)), dim = p)
## Y <- Theta + E
## t_out <- tgreedy(Y)
