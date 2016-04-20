library(flashr)
context("Known Factors")

test_that("tflash returns known factors when given known factors",{
    set.seed(231)
    p <- c(11, 13, 17)

    E <- array(rnorm(prod(p)), dim = p)
    X <- rnorm(p[1])
    beta1 <- rnorm(p[2])
    beta2 <- rnorm(p[3])
    Y <- outer(outer(X, beta1, "*"), beta2, "*") + E

    tout <- tflash(Y = Y, known_factors = list(X), known_modes = 1)
    expect_equal(X, tout$post_mean[[1]])

    tout <- tflash(Y = Y, known_factors = list(X, beta1), known_modes = c(1, 2))
    expect_equal(X, tout$post_mean[[1]])
    expect_equal(beta1, tout$post_mean[[2]])
}
)

test_that("tgreedy returns known factors when given known factors",{
    set.seed(101)
    n <- 11
    p <- 21

    E <- matrix(rnorm(n * p), nrow = n)
    X <- matrix(rnorm(n * 2), nrow = n)
    beta <- matrix(rnorm(p * 2), ncol = 2)
    Y <- X %*% t(beta) + E

    tout <- tgreedy(Y = Y, known_factors = list(X), known_modes = 1)

    expect_equal(X[, 1], tout$factor_list[[1]][, 1])
    expect_equal(X[, 2], tout$factor_list[[1]][, 2])
}
)


test_that("tbackfit returns known factors when given known factors",{
    set.seed(618)
    n <- 11
    p <- 21

    E <- matrix(rnorm(n * p), nrow = n)
    X <- matrix(rnorm(n * 2), nrow = n)
    beta <- matrix(rnorm(p * 2), ncol = 2)
    Y <- X %*% t(beta) + E

    gtout <- tgreedy(Y = Y, known_factors = list(X), known_modes = 1)
    tout <- tbackfitting(Y = Y, factor_list = gtout$factor_list,
                         sigma_est = gtout$sigma_est, known_factors = list(X),
                         known_modes = 1)
    expect_equal(X[, 1], tout$factor_list[[1]][, 1])
    expect_equal(X[, 2], tout$factor_list[[1]][, 2])
}
)

test_that("different number of known factors in each mode works ok",{
    set.seed(618)
    p <- c(11, 13, 17)
    u <- list()
    u[[1]] <- rnorm(p[1])
    u[[2]] <- rnorm(p[2])
    u[[3]] <- rnorm(p[3])
    v <- list()
    v[[1]] <- rnorm(p[1])
    v[[2]] <- rnorm(p[2])
    v[[3]] <- rnorm(p[3])
    
    Theta <- form_outer(u) + form_outer(v)
    E <- array(rnorm(prod(p)), dim = p)
    Y <- Theta + E
    
    known_factors = list(cbind(u[[1]], v[[1]]), u[[2]])
    tout <- tgreedy(Y = Y, known_factors = known_factors, known_modes = c(1, 2))
    expect_equal(tout$factor_list[[1]][, 1], u[[1]])
    expect_equal(tout$factor_list[[2]][, 1], u[[2]])
    expect_equal(tout$factor_list[[1]][, 2], v[[1]])
}
)
