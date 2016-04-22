library(flashr)
context("Kronecker Structured Variance")

test_that("Fitting a Kronecker structured variance model doesn't throw an error",{
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

    tout <- tflash_kron(Y = Y)
}
)

test_that("diag_mle returns homoscedastic modes when specified", {
    p <- c(11, 13, 17)
    R <- array(rnorm(prod(p)), dim = p)
    
    mout <- diag_mle(R, homo_modes = 1)
    expect_equal(mout[[1]] / mout[[1]][1], rep(1, length = p[1]))
    
    mout <- diag_mle(R, homo_modes = 2)
    expect_equal(mout[[2]] / mout[[2]][1], rep(1, length = p[2]))
    
    mout <- diag_mle(R, homo_modes = 3)
    expect_equal(mout[[3]] / mout[[3]][1], rep(1, length = p[3]))

    mout <- diag_mle(R, homo_modes = c(1, 3))
    expect_equal(mout[[1]] / mout[[1]][1], rep(1, length = p[1]))
    expect_equal(mout[[3]] / mout[[3]][1], rep(1, length = p[3]))
    
    mout <- diag_mle(R, homo_modes = c(1, 2))
    expect_equal(mout[[1]] / mout[[1]][1], rep(1, length = p[1]))
    expect_equal(mout[[2]] / mout[[2]][1], rep(1, length = p[2]))
    
    mout <- diag_mle(R, homo_modes = c(2, 3))
    expect_equal(mout[[2]] / mout[[2]][1], rep(1, length = p[2]))
    expect_equal(mout[[3]] / mout[[3]][1], rep(1, length = p[3]))
}
)
