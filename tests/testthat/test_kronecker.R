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
