library(flashr)
context("tflash vs flash")

test_that("known factors version of tflash and flash give equivalent results",{
    set.seed(211)
    n <- 10
    p <- 100
    k <- 5
    q <- 1
    
    pi_vals <- c(0.5, 0.5)
    tau_seq <- c(0, 1)
    
    X <- matrix(rnorm(n * 1), nrow = n)
    beta <- succotashr::draw_beta(pi_vals = pi_vals, tau_seq = tau_seq, p = p)
    E <- matrix(rnorm(n * p), nrow = n)
    Y <- X %*% t(beta) + E
    
    fout <- flash(Y = t(Y), factor_value = c(X), fix_factor = TRUE)
    tout <- tflash(Y = Y, known_factors = list(X), known_modes = 1)
    
    expect_equal(fout$l, tout$post_mean[[2]], tolerance = 10 ^ -2)
}
)


test_that("greedy and tgreedy give same results", {
    set.seed(89)
    n <- 100
    p <- 100
    r <- 2
    A <- matrix(rnorm(r * n), nrow = n)
    B <- matrix(rnorm(r * p), ncol = p)
    Theta <- A %*% B
    E <- matrix(rnorm(n * p), nrow = n)
    Y <- Theta + E

    trash <- capture.output(fout <- flashr::greedy(Y, K = 10))
    tout <- flashr::tgreedy(Y, k = 10)

    fl1scale <- fout$l[,1] / sqrt(sum(fout$l[,1] ^ 2))
    fl2scale <- fout$l[,2] / sqrt(sum(fout$l[,2] ^ 2))
    ff1scale <- fout$f[,1] / sqrt(sum(fout$f[,1] ^ 2))
    ff2scale <- fout$f[,2] / sqrt(sum(fout$f[,2] ^ 2))

    tl1scale <- tout$factor_list[[1]][,1] / sqrt(sum(tout$factor_list[[1]][,1] ^ 2))
    tl2scale <- tout$factor_list[[1]][,2] / sqrt(sum(tout$factor_list[[1]][,2] ^ 2))
    tf1scale <- tout$factor_list[[2]][,1] / sqrt(sum(tout$factor_list[[2]][,1] ^ 2))
    tf2scale <- tout$factor_list[[2]][,2] / sqrt(sum(tout$factor_list[[2]][,2] ^ 2))
    
    expect_equal(abs(fl1scale), abs(tl1scale), tolerance = 10 ^ -2)
    expect_equal(abs(fl2scale), abs(tl2scale), tolerance = 10 ^ -2)
    expect_equal(abs(ff1scale), abs(tf1scale), tolerance = 10 ^ -2)
    expect_equal(abs(ff2scale), abs(tf2scale), tolerance = 10 ^ -2)
    
}
)


test_that("backfitting is same between flash and tflash", {
    skip("bugs in backfitting")

    set.seed(9489)
    n <- 10
    p <- 20
    r <- 2
    A <- matrix(rnorm(r * n), nrow = n)
    B <- matrix(rnorm(r * p), ncol = p)
    Theta <- A %*% B
    E <- matrix(rnorm(n * p), nrow = n)
    Y <- Theta + E
    
    foutg <- flashr::greedy(Y, K = 10)
    fout <- flashr::backfitting(Y, intial_list = list(l = matrix(foutg$l), f = matrix(foutg$f)))
    toutg <- flashr::tgreedy(Y, k = 10)
    tout <- flashr::tbackfitting(Y = Y, factor_list = toutg$factor_list, sigma_est = toutg$sig_vec)

    fl1scale <- fout$l[,1] / sqrt(sum(fout$l[,1] ^ 2))
    fl2scale <- fout$l[,2] / sqrt(sum(fout$l[,2] ^ 2))
    ff1scale <- fout$f[,1] / sqrt(sum(fout$f[,1] ^ 2))
    ff2scale <- fout$f[,2] / sqrt(sum(fout$f[,2] ^ 2))

    tl1scale <- tout$factor_list[[1]][,1] / sqrt(sum(tout$factor_list[[1]][,1] ^ 2))
    tl2scale <- tout$factor_list[[1]][,2] / sqrt(sum(tout$factor_list[[1]][,2] ^ 2))
    tf1scale <- tout$factor_list[[2]][,1] / sqrt(sum(tout$factor_list[[2]][,1] ^ 2))
    tf2scale <- tout$factor_list[[2]][,2] / sqrt(sum(tout$factor_list[[2]][,2] ^ 2))
    
    expect_equal(abs(fl1scale), abs(tl1scale), tolerance = 10 ^ -2)
    expect_equal(abs(fl2scale), abs(tl2scale), tolerance = 10 ^ -2)
    expect_equal(abs(ff1scale), abs(tf1scale), tolerance = 10 ^ -2)
    expect_equal(abs(ff2scale), abs(tf2scale), tolerance = 10 ^ -2)
}
)
