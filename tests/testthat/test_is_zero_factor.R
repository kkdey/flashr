library(flashr)
context("Zero Factors")

test_that("is_zero_factor returns correct result",{
    l=matrix(rnorm(500),nrow=100,ncol=5)
    l[,1]=rep(0,100)
    l[,3]=rep(0,100)
    expect_equal(is_zero_factor(l),c(TRUE,FALSE,TRUE,FALSE,FALSE))
}
)
