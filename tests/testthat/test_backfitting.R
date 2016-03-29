## library(flashr)
## context("Backfitting")

## test_that("backfitting produces same results for vector as matrix", {
##     set.seed(10)
##     l = rnorm(10)
##     f=rnorm(5)
##     Y = l %*% t(f) + rnorm(10*5)
##     temp=backfitting(Y,l,f)
##     temp2 = backfitting(Y,matrix(l,ncol=1),matrix(f,ncol=1))
##     expect_equal(temp,temp2)
## })

## test_that("backfitting works when results goes to 0", {
##     set.seed(10)
##     l = rnorm(10)
##     f=rnorm(5)
##     Y = matrix(rnorm(10*5),nrow=10)
##     temp = backfitting(Y,l,f)
##     expect_equal(temp$l,matrix(0,nrow=10,ncol=0))
## })
