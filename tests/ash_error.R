

### ash error

ll <- get(load("ash_error_list.rda"))

m <- ll$m
betahat <- ll$betahat
sebetahat <- ll$sebetahat

ashr::ash(betahat, sebetahat)

library(ashr)
comp_postmean(m, betahat, sebetahat)
