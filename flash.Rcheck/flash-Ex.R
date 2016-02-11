pkgname <- "flash"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('flash')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
cleanEx()
nameEx("flash_VEM")
### * flash_VEM

flush(stderr()); flush(stdout())

### Name: flash_VEM
### Title: Multivariate Adaptive Shrinkage (original version)
### Aliases: flash_VEM

### ** Examples

NULL



### * <FOOTER>
###
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
