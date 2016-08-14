
Repo for flashr R package: Factor Loading Adaptive SHrinkage for matrix and tensor data.

## Installation


To install dependencies, run in R:

``` r
install.packages(c("irlba", "tensr", "devtools"))
devtools::install_github("stephens999/ashr", ref = "uni")
```

The "uni" branch of the ashr package is required to use any mixing distribution other than a normal.

Because this is currently a private repo, to install flashr you will need a private access token (PAT) which you can generate here: <https://github.com/settings/tokens>

Then you can run in R:

``` r
devtools::install_github("stephenslab/flashr", auth_token = "xxx")
```

where you replace `xxx` with your PAT.

To get help ??flashr

```
# flashr
```

Factor Loading with Adaptive Shrinkage in R
