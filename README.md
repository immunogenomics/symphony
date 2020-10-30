
<!-- README.md is generated from README.Rmd. Please edit that file -->

# symphony

<!-- badges: start -->

<!-- badges: end -->

Efficient single-cell reference mapping with Symphony

## Installation

Install the current version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("immunogenomics/symphony")
```

## Reference building

### Option 1: Starting from genes by cells matrix

``` r
library(symphony)
# Build reference
reference = symphony::buildReference(
    ref_exp,                 # reference cell genes by cells matrix
    ref_metadata,            # reference metadata
    vars = c('donor'),       # variable(s) to integrate over
    K = 100,                 # number of Harmony clusters
    verbose = TRUE,          # show output
    do_umap = TRUE,          # run UMAP and save UMAP model
    do_normalize = FALSE,    # within-cell normalize the expression matrix
    vargenes_method = 'vst', # 'vst' or 'mvp'
    topn = 2000,             # number of variable genes to use
    d = 20,                  # number of PCs
    save_uwot_path = '/absolute/path/uwot_model_1' # filepath to save UMAP model
)
```

### Option 2: Starting from existing Harmony object

TBD

What is special about using `README.Rmd` instead of just `README.md`?
You can include R chunks like so:

``` r
summary(cars)
#>      speed           dist       
#>  Min.   : 4.0   Min.   :  2.00  
#>  1st Qu.:12.0   1st Qu.: 26.00  
#>  Median :15.0   Median : 36.00  
#>  Mean   :15.4   Mean   : 42.98  
#>  3rd Qu.:19.0   3rd Qu.: 56.00  
#>  Max.   :25.0   Max.   :120.00
```

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date.

You can also embed plots, for example:

<img src="man/figures/README-pressure-1.png" width="100%" />

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub\!
