
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Symphony

<!-- badges: start -->

<!-- badges: end -->

Efficient and precise single-cell reference atlas mapping with Symphony

# Installation

Install the current version of Symphony from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("immunogenomics/symphony")
```

### Installation notes:

  - You may need to install the latest version of devtools (because of
    the recent GitHub change from “master” to “main” terminology, which
    can cause `install_github` to fail).
  - You may also need to install the lastest version of Harmony:

<!-- end list -->

``` r
devtools::install_github("immunogenomics/harmony")
```

# Usage/Demos

## Quick start

Check out the [quick start
tutorial](https://github.com/immunogenomics/symphony/blob/main/vignettes/pbmcs_tutorial.ipynb/).

## Reference building

### Option 1: Starting from reference genes by cells matrix

This function performs all steps of the reference building pipeline
including variable gene selection, scaling, PCA, Harmony, and Symphony
compression.

``` r
# Build reference
reference = symphony::buildReference(
    ref_exp,                 # reference genes by cells matrix
    ref_metadata,            # dataframe with cell metadata
    vars = c('donor'),       # variable(s) to integrate over
    K = 100,                 # number of Harmony clusters
    verbose = TRUE,          # display output?
    do_umap = TRUE,          # run UMAP and save UMAP model to file?
    do_normalize = FALSE,    # normalize the expression matrix?
    vargenes_method = 'vst', # 'vst' or 'mvp'
    topn = 2000,             # number of variable genes to use
    d = 20,                  # number of dimensions for PCA
    save_uwot_path = '/absolute/path/uwot_model_1' # filepath to save UMAP model
)
```

### Option 2: Starting from existing Harmony object

This function compresses an existing Harmony object into a Symphony
reference that enables query mapping. We recommend this option if you
would like your code to be more modular and flexible.

``` r

# Run Harmony to integrate the reference cells
ref_harmObj = harmony::HarmonyMatrix(
        data_mat = t(Z_pca_ref),   # starting embedding (e.g. PCA, CCA) of cells
        meta_data = ref_metadata,  # dataframe with cell metadata
        theta = c(2),              # cluster diversity enforcement
        vars_use = c('donor'),     # variable to integrate out
        nclust = 100,              # number of clusters in Harmony model
        max.iter.harmony = 10,
        return_object = TRUE,      # set to TRUE to return the full Harmony object
        do_pca = FALSE             # do not recompute PCs
)

# Build Symphony reference
reference = buildReferenceFromHarmonyObj(
        ref_harmObj,            # output object from HarmonyMatrix()
        ref_metadata,           # dataframe with cell metadata
        vargenes_means_sds,     # gene names, means, and std devs for scaling
        loadings,               # genes x PCs
        verbose = TRUE,         # display output?
        do_umap = TRUE,         # run UMAP and save UMAP model to file?
        save_uwot_path = '/absolute/path/uwot_model_1' # filepath to save UMAP model)
```

Note that `vargenes_means_sds` requires column names `c('symbol',
'mean', 'stddev')` (see [tutorial
example](https://github.com/immunogenomics/symphony/blob/main/vignettes/pbmcs_tutorial.ipynb/)).

## Query mapping

Once you have a prebuilt reference (e.g. loaded from a saved .rds R
object), you can map new query cells onto it starting from query gene
expression.

``` r
# Map query
query = mapQuery(query_exp, query_metadata, reference, do_normalize = FALSE)
```

`query$Z` contains the harmonized query feature embedding.

If your query itself has multiple sources of batch variation you would
like to integrate over (e.g. technology, donors, species), you can
specify them in the `vars` parameter.

``` r
# Map query
query = mapQuery(query_exp, query_metadata, vars = c('donor', 'technology') reference, do_normalize = FALSE)
```

# Reproducing results from manuscript

Code to reproduce Symphony results from the Kang et al. manuscript will
be made available on github.com/immunogenomics/referencemapping.
