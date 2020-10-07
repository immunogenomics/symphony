## Function for building a Symphony reference from a Harmony object.

###### Prior to running this function, ensure that the following steps have been performed:
###### Normalization
# exp_ref <- singlecellmethods::normalizeData(exp_ref, 1e4, 'log')

###### Variable gene selection
# var_genes = vargenes_vst(exp_ref, topn = 2000)

###### Subset gene expression matrix by variable genes
# vargenes_means_sds <- tibble(symbol = var_genes, mean = Matrix::rowMeans(exp_ref))
# vargenes_means_sds$stddev <- singlecellmethods::rowSDs(exp_ref, vargenes_means_sds$mean)
# exp_ref <- exp_ref[var_genes, ]
# vargenes_means_sds <- tibble(symbol = var_genes, mean = Matrix::rowMeans(exp_ref))
# vargenes_means_sds$stddev <- singlecellmethods::rowSDs(exp_ref, vargenes_means_sds$mean)

###### Scaling
# exp_ref_scaled <- singlecellmethods::scaleDataWithStats(exp_ref, vargenes_means_sds$mean, vargenes_means_sds$stddev, 1)

###### PCA
# s = irlba(exp_ref_scaled, nv = 20)

## Run Harmony
buildReferenceFromHarmonyObj <- function(harmony_obj, #output object from HarmonyMatrix()
                           metadata,
                           vargenes_means_sds,
                           pca_loadings,           # genes x PCs
                           verbose = TRUE, 
                           do_umap = TRUE, 
                           save_uwot_path = NULL) {
    
    set.seed(111) # for reproducibility
    
    if (verbose) message('save metadata, vargenes (S), and loadings (U)')
    res <- list(meta_data = metadata)
    res$vargenes = vargenes_means_sds
    res$loadings = pca_loadings
    
    if(verbose) message('Save R, Z_orig, Z_corr, and betas from Harmony object')
    res$R <- harmony_obj$R
    res$Z_orig <- harmony_obj$Z_orig
    res$Z_corr <- harmony_obj$Z_corr
    res$betas <- harmony:::moe_ridge_get_betas(harmony_obj)
    
    if(verbose) message('Calculate final L2 normalized reference centroids (Y_cos)')
    res$centroids <- t(singlecellmethods:::cosine_normalize_cpp(harmony_obj$R %*% t(harmony_obj$Z_corr), 1))
    
    if(verbose) message('Calculate reference compression terms (Nr and C)')
    res$cache <- compute_ref_cache(res$R, res$Z_corr)

    if (do_umap) {
        if (verbose) message('start UMAP')
        res$umap <- uwot::umap(
            t(res$Z_corr), n_neighbors = 30, learning_rate = 0.5, init = "laplacian", 
            metric = 'cosine', fast_sgd = TRUE,
            min_dist = .1, n_threads = 4, ret_model = TRUE
        )
        
        # Since the nn-index component of the uwot model is not able to be saved as an 
        # RDS object, we save the uwot model entirely at a user-defined path in order to
        # later use transform_umap to project the query cells into the umap embedding.
        if (!is.null(save_uwot_path)) {
            model = uwot::save_uwot(res$umap, file = save_uwot_path, unload = FALSE, verbose = FALSE)
            res$save_uwot_path = save_uwot_path
            if (verbose) message(paste('saved uwot model'))
        }
    }
    if (verbose) message('Finished nicely.')
    return(res)
}
