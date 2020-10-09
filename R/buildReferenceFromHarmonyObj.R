#' Function for building a Symphony reference from a Harmony object. Useful if you would like your
#' code to be more modular. Note that you must have your saved vargenes_means_sds and PCA loadings.
#'
#' @param harmony_obj Harmony object (output from HarmonyMatrix())
#' @param metadata Reference cell metadata (cells by attributes)
#' @param vargenes_means_sds Variable genes in dataframe with columns ('symbol', 'mean', 'stddev')
#' @param pca_loadings Gene loadings from PCA (= irlba(ref_exp_scaled, nv = 20)$u)
#' @param verbose Verbose output
#' @param do_umap Perform UMAP visualization on harmonized reference embedding
#' @param save_uwot_path Absolute path to save the uwot model (if do_umap is TRUE)
#' 
#' @export
buildReferenceFromHarmonyObj <- function(harmony_obj,
                           metadata,
                           vargenes_means_sds,
                           pca_loadings,           # genes x PCs
                           verbose = TRUE, 
                           do_umap = TRUE, 
                           save_uwot_path = NULL) {
    
    set.seed(111) # for reproducibility
    
    if (verbose) message('Save metadata, vargenes (S), and loadings (U)')
    res = list(meta_data = metadata)
    res$vargenes = vargenes_means_sds
    res$loadings = pca_loadings
    
    if(verbose) message('Save R, Z_orig, Z_corr, and betas from Harmony object')
    res$R = harmony_obj$R
    res$Z_orig = harmony_obj$Z_orig
    res$Z_corr = harmony_obj$Z_corr
    res$betas = harmony::moe_ridge_get_betas(harmony_obj)
    
    if(verbose) message('Calculate final L2 normalized reference centroids (Y_cos)')
    res$centroids = t(cosine_normalize_cpp(harmony_obj$R %*% t(harmony_obj$Z_corr), 1))
    
    if(verbose) message('Calculate reference compression terms (Nr and C)')
    res$cache = compute_ref_cache(res$R, res$Z_corr)

    if (do_umap) {
        if (verbose) message('UMAP')
        res$umap = uwot::umap(
            t(res$Z_corr), n_neighbors = 30, learning_rate = 0.5, init = "laplacian", 
            metric = 'cosine', fast_sgd = FALSE, n_sgd_threads = 1, # for reproducibility
            min_dist = .1, n_threads = 4, ret_model = TRUE
        )
        
        # Since the nn-index component of the uwot model is not able to be saved as an 
        # object, we save the uwot model at a user-defined path.
        if (!is.null(save_uwot_path)) { # TODO: check if valid file path
            model = uwot::save_uwot(res$umap, file = save_uwot_path, unload = FALSE, verbose = FALSE)
            res$save_uwot_path = save_uwot_path
            if (verbose) message(paste('Saved uwot model'))
        }
    }
    if (verbose) message('Finished nicely.')
    return(res)
}
