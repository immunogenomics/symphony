#' Function for building a Symphony reference from a Harmony object. Useful if you would like your
#' code to be more modular. Note that you must have saved vargenes_means_sds and PCA loadings.
#'
#' @param harmony_obj Harmony object (output from HarmonyMatrix())
#' @param metadata Reference cell metadata (cells by attributes)
#' @param vargenes_means_sds Variable genes in dataframe with columns named ('symbol', 'mean', 'stddev')
#' @param pca_loadings Gene loadings from PCA (e.g. irlba(ref_exp_scaled, nv = 20)$u)
#' @param verbose Verbose output
#' @param do_umap Perform UMAP visualization on harmonized reference embedding
#' @param save_uwot_path Absolute path to save the uwot model (if do_umap is TRUE)
#' @param umap_min_dist UMAP parameter (see uwot documentation for details)
#' @param seed Random seed
#' @return Symphony reference object. Integrated embedding is stored in the $Z_corr slot. Other slots include
#' cell-level metadata ($meta_data), variable genes means and standard deviations ($vargenes),
#' loadings from PCA or other dimensional reduction such as CCA ($loadings), original PCA embedding ($Z_orig), 
#' reference compression terms ($cache), betas from Harmony integration ($betas), cosine-normalized soft cluster centroids ($centroids), 
#' centroids in PC space ($centroids_pc), and optional umap coordinates ($umap$embedding).
#' @export
buildReferenceFromHarmonyObj <- function(harmony_obj,
                           metadata,
                           vargenes_means_sds,
                           pca_loadings,           # genes x PCs
                           verbose = TRUE, 
                           do_umap = TRUE, 
                           save_uwot_path = NULL,
                           umap_min_dist = 0.1,
                           seed = 111) {
    
    set.seed(seed) # for reproducibility
    
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

    # Add row and column names
    colnames(res$Z_orig) = row.names(metadata)
    rownames(res$Z_orig) = paste0("PC_", seq_len(nrow(res$Z_corr)))
    colnames(res$Z_corr) = row.names(metadata)
    rownames(res$Z_corr) = paste0("harmony_", seq_len(nrow(res$Z_corr)))
    
    # Compute centroids in harmony PC space
    cluster_sizes = res$cache[[1]] %>% as.matrix()
    centroid_sums = t(res$Z_corr %*% t(res$R)) %>% as.data.frame()
    centroids_pc = sweep(centroid_sums, 1, cluster_sizes, "/")
    colnames(centroids_pc) = paste0("harmony_", seq_len(nrow(res$Z_corr)))
    rownames(centroids_pc) = paste0("centroid_", seq_len(nrow(res$R)))
    res$centroids_pc = centroids_pc
    
    if (do_umap) {
        if (verbose) message('UMAP')
        umap = uwot::umap(
            t(res$Z_corr), n_neighbors = 30, learning_rate = 0.5, init = "laplacian", 
            metric = 'cosine', fast_sgd = FALSE, n_sgd_threads = 1, # for reproducibility
            min_dist = umap_min_dist, n_threads = 4, ret_model = TRUE
        )
        res$umap$embedding = umap$embedding
        colnames(res$umap$embedding) = c('UMAP1', 'UMAP2')
        
        # Since the nn-index component of the uwot model is not able to be saved as an 
        # object, we save the uwot model at a user-defined path.
        if (!is.null(save_uwot_path)) {
            
            # If file already exists, delete it (otherwise will result in an error)
            if (file.exists(save_uwot_path)) {
                if (verbose) message(paste('File already exists at that path... overwriting...'))
                file.remove(save_uwot_path)
            }
            
            model = uwot::save_uwot(umap, file = save_uwot_path, unload = FALSE, verbose = FALSE)
            res$save_uwot_path = save_uwot_path
            if (verbose) message(paste('Saved uwot model'))
        }
    }
    if (verbose) message('Finished nicely.')
    return(res)
}
