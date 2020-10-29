#' Function for building a Symphony reference starting from expression matrix
#'
#' @param exp_ref Reference gene expression (genes by cells)
#' @param metadata_ref Reference cell metadata (cells by attributes)
#' @param vars Reference variables to Harmonize over e.g. c('donor', 'technology')
#' @param K Number of centroids
#' @param verbose Verbose output
#' @param do_umap Perform UMAP visualization on harmonized reference embedding
#' @param do_normalize Perform log(CP10K) normalization
#' @param vargenes_method Variable gene selection method (either 'vst' or 'mvp')
#' @param topn Number of variable genes to subset by
#' @param tau Tau parameter for Harmony step
#' @param theta Theta parameter(s) for Harmony step
#' @param save_uwot_path Absolute path to save the uwot model (if do_umap is TRUE)
#' @param d Number of PC dimensions
#' @param additional_genes Any custom genes (e.g. marker genes) to include in addition to variable genes
#' 
#' @import data.table
#' @import tibble
#' @import irlba
#' @importFrom rlang .data
#'
#' @export
buildReference <- function(exp_ref,                   # Genes x cells
                           metadata_ref, 
                           vars = NULL,   
                           K = 50,                    # Number of soft clusters for Harmony
                           verbose = FALSE, 
                           do_umap = TRUE,
                           do_normalize = TRUE, 
                           vargenes_method = 'vst',   # vst or mvp
                           topn = 2000, 
                           tau = 0, 
                           theta = 2,
                           save_uwot_path = NULL,     # Path to save uwot model (use absolute path)
                           d = 20,
                           additional_genes = NULL) {
    
    set.seed(111) # for reproducible soft k-means and UMAP
    
    res = list(meta_data = metadata_ref)
        
    if (do_normalize) {
        if (verbose) message('Normalizing')
        exp_ref = normalizeData(exp_ref, 1e4, 'log')
    } 
    
    if (verbose) message('Finding variable genes')
    if (vargenes_method == 'mvp') {
        vargenes_df = findVariableGenes(exp_ref, rep('A', ncol(exp_ref)), num.bin = 20)
        var_genes = unique(data.table(vargenes_df)[, head(.SD[order(-.data$gene_dispersion_scaled)], topn), 
                                                        by = .data$group][, .data$symbol])
    } else if (vargenes_method == 'vst') {
        var_genes = vargenes_vst(exp_ref, topn = topn)
    } else {
        message("Invalid variable gene selection method. Options are 'vst' or 'mvp'.")
    }
    
    if(!is.null(additional_genes)) { # Add any additional genes
        var_genes = union(var_genes, additional_genes)
    }

    exp_ref = exp_ref[var_genes, ] # Subset gene expression matrix by the desired genes
        
    if (verbose) message('Scaling and PCA')
    vargenes_means_sds = tibble(symbol = var_genes, mean = Matrix::rowMeans(exp_ref))
    vargenes_means_sds$stddev = rowSDs(exp_ref, vargenes_means_sds$mean)
        
    # Scale data
    exp_ref_scaled = scaleDataWithStats(exp_ref, vargenes_means_sds$mean, vargenes_means_sds$stddev, 1)
    
    # PCA
    s = irlba::irlba(exp_ref_scaled, nv = d)
    Z_pca_ref = diag(s$d) %*% t(s$v) # [PCs by cells]
    res$loadings = s$u
    res$vargenes = vargenes_means_sds
    
    if (!is.null(vars)) {
        if (verbose) message('start Harmony')
        
        # Run Harmony to harmonize the reference
        ref_harmObj = harmony::HarmonyMatrix(
            data_mat = t(Z_pca_ref),  ## PCA embedding matrix of cells
            meta_data = metadata_ref, ## dataframe with cell labels
            theta = theta,            ## cluster diversity enforcement
            tau = tau,
            vars_use = vars,          ## variable to integrate over
            nclust = K,               ## number of clusters in Harmony model
            max.iter.harmony = 20,
            return_object = TRUE,     ## return the full Harmony model object
            do_pca = FALSE            ## do not recompute PCs
        )

        res$centroids <- t(cosine_normalize_cpp(ref_harmObj$R %*% t(ref_harmObj$Z_corr) , 1))
        res$R <- ref_harmObj$R
        res$betas <- harmony::moe_ridge_get_betas(ref_harmObj)
        res$Z_orig <- Z_pca_ref
        res$Z_corr <- ref_harmObj$Z_corr
        res$K <- K
        res$d <- d
    } else {
        clust_res <- soft_kmeans(Z_pca_ref, K)
        res$centroids <- clust_res$Y
        res$R <- clust_res$R
        res$betas <- NULL
        res$obj <- NULL
        res$Z_orig <- Z_pca_ref
        res$Z_corr <- Z_pca_ref
    }
    
    # Compute reference compression terms
    res$cache = compute_ref_cache(res$R, res$Z_corr)

    if (do_umap) {
        if (verbose) message('UMAP')
        res$umap <- uwot::umap(
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

    return(res)
}