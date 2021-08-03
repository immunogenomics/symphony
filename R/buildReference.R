#' Function for building a Symphony reference starting from expression matrix
#'
#' @param exp_ref Reference gene expression (genes by cells)
#' @param metadata_ref Reference cell metadata (cells by attributes)
#' @param vars Reference variables to Harmonize over e.g. c('donor', 'technology')
#' @param K Number of soft cluster centroids in model
#' @param verbose Verbose output
#' @param do_umap Perform UMAP visualization on harmonized reference embedding
#' @param do_normalize Perform log(CP10K+1) normalization
#' @param vargenes_method Variable gene selection method (either 'vst' or 'mvp')
#' @param vargenes_groups Name of metadata column specifying groups for variable gene selection. If not NULL, calculate topn variable genes in each group separately, then pool
#' @param topn Number of variable genes to subset by
#' @param tau Tau parameter for Harmony step
#' @param theta Theta parameter(s) for Harmony step
#' @param save_uwot_path Absolute path to save the uwot model (used if do_umap is TRUE)
#' @param d Number of PC dimensions
#' @param additional_genes Any custom genes (e.g. marker genes) to include in addition to variable genes
#' @param umap_min_dist umap parameter (see uwot documentation for details)
#' @param seed Random seed
#' 
#' @import data.table
#' @import tibble
#' @import irlba
#' @importFrom rlang .data
#' @return Symphony reference object. Integrated embedding is stored in the $Z_corr slot. Other slots include
#' cell-level metadata ($meta_data), variable genes means and standard deviations ($vargenes),
#' loadings from PCA ($loadings), original PCA embedding ($Z_orig), reference compression terms ($cache), 
#' betas from Harmony integration ($betas), cosine normalized soft cluster centroids ($centroids), 
#' centroids in PC space ($centroids_pc), and optional umap coordinates ($umap$embedding).
#' 
#' @export
buildReference <- function(exp_ref,                   # genes x cells
                           metadata_ref,              # cells x metadata fields
                           vars = NULL,               # metadata variables to Harmonize over
                           K = 100,                    # number of soft clusters for Harmony
                           verbose = FALSE,           # verbose output
                           do_umap = TRUE,            # run umap on reference cells
                           do_normalize = TRUE,       # run log(CP10K+1) normalization
                           vargenes_method = 'vst',   # vst or mvp
                           vargenes_groups = NULL,    # metadata column specifying groups for vargene selection
                           topn = 2000,               # number of variable genes (per group)
                           tau = 0,                   # Harmony parameter
                           theta = 2,                 # Harmony parameter
                           save_uwot_path = NULL,     # Path to save uwot model (use absolute path)
                           d = 20,                    # number of dimensions for PCs
                           additional_genes = NULL,   # vector of any additional genes beyond vargenes to include
                           umap_min_dist = 0.1,       # umap parameter
                           seed = 111) {
    
    set.seed(seed) # for reproducible soft k-means and UMAP
    
    res = list(meta_data = metadata_ref)

    if (do_normalize) {
        if (verbose) message('Normalizing')
        exp_ref = normalizeData(exp_ref, 1e4, 'log')
    }
    
    if (verbose) message('Finding variable genes using ', vargenes_method, ' method')
    if (vargenes_method == 'mvp') {
        if (is.null(vargenes_groups)) {
            vargenes_df = findVariableGenes(exp_ref, rep('A', ncol(exp_ref)), num.bin = 20)
        } else { # groups specified
            vargenes_df = findVariableGenes(exp_ref, groups = as.character(metadata_ref[[vargenes_groups]]), 
                                            num.bin = 20)
        }
        var_genes = unique(data.table(vargenes_df)[, head(.SD[order(-.data$gene_dispersion_scaled)], topn), by = .data$group][, .data$symbol])
    } else if (vargenes_method == 'vst') {
        if (is.null(vargenes_groups)) {
            var_genes = vargenes_vst(exp_ref, topn = topn)
        } else { # groups specified
            var_genes = vargenes_vst(exp_ref, groups = as.character(metadata_ref[[vargenes_groups]]), topn = topn)
        }
    } else {
        message("Invalid variable gene selection method. Options are 'vst' or 'mvp'.")
    }
    
    # Add in any additional genes
    if(!is.null(additional_genes)) { 
        if (verbose) message('Adding ', length(additional_genes), ' additional genes')
        var_genes = union(var_genes, additional_genes)   
    }
    if (verbose) message('Total ' , length(var_genes), ' genes for downstream steps')

    # Subset gene expression matrix by the desired genes
    exp_ref = exp_ref[var_genes, ]
        
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
    
    # Run Harmony integration
    if (!is.null(vars)) {
        if (verbose) message('Running Harmony integration')
        
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
            do_pca = FALSE,           ## do not recompute PCs
            verbose = verbose
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
        res$Z_orig <- Z_pca_ref
        res$Z_corr <- Z_pca_ref
    }
    
    # Add row and column names
    colnames(res$Z_orig) = row.names(metadata_ref)
    rownames(res$Z_orig) = paste0("PC_", seq_len(nrow(res$Z_corr)))
    colnames(res$Z_corr) = row.names(metadata_ref)
    rownames(res$Z_corr) = paste0("harmony_", seq_len(nrow(res$Z_corr)))
    
    # Compute reference compression terms
    if (verbose) message('Computing reference compression terms')
    res$cache = compute_ref_cache(res$R, res$Z_corr)
    
    # Compute centroids in harmony PC space
    cluster_sizes = res$cache[[1]] %>% as.matrix()
    centroid_sums = t(res$Z_corr %*% t(res$R)) %>% as.data.frame()
    centroids_pc = sweep(centroid_sums, 1, cluster_sizes, "/")
    colnames(centroids_pc) = paste0("harmony_", seq_len(nrow(res$Z_corr)))
    rownames(centroids_pc) = paste0("centroid_", seq_len(nrow(res$R)))
    res$centroids_pc = centroids_pc

    if (do_umap) {
        if (verbose) message('Running UMAP')
        umap <- uwot::umap(
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
    return(res)
}
