#' Function for mapping query cells to a Symphony reference
#'
#' @param exp_query Query gene expression (genes by cells)
#' @param metadata_query Query metadata (cells by attributes)
#' @param ref_obj Reference object as returned by Symphony buildReference()
#' @param vars Query batch variable(s) to integrate over (column names in metadata)
#' @param verbose Verbose output
#' @param do_normalize Perform log(CP10K) normalization on query expression
#' @param do_umap Perform umap projection into reference UMAP (if reference includes a uwot model)
#' @param sigma Fuzziness parameter for soft clustering (sigma = 1 is hard clustering)
#' 
#' @import utils
#' @importFrom magrittr "%>%"
#' @importFrom Matrix Matrix
#' 
#' @export
mapQuery = function(exp_query, 
                     metadata_query, 
                     ref_obj,          # From Symphony reference building
                     vars = NULL,      # Query batch variables to harmonize over
                     verbose = TRUE,
                     do_normalize = TRUE,
                     do_umap = TRUE,
                     sigma = 0.1) { 
    
    if (do_normalize) {
        if (verbose) message('Normalizing')
        exp_query = normalizeData(exp_query, 1e4, 'log')
    }
    
    ## Synchronize and scale query genes
    if (verbose) message('Scaling and synchronizing query gene expression')
    
    # Find shared genes between reference and query
    idx_shared_genes = which(ref_obj$vargenes$symbol %in% rownames(exp_query))
    shared_genes = ref_obj$vargenes$symbol[idx_shared_genes]
    
    # Subset and scale the query cells by reference means and standard deviations
    exp_query_scaled = scaleDataWithStats(exp_query[shared_genes, ],
                                          ref_obj$vargenes$mean[idx_shared_genes],
                                          ref_obj$vargenes$stddev[idx_shared_genes], 1)
    
    # To add rows of zeros for missing genes, start with full matrix of zeroes
    exp_query_scaled_sync = matrix(0, nrow = length(ref_obj$vargenes$symbol), ncol = ncol(exp_query))  

    # Rows get filled with exp_query_scaled values, leaving rows of 0s where appropriate
    exp_query_scaled_sync[idx_shared_genes, ] = exp_query_scaled
    rownames(exp_query_scaled_sync) = ref_obj$vargenes$symbol
    colnames(exp_query_scaled_sync) = colnames(exp_query)
    
    if (verbose) message('Project') 
    ### 1. Project into PCs using reference loadings
    Z_pca_query = t(ref_obj$loadings) %*% exp_query_scaled_sync
    
    if (verbose) message('Cluster')
    ### 2. Soft cluster assignment
    Z_pca_query_cos = cosine_normalize_cpp(Z_pca_query, 2)
    R_query = soft_cluster(ref_obj$centroids, Z_pca_query_cos, sigma)
    
    if (verbose) message('Correct')
    ### 3. Correction
    
    # Make query design matrix
    if (!is.null(vars)) {
        design = droplevels(metadata_query)[,vars] %>% as.data.frame()
        
        onehot = design %>% 
            purrr::map(function(.x) {
                if (length(unique(.x)) == 1) { # Special case if factor only has 1 level
                    rep(1, length(.x))
                } else {
                    stats::model.matrix(~0 + .x)
                }
            }) %>% purrr::reduce(cbind)
        
        Xq = cbind(1, intercept = onehot) %>% t()
    } else { 
        # If no batches specified, treat all query cells as 1 batch
        Xq = Matrix(rbind(rep(1, ncol(Z_pca_query)), rep(1, ncol(Z_pca_query))), sparse = TRUE)
    }
    
    # Mixture of experts correction (calls cpp code)
    Zq_corr = moe_correct_ref(as.matrix(Z_pca_query), 
                              as.matrix(Xq), 
                              as.matrix(R_query), 
                              as.matrix(ref_obj$cache[[1]]), 
                              as.matrix(ref_obj$cache[[2]]))
    
    if (verbose) message('UMAP')
    ## UMAP projection of query if the reference uwot model is present
    umap_query = NULL
    
    if (do_umap & !is.null(ref_obj$save_uwot_path)) {
        ref_umap_model = uwot::load_uwot(ref_obj$save_uwot_path, verbose = FALSE)
        umap_query = uwot::umap_transform(t(Zq_corr), ref_umap_model)
    }
    
    if (verbose) message('All done!')
    return(list(Z = Zq_corr, Zq_pca = Z_pca_query, R = R_query, Xq = Xq, 
                umap = umap_query, meta_data = metadata_query))
}