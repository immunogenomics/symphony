#' Calculates the k-NN correlation, which measures how well the sorted ordering of k nearest reference
#' neighbors in a gold standard embedding correlate with the ordering for the same reference cells in 
#' an alternative embedding (i.e. from reference mapping).
#' NOTE: it is very important for the order of reference cells (cols) in gold_ref matches that of alt_ref
#' (same for matching columns of gold_query and alt_query).
#'
#' @param gold_ref Reference cells in gold standard embedding (PCs by cells)
#' @param alt_ref Reference cells in alternative embedding (PCs by cells)
#' @param gold_query Query cells in gold standard embedding (PCs by cells)
#' @param alt_query Query cells in alternative embedding (PCs by cells)
#' @param k Number of reference neighbors to use for kNN-correlation calculation
#' 
#' @import RANN
#' @return Vector of k-NN correlations for query cells
#' @export
calcknncorr = function(gold_ref, alt_ref, gold_query, alt_query, k = 500) {
    message('Note: This function assumes that ordering of cells (cols) between gold and alt embeddings match')
    # Calculate the query cells' k nearest reference neighbors in the gold standard embedding
    nn_in_gold = nn2(t(gold_ref), t(gold_query), k = k)

    corrs = numeric(ncol(gold_query)) # initialize results
    for (i in 1:nrow(nn_in_gold$nn.idx)) { # for each query cell
        neighbors_idx = nn_in_gold$nn.idx[i, ]
    
        # Get position of anchor cell in gold embedding
        query_anchor_gold = gold_query[, i] %>% # duplicate k times
            matrix(nrow = nrow(gold_query), ncol = k, byrow = FALSE)
    
        # Get position of anchor cell in alternate embedding
        query_anchor_alt = alt_query[, i] %>% # duplicate k times
            matrix(nrow = nrow(gold_query), ncol = k, byrow = FALSE)
    
        # Get positions for its nearest neighbors in the gold embedding
        query_neighbors_gold = gold_ref[, neighbors_idx]
    
        # Get positions for its nearest neighbors in the alt embedding
        query_neighbors_alt = alt_ref[, neighbors_idx]
    
        # Calculate distance between anchor cell and each neighbor in gold embedding
        distances_gold = sqrt(colSums((query_anchor_gold - query_neighbors_gold)**2))
    
        # Calculate distance between anchor cell and each neighbor in alt embedding
        distances_alt = sqrt(colSums((query_anchor_alt - query_neighbors_alt)**2))
    
        # Calculate Spearman correlation between the two distance vectors
        corrs[i] = cor(distances_gold, distances_alt, method = 'spearman')
    }
    return(corrs)
}

#' Calculates the k-NN correlation within the query cells only, which measures how well the sorted 
#' ordering of k nearest query neighbors in a query de novo PCA embedding correlate with the ordering 
#' for the cells in the reference mapping embedding.
#'
#' @param query Query object (returned from mapQuery)
#' @param var Query metadata batch variable (PCA is calculated within each batch separately); if NULL, do not split by batch
#' @param k Number of neighbors to use for kNN-correlation calculation
#' @param topn number of variable genes to calculate within each query batch for query PCA
#' @param d number of dimensions for query PCA within each query batch
#' @param distance either 'euclidean' or 'cosine'
#' 
#' @import RANN
#' @return Vector of within-query k-NN correlations for query cells
#' @export
calcknncorrWithinQuery = function(query, var = NULL, k = 100, topn = 2000, d = 20, distance = 'euclidean') {
    corrs = numeric(nrow(query$meta_data)) # initialize results
    
    if (!is.null(var)) {
        for (batch in unique(query$meta_data[[var]])) { # for each batch
            message(paste0('Calculating k-NN correlation within query batch ', batch))
        
            batch_idx = which(query$meta_data[[var]] == batch)
            query$exp = Matrix(query$exp, sparse = TRUE)
            query_exp_batch = query$exp[, batch_idx]
        
            Z_pca = runPCAQueryAlone(query_exp_batch, topn = topn, d = d)
            Z_mapping = query$Z[, batch_idx]
        
            # Calculate correlation & save results
            corrs[batch_idx] = calcknncorrWithinQueryBatch(Z_pca, Z_mapping, k, distance = distance)
        }
        
    } else {
        message('No batch var specified. Treating query as 1 batch.')
        Z_pca = runPCAQueryAlone(query$exp, topn = topn, d = d)
        Z_mapping = query$Z
        
        # Calculate correlation & save results
        corrs = calcknncorrWithinQueryBatch(Z_pca, Z_mapping, k, distance = distance)
    }
    return(corrs)
}

# Non-exported (called by calcknncorrWithinQuery above)
# Calculates the k-NN correlation within a query batch.
calcknncorrWithinQueryBatch = function(Z_pca, Z_mapping, k, distance) {
    if (!identical(dim(Z_pca), dim(Z_mapping))) {
        stop('Error: PCA and mapping embeddings have different dimensions')
    }
    
    if (distance == 'cosine') { # L2 normalize
        Z_pca = Z_pca %>% cosine_normalize_cpp(2)
        Z_mapping = Z_mapping %>% cosine_normalize_cpp(2)
    }
    
    # Calculate nearest neighbors in query PCA space
    nn_in_query_pca = nn2(t(Z_pca), t(Z_pca), k = k + 1) # k+1 because do not count itself
    rownames(nn_in_query_pca$nn.idx) = nn_in_query_pca$nn.idx[, 1]
    nn_in_query_pca$nn.idx = nn_in_query_pca$nn.idx[, -1] # do not count itself
    
    corrs_batch = numeric(nrow(nn_in_query_pca$nn.idx))
    
    for (i in 1:ncol(Z_pca)) { # For each anchor cell, calculate k-NN-corr
        if(ncol(Z_pca) < k + 1) {
            k = ncol(Z_pca) - 1
            message(paste('Warning: Batch has too few cells. Using k =', k , 'instead.'))
        }
        
        neighbors_idx = nn_in_query_pca$nn.idx[i, ] # neighbors are defined in query PC space within each batch
        
        # Get position of anchor cell in query PCA embedding
        query_anchor_pca = Z_pca[, i] %>% 
                                matrix(nrow = length(Z_pca[, i]), ncol = k, byrow = FALSE)
        
        # Get position of anchor cell in mapping embedding
        query_anchor_mapping = Z_mapping[, i] %>% matrix(nrow = length(Z_mapping[, i]), ncol = k, byrow = FALSE)
        
        # Get positions for its nearest neighbors in the query PCA embedding
        query_neighbors_pca = Z_pca[, neighbors_idx] # [20 x k]
    
        # Get positions for its nearest neighbors in the mapping embedding
        query_neighbors_mapping = Z_mapping[, neighbors_idx] # [20 x k]
    
        # Calculate Euclidean distance between anchor cell and each neighbor in query PCA embedding
        distances_pca = sqrt(colSums((query_anchor_pca - query_neighbors_pca)**2))
    
        # Calculate Euclidean distance between anchor cell and each neighbor in mapping embedding
        distances_mapping = sqrt(colSums((query_anchor_mapping - query_neighbors_mapping)**2))
    
        # Calculate Spearman correlation between the two distance vectors
        corrs_batch[i] = cor(distances_mapping, distances_pca, method = 'spearman')
    }
    return(corrs_batch)
}

#' Runs a standard PCA pipeline on query (1 batch). Assumes query_exp is already normalized.
#'
#' @param query_exp Query expression matrix (genes x cells)
#' @param topn Number of variable genes to use
#' @param d Number of dimensions
#' @param seed random seed
#'
#' @import irlba
#' @return A matrix of PCs by cells
#' @export
runPCAQueryAlone = function(query_exp, topn = 2000, d = 20, seed = 1) {
    # Subset by variable genes
    vargenes = vargenes_vst(query_exp, topn = topn)
    vargenes_exp = query_exp[vargenes, ]

    vargenes_means_sds = tibble(symbol = vargenes, mean = Matrix::rowMeans(vargenes_exp))
    vargenes_means_sds$stddev <- rowSDs(vargenes_exp, vargenes_means_sds$mean)
        
    # Scale data
    exp_scaled <- scaleDataWithStats(vargenes_exp, vargenes_means_sds$mean, vargenes_means_sds$stddev, 1)
    
    # Run SVD, save loadings
    set.seed(seed)
    s = irlba(exp_scaled, nv = d)
    Z_pca = diag(s$d) %*% t(s$v) # [PCs by cells]
    return(Z_pca)
}
