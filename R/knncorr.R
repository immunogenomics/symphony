#' Calculates the k-NN correlation, which measures how well the sorted ordering of k nearest reference
#' neighbors in a gold standard embedding correlate with the ordering for the same reference cells in 
#' an alternative embedding (i.e. from reference mapping).
#'
#' @param gold_ref Reference cells in gold standard embedding [pcs x cells]
#' @param alt_ref Reference cells in alternative embedding [pcs x cells]
#' @param gold_query Query cells in gold standard embedding [pcs x cells]
#' @param alt_query Query cells in alternative embedding [pcs x cells]
#' @param k Number of reference neighbors to use for kNN-correlation calculation
#' 
#' @import RANN
#' 
#' @export
calcknncorr = function(gold_ref, alt_ref, gold_query, alt_query, k = 500) {
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