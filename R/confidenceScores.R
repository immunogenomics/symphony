#' ---- Per-cell Confidence Score ----
#' Calculates the weighted Mahalanobis distance for the query cells to reference clusters. Returns a vector
#' of distance scores, one per query cell. Higher distance indicates less confidence.
#'
#' @param reference Reference object as returned by Symphony buildReference()
#' @param query Query object as returned by Symphony mapQuery()
#' @param Z_orig Define reference distribution using original PCA embedding or harmonized PC embedding
#' @param metric Uses Mahalanobis by default, but added as a parameter for potential future use
#' 
#' @import utils
#' @import stats
#' @export
calcPerCellConfidence = function(reference, query, Z_orig = FALSE, metric = 'mahalanobis') {
    
    ### Calculate (weighted) covariance matrix and centroid for all k soft clusters
    
    # initialize a d x k matrix
    center_ks = matrix(rep(0, len = ncol(reference$centroids) * nrow(reference$Z_corr)), nrow = nrow(reference$Z_corr))
    # initialize k * (d * d) tensor
    cov_ks = list()
    
    # Calculate 
    for (k in 1:ncol(reference$centroids)) {
        if (Z_orig) {
            cov_k = cov.wt(t(reference$Z_orig), wt = reference$R[k,])
        } else {
            cov_k = cov.wt(t(reference$Z_corr), wt = reference$R[k,])
        }
        cov_ks[[k]] = cov_k$cov # covariance matrix
        center_ks[, k] = cov_k$center # centroid in hPC space (d x 1)
    }
    
    # Calculate the Mahalanobis distance from each query cell to all centroids
    mah_dist_ks = matrix(rep(0, len = ncol(query$Zq_pca) * ncol(reference$centroids)), nrow = ncol(query$exp))
    for (k in 1:ncol(reference$centroids)) {
        mah_dist_ks[, k] = sqrt(mahalanobis(x = t(query$Zq_pca), center = center_ks[, k], cov = cov_ks[[k]]))
    }
    
    # Return the per-cell score, which is the average of the distances weighted by the clusters the cell belongs to
    maha = rowSums(mah_dist_ks * t(query$R))
    return(maha)
}

#' ---- Per-cluster Confidence Score ----
#' Calculates the Mahalanobis distance from user-defined query clusters to their nearest
#' reference centroid after initial projection into reference PCA space. 
#' All query cells in a cluster get the same score. Higher distance indicates less confidence.
#'
#' @param reference Reference object as returned by Symphony buildReference()
#' @param query Query object as returned by Symphony mapQuery()
#' @param query_cluster_labels Vector of user-defined labels denoting clusters / putative novel cell type to calculate the score for
#' @param metric Uses Mahalanobis by default, but added as a parameter for potential future use
#' 
#' @import utils 
#' @import stats
#' @export
calcPerClusterConfidence = function(reference, query, query_cluster_labels, metric = 'mahalanobis') {
    
    query_cluster_labels = as.character(query_cluster_labels)
    query_cluster_labels_unique = unique(query_cluster_labels)
    num_clusters = length(query_cluster_labels_unique)
    message('Calculating mapping confidence for ', num_clusters, ' query clusters')
    
    ### Calculate the Mahalobinis distance between each query cluster and it's nearest reference centroid
    # c denotes the number of user-defined query clusters
    
    # initialize a d x k matrix
    center_cs = matrix(rep(0, len = nrow(query$Zq_pca) * num_clusters), nrow(query$Zq_pca)) # init
    # initialize c * (d * d) tensor
    cov_cs = list()
    colnames(center_cs) = query_cluster_labels_unique
    
    # Calculate query cluster centroid and covariances in PC space
    for (c in 1:num_clusters) {
        cluster_idx = which(query_cluster_labels == query_cluster_labels_unique[c])
        cluster_Zq_pca = query$Zq_pca[, cluster_idx]
        cov_cs[[c]] = cov(t(cluster_Zq_pca))
        center_cs[, c] = rowMeans(cluster_Zq_pca)
    }
    
    ## Find nearest reference cluster centroid
    nearest_centroid_idx = max.col(t(center_cs) %*% t(reference$centroids_pc))
    centroid_closest = reference$centroids_pc[nearest_centroid_idx, ]

    # Calculate Mahalanobis distance from query cluster to nearest reference centroid
    mah_dist_cs = as.data.frame(matrix(rep(0, len = num_clusters * 2), nrow = num_clusters)) # init
    colnames(mah_dist_cs) = c('query_cluster', 'distance_score') # init
    mah_dist_cs$query_cluster = query_cluster_labels_unique
    
    for (c in 1:num_clusters) {
        # if the number of cells in a query cluster is less than d, we add a ridge to the diagonal
        cluster_size = length(which(query_cluster_labels == query_cluster_labels_unique[c]))
        if (cluster_size <  nrow(query$Z)) {
            message('(Warning) cluster contains too few cells to estimate confidence: ', query_cluster_labels_unique[c])
            mah_dist_cs$distance_score[c] = NA
        } else {
            cov = cov_cs[[c]] +  1 * diag(nrow(query$Z)) # ridge to help stabilize numerical estimates
            mah_dist_cs$distance_score[c] = mahalanobis(x = centroid_closest[c,], center = center_cs[,c], cov = cov)
        }
    }

    return(mah_dist_cs)
}