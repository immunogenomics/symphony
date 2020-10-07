## Function for mapping query cells to reference
mapQuery <- function(exp_query, 
                     metadata_query, 
                     ref_obj,          # from Symphony reference building
                     vars = NULL,      # query batch variables to harmonize over
                     verbose = TRUE,
                     do_normalize = TRUE,
                     do_umap = TRUE) { 
    
    if (do_normalize) {
        if (verbose) message('Normalizing')
        exp_query <- singlecellmethods::normalizeData(exp_query, 1e4, 'log')
    }
    
    ## Synchronize and scale query genes
    if (verbose) message('Scaling query gene expression')
    
    # Find shared genes between reference and query
    idx_shared_genes = which(ref_obj$vargenes$symbol %in% rownames(exp_query))
    shared_genes = reference$vargenes$symbol[idx_shared_genes]
    
    # Subset and scale the query cells by reference means and standard deviations
    exp_query_scaled = singlecellmethods::scaleDataWithStats(
        exp_query[shared_genes, ], ref_obj$vargenes$mean[idx_shared_genes], ref_obj$vargenes$stddev[idx_shared_genes], 1)
    
    # To add rows of zeros for missing genes, start with full matrix of zeroes
    exp_query_scaled_sync = matrix(0, nrow = length(ref_obj$vargenes$symbol), ncol = ncol(exp_query))  

    # Rows getting filled with exp_query_scaled values, leaving rows of 0s where appropriate
    exp_query_scaled_sync[idx_shared_genes,] = exp_query_scaled
    rownames(exp_query_scaled_sync) = ref_obj$vargenes$symbol
    colnames(exp_query_scaled_sync) = colnames(exp_query)
    
    if (verbose) message('Project')
    ### 1. Project into PCs using ref loadings
    Z_pca_query = t(ref_obj$loadings) %*% exp_query_scaled_sync
    
    if (verbose) message('Cluster')
    ### 2. Soft cluster assignment
    Z_pca_query_cos <- singlecellmethods:::cosine_normalize_cpp(Z_pca_query, 2)
    R_query = soft_cluster(ref_obj$centroids, Z_pca_query_cos, 0.1)
    
    if (verbose) message('Correct')
    ### 3. Correct
    
    # Make query design matrix
    if (!is.null(vars)) {
        design = droplevels(metadata_query)[,vars] %>% as.data.frame()
        
        onehot = design %>% 
            map(function(.x) {
                if (length(unique(.x)) == 1) { # special case if factor only has 1 level
                    rep(1, length(.x))
                } else {
                    model.matrix(~0 + .x)
                }
            }) %>% reduce(cbind)
        
        Xq = cbind(1, intercept = onehot) %>% t()
    } else { 
        #Treat query as 1 batch
        Xq = Matrix(rbind(rep(1, ncol(Z_pca_query)), rep(1, ncol(Z_pca_query))), sparse = TRUE)
    }
    
    # Mixture of experts correction (calls cpp code)
    Zq_corr = moe_correct_ref(as.matrix(Z_pca_query), 
                              as.matrix(Xq), 
                              as.matrix(R_query), 
                              as.matrix(ref_obj$cache[[1]]), 
                              as.matrix(ref_obj$cache[[2]]))
    
    if (verbose) message('UMAP')
    ## UMAP projection of query if the reference umap is present
    umap_query = NULL
    if (do_umap & !is.null(ref_obj$save_uwot_path)) {
        ref_umap_model = uwot::load_uwot(ref_obj$save_uwot_path, verbose = FALSE)
        umap_query = uwot::umap_transform(t(Zq_corr), ref_umap_model)
    }
    
    if (verbose) message('All done!')
    return(list(Z = Zq_corr, Zq_pca = Z_pca_query, R = R_query, Xq = Xq, 
                umap = umap_query, meta_data = metadata_query))
}

#Note: this will need to be fixed if we're including merging redundant clusters in ref building
probPredict = function(query_obj, ref_obj) {
    ## Predict cell type using probabilistic method
    type_probs <- crossprod(query_obj$R, ref_obj$cluster_annotations)
    cell_type <- colnames(type_probs)[apply(type_probs, 1, which.max)]
    type_thresh <- 0.2
    cell_type[which(apply(type_probs, 1, max) < type_thresh)] <- 'unassigned'
    query_obj$meta_data$cell_type_pred_prob <- as.factor(cell_type)
    return(query_obj)
}

## Predict cell annotations using knn method
knnPredict <- function(query_obj, ref_obj,
                       train_labels, # cell labels for knn classification 
                       k = 5) {
    knn_pred = class::knn(t(ref_obj$Z_corr), t(query_obj$Z), train_labels, k = k)
    query_obj$meta_data$cell_type_pred_knn = knn_pred
    return(query_obj)
}

## Predict cell type using knn method with cos distance
knnPredictCos <- function(query_obj, ref_obj,
                       train_labels, # cell labels for knn classification 
                       k = 5) {
    Z_ref_cos = singlecellmethods:::cosine_normalize_cpp(ref_obj$Z_corr, 2)
    Z_query_cos = singlecellmethods:::cosine_normalize_cpp(query_obj$Z, 2)
    knn_pred = class::knn(t(Z_ref_cos), t(Z_query_cos), train_labels, k = k)
    query_obj$meta_data$cell_type_pred_knn_cos = knn_pred
    return(query_obj)
}

# Function for building a reference.
buildReference <- function(exp_ref,       #[genes x cells]
                           metadata_ref, 
                           vars = NULL,   # variables to Harmonize over e.g. c('donor', 'technology')
                           #cell_types,    # reference annotations (removed 9/27/20)
                           K = 50,        # number of soft clusters for Harmony
                           verbose = FALSE, 
                           do_umap = TRUE, # run UMAP on the reference cells?
                           weightedPCA = FALSE,
                           do_normalize = TRUE, 
                           pca_function = 'irlba', # use 'irlba' rather than 'rsvd'
                           vargenes_method = 'vst', # vst or mvp
                           topn = 2000, # number of variable genes to subset by
                           tau = 0,     # tau for Harmony
                           theta = 2,   # theta for Harmony
                           save_uwot_path = NULL, # path to save uwot model (use absolute path)
                           d = 20,                # number of dimensions
                           additional_genes = NULL) { # custom genes (e.g. markers) to include beyond variable genes
    
    set.seed(111) # for reproducible soft k-means and UMAP
    
    res <- list(meta_data = metadata_ref)
        
    if (do_normalize) {
        if (verbose) message('start normalizing')
        exp_ref <- singlecellmethods::normalizeData(exp_ref, 1e4, 'log')
    } 
    
    if (verbose) message('start finding variable genes')
    if (vargenes_method == 'mvp') {
        vargenes_df <- singlecellmethods::findVariableGenes(exp_ref, rep('A', ncol(exp_ref)), num.bin = 20)
        var_genes <- unique(data.table(vargenes_df)[, head(.SD[order(-gene_dispersion_scaled)], topn), 
                                                        by = group][, symbol])
    } else if (vargenes_method == 'vst') {
        var_genes = vargenes_vst(exp_ref, topn = topn)
    } else {
        message('Invalid variable gene selection method. Options are vst or mvp.')
    }
    
    if(!is.null(additional_genes)) { # Add additional genes
        var_genes = union(var_genes, additional_genes)
    }

    exp_ref <- exp_ref[var_genes, ] # Subset gene expression matrix by the desired genes
    
    if (weightedPCA) {
        if (verbose) message('start weighted scaling and weighted PCA')
        s = singlecellmethods::weighted_pca(as(exp_ref, 'dgCMatrix'), metadata_ref$cell_type)
        res$vargenes = s$vargenes # weighted_pca returns weighted means and sds for scaling query
        res$loadings = s$loadings
        Z_pca_ref = t(s$embeddings)
        
    } else { #regular PCA
        
        if (verbose) message('start scaling and PCA')
        vargenes_means_sds <- tibble(symbol = var_genes, mean = Matrix::rowMeans(exp_ref))
        vargenes_means_sds$stddev <- singlecellmethods::rowSDs(exp_ref, vargenes_means_sds$mean)
        
        # Scale data
        exp_ref_scaled <- singlecellmethods::scaleDataWithStats(exp_ref, vargenes_means_sds$mean,
                                                    vargenes_means_sds$stddev, 1)
        # PCA
        if (pca_function == 'rsvd') {
            s <- rsvd::rsvd(exp_ref_scaled, k = d)
        } else if (pca_function == 'svd') {
            s = svd(exp_ref_scaled)
        } else if (pca_function == 'irlba') {
            s = irlba(exp_ref_scaled, nv = d)
        } else {
            message('Invalid PCA method. Options are rsvd, svd, or irlba.')
        }
        
        Z_pca_ref = diag(s$d) %*% t(s$v) # [pcs by cells]
        res$loadings <- s$u
        res$vargenes <- vargenes_means_sds
    }
    
    if (!is.null(vars)) {
        if (verbose) message('start Harmony')
        
        # Run Harmony to harmonize the reference
        ref_harmObj = harmony::HarmonyMatrix(
            data_mat = t(Z_pca_ref), ## PCA embedding matrix of cells
            meta_data = metadata_ref, ## dataframe with cell labels
            theta = theta, ## cluster diversity enforcement
#           sigma = 0.1,
            tau = tau,
            vars_use = vars, ## variable to integrate out
            nclust = K, ## number of clusters in Harmony model
            max.iter.harmony = 20,
            return_object = TRUE, ## return the full Harmony model object
            do_pca = FALSE ## don't recompute PCs
        )

        res$centroids <- t(singlecellmethods:::cosine_normalize_cpp(ref_harmObj$R %*% t(ref_harmObj$Z_corr) , 1))
        res$R <- ref_harmObj$R
        res$betas <- harmony:::moe_ridge_get_betas(ref_harmObj)
        #res$obj <- ref_harmObj # pointer ends up being undefined
        res$Z_orig <- Z_pca_ref
        res$Z_corr <- ref_harmObj$Z_corr
        res$K <- K
        res$d <- d
    } else {
        clust_res <- singlecellmethods::soft_kmeans(Z_pca_ref, K)
        res$centroids <- clust_res$Y
        res$R <- clust_res$R
        res$betas <- NULL
        res$obj <- NULL
        res$Z_orig <- Z_pca_ref
        res$Z_corr <- Z_pca_ref
    }
    
    # Compute reference compression terms (calls cpp code)
    res$cache <- compute_ref_cache(res$R, res$Z_corr)

    if (do_umap) {
        if (verbose) message('start UMAP')
        res$umap <- uwot::umap(
            t(res$Z_corr), n_neighbors = 30, learning_rate = 0.5, init = "laplacian", 
            metric = 'cosine', fast_sgd = FALSE, n_sgd_threads = 1, # for reproducibility
            min_dist = .1, n_threads = 4, ret_model = TRUE
        )
        
        # Since the nn-index component of the uwot model is not able to be saved as an 
        # RDS object, we save the uwot model entirely at a user-defined path.
        if (!is.null(save_uwot_path)) {
            model = uwot::save_uwot(res$umap, file = save_uwot_path, unload = FALSE, verbose = FALSE)
            res$save_uwot_path = save_uwot_path
            if (verbose) message(paste('saved uwot model'))
        }
    }
    
    ## Save cell type labels 
    #if (verbose) message('annotate clusters')
    #cell_types <- factor(cell_types)
    #cluster_annotations <- res$R %*% model.matrix(~0 + cell_types)
    #cluster_annotations <- diag(1 / rowSums(cluster_annotations)) %*% cluster_annotations
    #colnames(cluster_annotations) <- gsub('cell_types', '', colnames(cluster_annotations))
    #res$cluster_annotations <- cluster_annotations

    return(res)
}
     

buildReferenceFromHarmonyObj <- function(harmony_obj, #output object from HarmonyMatrix()
                           metadata,
                           vargenes_means_sds,     # gene names, means, and std devs for scaling
                           pca_loadings,           # genes x PCs
                           verbose = TRUE, 
                           do_umap = TRUE,         # run reference umap?
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
        # RDS object, we save the uwot model entirely at a user-defined path.
        if (!is.null(save_uwot_path)) {
            model = uwot::save_uwot(res$umap, file = save_uwot_path, unload = FALSE, verbose = FALSE)
            res$save_uwot_path = save_uwot_path
            if (verbose) message(paste('saved uwot model'))
        }
    }
    return(res)
}



# Function for evaluating F1 by cell type, modified from benchmarking paper Abdelaal et al. 2019
evaluate <- function(true, predicted){
  "
  Parameters
  ----------
  TrueLabelsPath: csv file with the true labels (format: one column, no index)
  PredLabelsPath: csv file with the predicted labels (format: one column, no index)
  Indices: which part of the csv file should be read (e.g. if more datasets are tested at the same time) (format: c(begin, end))
  
  Returns
  -------
  Conf: confusion matrix
  MedF1 : median F1-score
  F1 : F1-score per class
  Acc : accuracy
  PercUnl : percentage of unlabeled cells
  PopSize : number of cells per cell type
  "
  
  true_lab <- unlist(true)
  pred_lab <- unlist(predicted)
  
  unique_true <- unlist(unique(true_lab))
  unique_pred <- unlist(unique(pred_lab))
  
  unique_all <- unique(c(unique_true,unique_pred))
  conf <- table(true_lab,pred_lab)
  pop_size <- rowSums(conf)
  
  #pred_lab = gsub('Node..','Node',pred_lab)
  
  conf_F1 <- table(true_lab,pred_lab,exclude = c('unassigned','Unassigned','Unknown','rand','Node','ambiguous','unknown'))

  F1 <- vector()
  sum_acc <- 0
  
  for (i in c(1:length(row.names(conf_F1)))){
    if(pop_size[row.names(conf_F1)[i]] == 0) {
        F1[i] = NA # F1 score is N/A
        next
    }
      
    findLabel = colnames(conf_F1) == row.names(conf_F1)[i]
    
    if(sum(findLabel)){
      prec <- conf_F1[i,findLabel] / colSums(conf_F1)[findLabel]
      rec <- conf_F1[i,findLabel] / rowSums(conf_F1)[i]
      if (prec == 0 | rec == 0){
        F1[i] = 0
      } else{
        F1[i] <- (2*prec*rec) / (prec + rec)
      }
      sum_acc <- sum_acc + conf_F1[i,findLabel]
    } else {
      F1[i] = 0
    }
    #print(paste(row.names(conf_F1)[i], ' --> ', sum(findLabel), 'F1 =', F1[i]))
  }
  
  #pop_size <- pop_size[pop_size > 0]
  
  #names(F1) <- names(pop_size)
  names(F1) = row.names(conf_F1)
  
  med_F1 <- median(na.omit(F1))
  
  total <- length(pred_lab)
  num_unlab <- sum(pred_lab == 'unassigned') + sum(pred_lab == 'Unassigned') + sum(pred_lab == 'rand') + sum(pred_lab == 'Unknown') + sum(pred_lab == 'unknown') + sum(pred_lab == 'Node') + sum(pred_lab == 'ambiguous')
  per_unlab <- num_unlab / total
  
  acc <- sum_acc/sum(conf_F1)
  
  result <- list(Conf = conf, MedF1 = med_F1, F1 = F1, Acc = acc, PercUnl = per_unlab, PopSize = pop_size)
  
  return(result)
}