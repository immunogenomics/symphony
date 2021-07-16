buildReferenceFromSeurat <- function(
    obj, assay = 'RNA', verbose = TRUE, save_umap = TRUE, save_uwot_path = NULL
) {
    if(!assay %in% c('RNA', 'SCT')) {
        stop('Only supported assays are RNA or SCT.')
    }
    res <- list()
    ## TODO: check that these objects are all correctly initialized
    res$Z_corr <- t(obj@reductions$harmony@cell.embeddings)
    res$Z_orig <- t(obj@reductions$pca@cell.embeddings)
    message('Saved embeddings')
    
    res$R <- t(obj@reductions$harmony@misc$R)
    message('Saved soft cluster assignments')
    
    if (assay == 'RNA') {
        vargenes_means_sds <- tibble(
            symbol = obj@assays[[assay]]@var.features, 
            mean = Matrix::rowMeans(obj@assays[[assay]]@data[obj@assays[[assay]]@var.features, ])
        )
        vargenes_means_sds$stddev <- rowSDs(
            obj@assays[[assay]]@data[obj@assays[[assay]]@var.features, ], 
            vargenes_means_sds$mean
        )
    } else if (assay == 'SCT') {
        vargenes_means_sds <- tibble(
            symbol = obj@assays[[assay]]@var.features, 
            mean = Matrix::rowMeans(obj@assays[[assay]]@scale.data[obj@assays[[assay]]@var.features, ])
        )
        asdgc = Matrix(obj@assays[[assay]]@scale.data[obj@assays[[assay]]@var.features, ], sparse = TRUE)
        vargenes_means_sds$stddev <- rowSDs(
            asdgc, 
            vargenes_means_sds$mean
        )
    }
    
    res$vargenes_means_sds <- vargenes_means_sds
    message('Saved variable gene information for ', nrow(vargenes_means_sds), ' genes.')
    
    res$loadings <- obj@reductions$pca@feature.loadings
    message('Saved PCA loadings.')
    
    res$meta_data <- obj@meta.data
    message('Saved metadata.')
    
    ## Check UMAP 
    if (save_umap) {
        if (is.null(save_uwot_path)) {
            error('Please provide a valid path to save_uwot_path in order to save uwot model.')
        }
        if (is.null(obj@reductions$umap@misc$model)) {
            error('uwot model not initialiazed in Seurat object. Please do RunUMAP with umap.method=\'uwot\', return.model=TRUE first.')
        }
        res$umap <- obj@reductions$umap@misc$model
        res$save_uwot_path <- save_uwot_path
        if (file.exists(res$save_uwot_path)) {
            file.remove(res$save_uwot_path)    
        }
        uwot::save_uwot(res$umap, save_uwot_path)
    }
    
    ## Build Reference! 
    if (verbose) 
        message("Calculate final L2 normalized reference centroids (Y_cos)")
    res$centroids = t(cosine_normalize_cpp(res$R %*% t(res$Z_corr), 1))
    if (verbose) 
        message("Calculate reference compression terms (Nr and C)")
    res$cache = compute_ref_cache(res$R, res$Z_corr)
    colnames(res$Z_orig) = row.names(res$metadata)
    rownames(res$Z_orig) = paste0("PC_", seq_len(nrow(res$Z_corr)))
    colnames(res$Z_corr) = row.names(res$metadata)
    rownames(res$Z_corr) = paste0("harmony_", seq_len(nrow(res$Z_corr)))
        
    if (verbose) 
        message("Finished nicely.")
    return(res)    
}

environment(buildReferenceFromSeurat) <- environment(symphony::buildReference)

RunHarmony.Seurat <- function(
  object,
  group.by.vars,
  reduction = 'pca',
  dims.use = NULL,
  theta = NULL,
  lambda = NULL,
  sigma = 0.1,
  nclust = NULL,
  tau = 0,
  block.size = 0.05,
  max.iter.harmony = 10,
  max.iter.cluster = 20,
  epsilon.cluster = 1e-5,
  epsilon.harmony = 1e-4,
  plot_convergence = FALSE,
  verbose = TRUE,
  reference_values = NULL,
  reduction.save = "harmony",
  assay.use = 'RNA',
  project.dim = TRUE,
  ...
) {
  if (reduction == "pca") {
    tryCatch(
      embedding <- Seurat::Embeddings(object, reduction = "pca"),
      error = function(e) {
        if (verbose) {
          message("Harmony needs PCA. Trying to run PCA now.")
        }
        tryCatch(
          object <- Seurat::RunPCA(
            object,
            assay = assay.use, verbose = verbose
          ),
          error = function(e) {
            stop("Harmony needs PCA. Tried to run PCA and failed.")
          }
        )
      }
    )
  } else {
    available.dimreduc <- names(methods::slot(object = object, name = "reductions"))
    if (!(reduction %in% available.dimreduc)) {
      stop("Requested dimension reduction is not present in the Seurat object")
    }
    embedding <- Seurat::Embeddings(object, reduction = reduction)
  }
  if (is.null(dims.use)) {
    dims.use <- seq_len(ncol(embedding))
  }
  dims_avail <- seq_len(ncol(embedding))
  if (!all(dims.use %in% dims_avail)) {
    stop("trying to use more dimensions than computed. Rereun dimension reduction
         with more dimensions or run Harmony with fewer dimensions")
  }
  if (length(dims.use) == 1) {
    stop("only specified one dimension in dims.use")
  }
  metavars_df <- Seurat::FetchData(object, group.by.vars)
    
  harmonyObject <- HarmonyMatrix(
    embedding,
    metavars_df,
    group.by.vars,
    FALSE,
    0,
    theta,
    lambda,
    sigma,
    nclust,
    tau,
    block.size,
    max.iter.harmony,
    max.iter.cluster,
    epsilon.cluster,
    epsilon.harmony,
    plot_convergence,
    TRUE,
    verbose,
    reference_values
  )

  harmonyEmbed <- t(as.matrix(harmonyObject$Z_corr))
  rownames(harmonyEmbed) <- row.names(embedding)
  colnames(harmonyEmbed) <- paste0(reduction.save, "_", seq_len(ncol(harmonyEmbed)))

  harmonyClusters <- t(harmonyObject$R)
  rownames(harmonyClusters) <- row.names(embedding)
  colnames(harmonyClusters) <- paste0('R', seq_len(ncol(harmonyClusters)))
  
  suppressWarnings({
    harmonydata <- Seurat::CreateDimReducObject(
      embeddings = harmonyEmbed,
      stdev = as.numeric(apply(harmonyEmbed, 2, stats::sd)),
      assay = assay.use,
      key = reduction.save,
      misc=list(R=harmonyClusters)
    )
  })

  object[[reduction.save]] <- harmonydata
  if (project.dim) {
    object <- Seurat::ProjectDim(
      object,
      reduction = reduction.save,
      overwrite = TRUE,
      verbose = FALSE
    )
  }
  return(object)
}

environment(RunHarmony.Seurat) <- environment(harmony::HarmonyMatrix)

RunUMAP2 <- function (object, reduction.key = "UMAP_", assay = NULL, reduction.model = NULL, 
    return.model = FALSE, umap.method = "uwot", n.neighbors = 30L, 
    n.components = 2L, metric = "cosine", n.epochs = NULL, learning.rate = 1, 
    min.dist = 0.3, spread = 1, set.op.mix.ratio = 1, local.connectivity = 1L, 
    repulsion.strength = 1, negative.sample.rate = 5, a = NULL, 
    b = NULL, uwot.sgd = FALSE, seed.use = 42, metric.kwds = NULL, 
    angular.rp.forest = FALSE, verbose = TRUE, ...) 
{
    CheckDots(...)
    if (!is.null(x = seed.use)) {
        set.seed(seed = seed.use)
    }
    if (umap.method != "umap-learn" && getOption("Seurat.warn.umap.uwot", 
        TRUE)) {
        warning("The default method for RunUMAP has changed from calling Python UMAP via reticulate to the R-native UWOT using the cosine metric", 
            "\nTo use Python UMAP via reticulate, set umap.method to 'umap-learn' and metric to 'correlation'", 
            "\nThis message will be shown once per session", 
            call. = FALSE, immediate. = TRUE)
        options(Seurat.warn.umap.uwot = FALSE)
    }
    if (umap.method == "uwot-learn") {
        warning("'uwot-learn' is deprecated. Set umap.method = 'uwot' and return.model = TRUE")
        umap.method <- "uwot"
        return.model <- TRUE
    }
    if (return.model) {
        if (verbose) {
            message("UMAP will return its model")
        }
        umap.method = "uwot"
    }
    if (inherits(x = object, what = "Neighbor")) {
        object <- list(idx = Indices(object), dist = Distances(object))
    }
    if (!is.null(x = reduction.model)) {
        if (verbose) {
            message("Running UMAP projection")
        }
        umap.method <- "uwot-predict"
    }
    umap.output <- switch(EXPR = umap.method, `umap-learn` = {
        if (!py_module_available(module = "umap")) {
            stop("Cannot find UMAP, please install through pip (e.g. pip install umap-learn).")
        }
        if (!is.null(x = seed.use)) {
            py_set_seed(seed = seed.use)
        }
        if (typeof(x = n.epochs) == "double") {
            n.epochs <- as.integer(x = n.epochs)
        }
        umap_import <- import(module = "umap", delay_load = TRUE)
        umap <- umap_import$UMAP(n_neighbors = as.integer(x = n.neighbors), 
            n_components = as.integer(x = n.components), metric = metric, 
            n_epochs = n.epochs, learning_rate = learning.rate, 
            min_dist = min.dist, spread = spread, set_op_mix_ratio = set.op.mix.ratio, 
            local_connectivity = local.connectivity, repulsion_strength = repulsion.strength, 
            negative_sample_rate = negative.sample.rate, a = a, 
            b = b, metric_kwds = metric.kwds, angular_rp_forest = angular.rp.forest, 
            verbose = verbose)
        umap$fit_transform(as.matrix(x = object))
    }, uwot = {
        if (metric == "correlation") {
            warning("UWOT does not implement the correlation metric, using cosine instead", 
                call. = FALSE, immediate. = TRUE)
            metric <- "cosine"
        }
        if (is.list(x = object)) {
            umap(X = NULL, nn_method = object, n_threads = nbrOfWorkers(), 
                n_components = as.integer(x = n.components), 
                metric = metric, n_epochs = n.epochs, learning_rate = learning.rate, 
                min_dist = min.dist, spread = spread, set_op_mix_ratio = set.op.mix.ratio, 
                local_connectivity = local.connectivity, repulsion_strength = repulsion.strength, 
                negative_sample_rate = negative.sample.rate, 
                a = a, b = b, fast_sgd = uwot.sgd, verbose = verbose, 
                ret_model = return.model)
        } else {
            umap(X = object, n_threads = nbrOfWorkers(), n_neighbors = as.integer(x = n.neighbors), 
                n_components = as.integer(x = n.components), 
                metric = metric, n_epochs = n.epochs, learning_rate = learning.rate, 
                min_dist = min.dist, spread = spread, set_op_mix_ratio = set.op.mix.ratio, 
                local_connectivity = local.connectivity, repulsion_strength = repulsion.strength, 
                negative_sample_rate = negative.sample.rate, 
                a = a, b = b, fast_sgd = uwot.sgd, verbose = verbose, 
                ret_model = return.model)
        }
    }, `uwot-predict` = {
        if (metric == "correlation") {
            warning("UWOT does not implement the correlation metric, using cosine instead", 
                call. = FALSE, immediate. = TRUE)
            metric <- "cosine"
        }
        if (is.null(x = reduction.model) || !inherits(x = reduction.model, 
            what = "DimReduc")) {
            stop("If running projection UMAP, please pass a DimReduc object with the model stored to reduction.model.", 
                call. = FALSE)
        }
        model <- Misc(object = reduction.model, slot = "model")
        if (length(x = model) == 0) {
            stop("The provided reduction.model does not have a model stored. Please try running umot-learn on the object first", 
                call. = FALSE)
        }
        if (is.list(x = object)) {
            uwot::umap_transform(X = NULL, nn_method = object, 
                model = model, n_threads = nbrOfWorkers(), n_epochs = n.epochs, 
                verbose = verbose)
        } else {
            umap_transform(X = object, model = model, n_threads = nbrOfWorkers(), 
                n_epochs = n.epochs, verbose = verbose)
        }
    }, stop("Unknown umap method: ", umap.method, call. = FALSE))
    if (return.model) {
#         umap.output$nn_index <- NULL
        umap.model <- umap.output
        umap.output <- umap.output$embedding
    }
    colnames(x = umap.output) <- paste0(reduction.key, 1:ncol(x = umap.output))
    if (inherits(x = object, what = "dist")) {
        rownames(x = umap.output) <- attr(x = object, "Labels")
    }
    else if (is.list(x = object)) {
        rownames(x = umap.output) <- rownames(x = object$idx)
    }
    else {
        rownames(x = umap.output) <- rownames(x = object)
    }
    umap.reduction <- CreateDimReducObject(embeddings = umap.output, 
        key = reduction.key, assay = assay, global = TRUE)
    if (return.model) {
        Misc(umap.reduction, slot = "model") <- umap.model
    }
    return(umap.reduction)
}


environment(RunUMAP2) <- environment(Seurat:::RunUMAP.default)

mapQuery <- function (exp_query, metadata_query, ref_obj, vars = NULL, verbose = TRUE, 
    do_normalize = TRUE, do_umap = TRUE, sigma = 0.1, return_type = c('symphony', 'Seurat')) 
{
    if (return_type == 'Seurat') {
        que <- Seurat::CreateSeuratObject(
            counts=exp_query,
            meta.data=metadata_query,
            assay='SymphonyQuery'
        )        
    }
    
    if (do_normalize) {
        if (verbose) 
            message("Normalizing")
        exp_query = normalizeData(exp_query, 10000, "log")
    }
    if (verbose) 
        message("Scaling and synchronizing query gene expression")
    idx_shared_genes = which(ref_obj$vargenes$symbol %in% rownames(exp_query))
    shared_genes = ref_obj$vargenes$symbol[idx_shared_genes]
    if (verbose) 
        message("Found ", length(shared_genes), " reference variable genes in query dataset")
    exp_query_scaled = scaleDataWithStats(exp_query[shared_genes, 
        ], ref_obj$vargenes$mean[idx_shared_genes], ref_obj$vargenes$stddev[idx_shared_genes], 
        1)
    exp_query_scaled_sync = matrix(0, nrow = length(ref_obj$vargenes$symbol), 
        ncol = ncol(exp_query))
    exp_query_scaled_sync[idx_shared_genes, ] = exp_query_scaled
    rownames(exp_query_scaled_sync) = ref_obj$vargenes$symbol
    colnames(exp_query_scaled_sync) = colnames(exp_query)
    if (verbose) 
        message("Project query cells using reference gene loadings")
    Z_pca_query = t(ref_obj$loadings) %*% exp_query_scaled_sync
    if (verbose) 
        message("Clustering query cells to reference centroids")
    Z_pca_query_cos = cosine_normalize_cpp(Z_pca_query, 2)
    R_query = soft_cluster(ref_obj$centroids, Z_pca_query_cos, 
        sigma)
    if (verbose) 
        message("Correcting query batch effects")
    if (!is.null(vars)) {
        design = droplevels(metadata_query)[, vars] %>% as.data.frame()
        onehot = design %>% purrr::map(function(.x) {
            if (length(unique(.x)) == 1) {
                rep(1, length(.x))
            }
            else {
                stats::model.matrix(~0 + .x)
            }
        }) %>% purrr::reduce(cbind)
        Xq = cbind(1, intercept = onehot) %>% t()
    }
    else {
        Xq = Matrix(rbind(rep(1, ncol(Z_pca_query)), rep(1, ncol(Z_pca_query))), 
            sparse = TRUE)
    }
    Zq_corr = moe_correct_ref(as.matrix(Z_pca_query), as.matrix(Xq), 
        as.matrix(R_query), as.matrix(ref_obj$cache[[1]]), as.matrix(ref_obj$cache[[2]]))
    colnames(Z_pca_query) = row.names(metadata_query)
    rownames(Z_pca_query) = paste0("PC_", seq_len(nrow(Zq_corr)))
    colnames(Zq_corr) = row.names(metadata_query)
    rownames(Zq_corr) = paste0("harmony_", seq_len(nrow(Zq_corr)))
    umap_query = NULL
    if (do_umap & !is.null(ref_obj$save_uwot_path)) {
        if (verbose) 
            message("UMAP")
        ref_umap_model = uwot::load_uwot(ref_obj$save_uwot_path, 
            verbose = FALSE)
        
        ## UMAP may have been learned on subset of columns
        umap_query = uwot::umap_transform(t(Zq_corr)[, 1:ref_umap_model$norig_col], ref_umap_model)
#         umap_query = uwot::umap_transform(t(Zq_corr), ref_umap_model)
        colnames(umap_query) = c("UMAP1", "UMAP2")
        rownames(umap_query) <- row.names(metadata_query)
    }
    if (verbose) 
        message("All done!")
    
    if (return_type == 'Seurat') {
        que@assays$SymphonyQuery@data <- exp_query
        que@assays$SymphonyQuery@scale.data <- exp_query_scaled_sync
        que[['pca']] <- Seurat::CreateDimReducObject(
            embeddings = t(Z_pca_query),
            loadings = ref_obj$loadings, 
            stdev = as.numeric(apply(Z_pca_query, 1, stats::sd)),
            assay = 'SymphonyQuery',
            key = 'pca_'
        )
        que[['harmony']] <- Seurat::CreateDimReducObject(
            embeddings = t(Zq_corr),
            stdev = as.numeric(apply(Zq_corr, 1, stats::sd)),
            assay = 'SymphonyQuery',
            key = 'harmony_',
            misc=list(R=R_query)
        )
        que <- Seurat::ProjectDim(que, reduction = 'harmony', overwrite = TRUE, verbose = FALSE)
        if (do_umap) {
            que[['umap']] <- Seurat::CreateDimReducObject(
                embeddings = umap_query,
                assay = 'SymphonyQuery',
                key = 'umap_'
            )            
        }
        return(que)
    } else if (return_type == 'symphony') {
        return(list(Z = Zq_corr, Zq_pca = Z_pca_query, R = R_query, 
            Xq = Xq, umap = umap_query, meta_data = metadata_query))
    } else {
        stop(glue('The return type = \"{return_type}\" is not available.'))
    }
    
}

environment(mapQuery) <- environment(symphony::mapQuery)

knnPredict.Seurat <- function(query_obj, ref_obj, label_transfer, k = 5) 
{
    if (!label_transfer %in% colnames(ref_obj$meta_data)) {
        stop('Label \"{label_transfer}\" is not available in the reference metadata.')
    }
    knn_pred <- class::knn(t(ref_obj$Z_corr), Embeddings(query, 'harmony'), 
        ref$meta_data[[label_transfer]], k = k)
    query_obj@meta.data[[label_transfer]] <- knn_pred
    return(query_obj)
}
