# Methods from singlecellmethods ---------------------------------

normalizeData <- function(A, scaling_factor = 1e4, method) {
    if(!'dgCMatrix' %in% class(A)) A <- as(A, "dgCMatrix")
    
    if (method == "log") {
        A@x <- A@x / rep.int(Matrix::colSums(A), diff(A@p))
        A@x <- scaling_factor * A@x
        A@x <- log(1 + A@x)
    } else if (method == "fft") {
        A@x <- A@x / rep.int(Matrix::colSums(A), diff(A@p))
        A@x <- scaling_factor * A@x
        A@x <- sqrt(A@x) + sqrt(1 + A@x)
    } else if (method == "geneCLR") {
        A@x <- as.numeric(normalizeCLR_dgc(A@x, A@p, A@i, ncol(A), nrow(A), 1))        
    } else if (method == "cellCLR") {
        A@x <- as.numeric(normalizeCLR_dgc(A@x, A@p, A@i, ncol(A), nrow(A), 2))
    } else {
        stop(sprintf("ERROR: method %s not implemented", method))
    }

    return(A)
}

scaleData <- function(A, margin = 1, thresh = 10) {
    A <- as(A, "dgCMatrix")
    
    if (margin != 1) A <- t(A)
    
    res <- scaleRows_dgc(A@x, A@p, A@i, ncol(A), nrow(A), thresh)
    if (margin != 1) res <- t(res)
    row.names(res) <- row.names(A)
    colnames(res) <- colnames(A)
    return(res)
}


#' Scale data with given mean and standard deviations
#'
#' @param A expression matrix (genes by cells)
#' @param mean_vec vector of mean values
#' @param sd_vec vector of standard deviation values
#' @param margin 1 for row-wise calculation
#' @param thresh threshold to clip max values
#' @return A matrix of scaled expression values.
#' @export
scaleDataWithStats <- function(A, mean_vec, sd_vec, margin = 1, thresh = 10) {
    if (!"dgCMatrix" %in% class(A))
        A <- as(A, "dgCMatrix")
    
    if (margin != 1) A <- t(A)
    
    res <- scaleRowsWithStats_dgc(A@x, A@p, A@i, mean_vec, sd_vec, 
                                  ncol(A), nrow(A), thresh)
    if (margin != 1) res <- t(res)
    row.names(res) <- row.names(A)
    colnames(res) <- colnames(A)
    return(res)
}

#' Calculate standard deviations by row
#'
#' @param A expression matrix (genes by cells)
#' @param row_means row means
#' @param weights weights for weighted standard dev calculation
#' @return A vector of row standard deviations 
#' @export
rowSDs <- function(A, row_means=NULL, weights=NULL) {
    if (is.null(row_means)) {
        row_means <- rowMeans(A, weights)
    }
    if (is.null(weights)) {
        res <- as.numeric(rowSDs_dgc(A@x, A@p, A@i, row_means, ncol(A), nrow(A), TRUE))
    } else {
        res <- as.numeric(rowSDsWeighted_dgc(A@x, A@p, A@i, row_means, weights, ncol(A), nrow(A), TRUE))
    }
    names(res) <- row.names(A)
    return(res)
}

rowMeans <- function(A, weights=NULL) {
    if (is.null(weights)) {
        res <- Matrix::rowMeans(A)
    } else {
        res <- as.numeric(rowMeansWeighted_dgc(A@x, A@p, A@i, weights, ncol(A), nrow(A)))
    }
    names(res) <- row.names(A)
    return(res)
}

rowVarsStd <- function(A, row_means, row_sds, vmax, weights=NULL) {
    if (is.null(weights)) {
        res <- as.numeric(rowVarSDs_dgc(A@x, A@p, A@i, row_means, row_sds, vmax, ncol(A), nrow(A), FALSE))
    } 
#     else {
#         res <- as.numeric(rowSDsWeighted_dgc(A@x, A@p, A@i, row_means, weights, ncol(A), nrow(A), TRUE))
#     }
    names(res) <- row.names(A)
    return(res)
}

rowVars <- function(A, row_means=NULL, weights=NULL) {
    if (is.null(row_means)) {
        row_means <- rowMeans(A, weights)
    }
    if (is.null(weights)) {
        res <- as.numeric(rowSDs_dgc(A@x, A@p, A@i, row_means, ncol(A), nrow(A), FALSE))
    } else {
        res <- as.numeric(rowSDsWeighted_dgc(A@x, A@p, A@i, row_means, weights, ncol(A), nrow(A), FALSE))
    }
    names(res) <- row.names(A)
    return(res)
}

## columns are observations
soft_kmeans <- function(X, k, w, max_iter=20, sigma=0.1) {
    message('WARNING: soft_kmeans fxn uses cosine distance only')
    Z <- cosine_normalize_cpp(X, 2)
    if (missing(w))
    Y <- stats::kmeans(t(Z), centers = k, iter.max = 25, nstart = 10)$centers %>% t() ## D x K
    res <- soft_kmeans_cpp(Y, Z, max_iter, sigma)
    return(res)
}

# Symphony utils ---------------------------------

#' Predict annotations of query cells from the reference using k-NN method
#'
#' @param query_obj Symphony query object
#' @param ref_obj Symphony reference object
#' @param train_labels vector of labels to train
#' @param k number of neighbors
#' @param save_as string that result column will be named in query metadata
#' @param confidence return k-NN confidence scores (proportion of neighbors voting for the predicted annotation)
#' @param seed random seed (k-NN has some stochasticity in the case of ties) 
#' @return Symphony query object, with predicted reference labels stored in the 'save_as' slot of the query$meta_data
#' @export
knnPredict <- function(query_obj, ref_obj,
                       train_labels,                   # cell labels for k-NN classification 
                       k = 5,                          # number of reference neighbors
                       save_as = 'cell_type_pred_knn', # metadata column to save result
                       confidence = TRUE,              # return prediction confidence
                       seed = 0) {                     # random seed
    set.seed(seed)
    if (confidence) {
        knn_pred = class::knn(t(ref_obj$Z_corr), t(query_obj$Z), train_labels, k = k, prob = TRUE)
        knn_prob = attributes(knn_pred)$prob
        query_obj$meta_data[save_as] = knn_pred
        query_obj$meta_data[paste0(save_as, '_prob')] = knn_prob
    } else {
        knn_pred = class::knn(t(ref_obj$Z_corr), t(query_obj$Z), train_labels, k = k, prob = FALSE)
        query_obj$meta_data[save_as] = knn_pred
    }    
    return(query_obj)
}

## Predict cell type using k-NN method with cos distance
knnPredictCos <- function(query_obj, ref_obj,
                       train_labels, # cell labels for knn classification 
                       k = 5) {
    Z_ref_cos = cosine_normalize_cpp(ref_obj$Z_corr, 2)
    Z_query_cos = cosine_normalize_cpp(query_obj$Z, 2)
    knn_pred = class::knn(t(Z_ref_cos), t(Z_query_cos), train_labels, k = k)
    query_obj$meta_data$cell_type_pred_knn_cos = knn_pred
    return(query_obj)
}


#' Function for evaluating F1 by cell type,
#' adapted from automated cell type identifiaction benchmarking paper (Abdelaal et al. Genome Biology, 2019)
#' @param true vector of true labels
#' @param predicted vector of predicted labels
#' @return A list of results with confusion matrix ($Conf), median F1-score ($MedF1), F1 scores per class ($F1), and accuracy ($Acc).
#' @export
evaluate <- function(true, predicted) {
  "  
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
