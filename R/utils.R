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