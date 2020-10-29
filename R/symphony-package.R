#' symphony
#' 
#' Efficient single-cell reference mapping
#' 
#' @name symphony
#' @docType package
#' @useDynLib symphony
#' @importFrom Rcpp loadModule
loadModule("symphony_module", TRUE)
NULL
