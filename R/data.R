#' Log(CP10k+1) normalized counts matrix (genes by cells) for 10x PBMCs dataset for vignette.
#' 
#' @format: Sparse matrix (dgCMatrix): dimensions 1,764 genes by 1,200 cells
#' 
#' @source \url{https://support.10xgenomics.com/single-cell-gene-expression/datasets}
"pbmcs_exprs_small"

#' Metadata for 10x PBMCs dataset for vignette.
#' 
#' @format: A data frame with 1,200 cells and 7 metadata fields. 
#' \describe{
#'   \item{cell_id}{unique cell ID}
#'   \item{donor}{dataset (3pv1, 3pv2, or 5p)}
#'   \item{nUMI}{number of UMIs}
#'   \item{nGene}{number of genes}
#'   \item{percent_mito}{percent mito genes}
#'   \item{cell_type}{cell type assigned in Symphony publication}
#'   \item{cell_type_broad}{cell subtype assigned in Symphony publication}
#'
#' }
#' 
#' @source \url{https://support.10xgenomics.com/single-cell-gene-expression/datasets}
"pbmcs_meta_small"
