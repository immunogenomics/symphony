#' Function to find variable genes using mean variance relationship method
#' 
#' @importFrom methods as
#' @importFrom stats loess median na.omit quantile
#' @importFrom rlang .data
#'
#' @param X expression matrix
#' @param groups vector of groups
#' @param min_expr min expression cutoff
#' @param max_expr max expression cutoff
#' @param min_dispersion min dispersion cutoff
#' @param max_dispersion max dispersion cutoff
#' @param num.bin number of bins to use for scaled analysis
#' @param binning.method how bins are computed
#' @param return_top_n returns top n genes 
#' @return A data.frame of variable genes
#' @export
findVariableGenes <- function(X, groups, min_expr = .1, max_expr = Inf, 
                               min_dispersion = 0, max_dispersion = Inf, 
                               num.bin = 20, binning.method = "equal_width", return_top_n = 0) {
    #https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
    group <- gene_mean <- symbol <- gene_dispersion <- NULL # prevents R CMD check note
    
    ## TODO: check that groups are 0 indexed
    groups <- factor(groups)
    groups_int <- as.integer(factor(groups)) - 1
    groups_table <- table(groups_int)
    
    ## initially compute means in non-log space, to use in vmr function below
    means_nonlog <- exp_mean(X@x, X@p, X@i, ncol(X), nrow(X), groups_int, groups_table)
    colnames(means_nonlog) <- levels(groups)
    
    vmr <- log_vmr(X@x, X@p, X@i, ncol(X), nrow(X), means_nonlog, groups_int, groups_table)    
    colnames(vmr) <- levels(groups)    

    ## transform means to logspace and join means and VMR  
    vargenes_df <- dplyr::inner_join(
        means_nonlog %>% log1p %>% as_tibble() %>% 
            cbind(symbol = row.names(X)) %>% 
            tidyr::gather(group, gene_mean, -symbol),
        vmr %>% as_tibble() %>% 
            cbind(symbol = row.names(X)) %>% 
            tidyr::gather(group, gene_dispersion, -symbol), 
        by = c("symbol", "group")
    )
    
    if (num.bin > 0) {
        if (binning.method == "equal_width") {
            .breaks <- num.bin
        }
        else if (binning.method == "equal_frequency") {
            .breaks <- c(-1, quantile(vargenes_df$gene_mean[vargenes_df$gene_mean > 0], probs = seq(0, 1, length.out = num.bin)))
        }
        else {
            stop(paste0("Invalid selection: '", binning.method, "' for 'binning.method'."))
        }
        
        vargenes_df <- data.table(vargenes_df)[
            , .data$the_bin := cut(.data$gene_mean, .breaks), by = .data$group
        ][]
        
        vargenes_df <- data.table(vargenes_df)[
            , .data$gene_dispersion_scaled := scale(.data$gene_dispersion), by = c('the_bin', 'group')
        ][]
        
        vargenes_df <- data.table(vargenes_df)[, .data$the_bin := NULL][]
    }
    
    vargenes_df <- vargenes_df %>% 
        dplyr::arrange(-.data$gene_dispersion) %>% 
        subset(.data$gene_mean >= min_expr & .data$gene_mean <= max_expr) %>%
        subset(.data$gene_dispersion >= min_dispersion & .data$gene_dispersion <= max_dispersion) 

    return(vargenes_df)
#     if (return_top_n > 0) {
#         vargenes_union <- unique(data.table(vargenes_df)[, head(.SD, return_top_n), by = group][, symbol])
#         return(vargenes_union)
#     } else {
#         return(vargenes_df)
#     }
    
}

#' Function to find variable genes using variance stabilizing transform (vst) method
#'
#' @param object expression matrix
#' @param groups finds variable genes within each group then pools
#' @param topn Return top n genes
#' @param loess.span Loess span parameter used when fitting the variance-mean relationship
#' @return A data.frame of variable genes, with means and standard deviations.
#' @export
vargenes_vst <- function(object, groups, topn, loess.span = 0.3) {
    clip.max <- sqrt(ncol(object))

    N <- ncol(object)
    if (missing(groups)) {
        groups <- rep('A', N)
    }
    
    res <- split(seq_len(N), groups) %>% lapply(function(idx) {
        object_group <- object[, idx]
        ## row means
        hvf.info <- data.frame(mean = Matrix::rowMeans(object_group))

        ## row vars
        hvf.info$variance <- rowVars(object_group, hvf.info$mean)

        ## initialize
        hvf.info$variance.expected <- 0
        hvf.info$variance.standardized <- 0

        not.const <- hvf.info$variance > 0

        ## loess curve fit 
        suppressWarnings({
            fit <- loess(formula = log10(variance) ~ log10(mean), 
                data = hvf.info[not.const, ], span = loess.span)            
        })

        ## extract fitted variance 
        hvf.info$variance.expected[not.const] <- 10^fit$fitted

        ## get row standard deviations after clipping
        hvf.info$variance.standardized <- rowVarsStd(
            object_group, 
            hvf.info$mean, 
            sqrt(hvf.info$variance.expected), 
            clip.max
        )

        hvf.info <- hvf.info %>% 
            tibble::rownames_to_column('symbol') %>% 
            dplyr::arrange(-.data$variance.standardized) %>% 
            tibble::rowid_to_column('rank') %>% 
            transform(group = unique(groups[idx]))

        return(hvf.info)        
    })
    
    
    if (missing(topn)) {
        ## MODE 1: return table 
        res <- Reduce(rbind, res) %>% 
            dplyr::select(.data$group, .data$symbol, .data$rank, .data$everything())

        if (length(unique(res$group)) == 1) {
            res$group <- NULL
        }
    } else {
        ## MODE 2: return genes
        res <- Reduce(union, lapply(res, function(x) head(x, topn)$symbol))
    }
    return(res)
}
