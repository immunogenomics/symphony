#' Function to plot reference, colored by cell type
#' 
#' @param reference Symphony reference object (must have UMAP stored)
#' @param as.density if TRUE, plot as density; if FALSE, plot as individual cells
#' @param bins for density, nbins parameter for stat_density_2d
#' @param bandwidth for density, bandwidth parameter for stat_density_2d
#' @param title Plot title
#' @param color.by metadata column name for phenotype labels
#' @param celltype.colors custom color mapping
#' @param show.legend Show cell type legend
#' @param show.labels Show cell type labels
#' @param show.centroids Plot soft cluster centroid locations
#' @import ggplot2
#' @import ggrastr
#' @import RColorBrewer
#' @import ggrepel
#' @import uwot
#' @return A ggplot object.
#' @export
plotReference = function(reference,              # Symphony reference object
                          as.density = TRUE,      # if FALSE, plot as individual cells
                          bins = 10,              # for density, nbins parameter for stat_density_2d
                          bandwidth = 1.5,        # for density, bandwidth parameter for stat_density_2d
                          title = 'Reference',    # Plot title
                          color.by = 'cell_type', # metadata column name for cell type labels
                          celltype.colors = NULL, # custom color mapping
                          show.legend = TRUE,     # Show cell type legend
                          show.labels = TRUE,     # Show cell type labels
                          show.centroids = FALSE) { # Plot soft cluster centroid locations
    
    if (is.null(reference$umap)) {
        stop('Error: umap slot is empty. UMAP was not saved for this reference!')
    }
    
    umap_labels = cbind(reference$meta_data, reference$umap$embedding)
    
    p = umap_labels %>%
            dplyr::sample_frac(1L) %>% # permute rows randomly
            ggplot(aes(x = UMAP1, y = UMAP2))
    
    if (as.density) {
        # Plot as density
        p = p + stat_density_2d(geom = 'polygon', aes(alpha = ..level.., fill = get(color.by)), 
                    contour_var = "ndensity", bins = bins, h = bandwidth)
        if (!is.null(celltype.colors)) { p = p + scale_fill_manual(values = celltype.colors) + 
                                           labs(fill = color.by)}
    } else { 
        # Plot as individual points
        p = p + geom_point_rast(aes(col = get(color.by)), size = 0.3, stroke = 0.2, shape = 16)
        if (!is.null(celltype.colors)) { p = p + scale_color_manual(values = celltype.colors) + labs(color = color.by)}
    }
    
    # Default formatting
    p = p + theme_bw() +
            labs(title = title) + 
            theme(plot.title = element_text(hjust = 0.5)) +
            theme(legend.position="bottom") +
            theme(legend.text = element_text(size=8), legend.title=element_text(size=12)) + 
            guides(colour = guide_legend(override.aes = list(size = 4))) + guides(alpha = 'none')

    if (show.centroids) {
        # Add centroid locations
        centroids = reference$Z_corr %*% t(reference$R)
        ref_umap_model = uwot::load_uwot(reference$save_uwot_path, verbose = FALSE)
        umap_centroids = uwot::umap_transform(t(centroids), ref_umap_model) %>% as.data.frame()
        colnames(umap_centroids) = c('UMAP1', 'UMAP2')
        
        p = p + geom_point(data = umap_centroids, aes(x = UMAP1, y = UMAP2, fill = 'centroid'), size = 0.3)
    }
    
    if (show.labels) {
        # Add cell type labels (at median coordinate per cell type)
        labels.cent = umap_labels %>% 
            dplyr::group_by_at(color.by) %>% #group_by_at takes variable column name
            dplyr::select(UMAP1, UMAP2) %>% 
            dplyr::summarize_all(median)
        
        p = p + ggrepel::geom_text_repel(data = labels.cent, aes(x= UMAP1, y = UMAP2, label = get(color.by)), 
                    segment.alpha = 0.5, segment.size = 0.2, box.padding = 0.01, color = 'black')
    }
    
    if (!show.legend) {
        p = p + theme(legend.position="none")
    }
    return(p)
}
