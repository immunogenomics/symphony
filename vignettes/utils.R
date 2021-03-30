fig.size <- function (height, width) {
    options(repr.plot.height = height, repr.plot.width = width)
}

# Colors for PBMCs
pbmc_colors = c("B" = "#66C2A5", 
              "DC" = "#FC8D62",
              "HSC" = "#8DA0CB",
              "MK" = "#E78AC3", 
              "Mono_CD14" = "#A6D854",
              "Mono_CD16" = "#f2ec72",
              "NK" = "#62AAEA", 
              "T_CD4" = "#D1C656",
              "T_CD8" = "#968763")

# Colors for pancreas
celltype.colors = c('alpha'="#ed2bb1",
                    'beta'="#239eb3",
                    'gamma'="#d1bfec",
                    'delta'= "#FF6347",
                    'stellate'="#11e38c",
                    'immune'="#812050",
                    'ductal'="#b2d27a",
                    'endothelial'="#4e2da6",
                    'acinar'="#f6bb86",
                    'schwann'="#115d52",
                    'epsilon'="#a1def0",
                    'mast'="#8fec2f")

# Colors for pancreas query donors (Baron et al., 2016)
querydonor.colors = c('human1' = '#b9dbf0',
                      'human2' = '#77a1ba',
                      'human3' = '#6c7ca8',
                      'human4' = '#364261',
                      'mouse1' = '#e68c8c',
                      'mouse2' = '#b35757')

# Colors for fetal liver hematopoeisis
group.colors = c(   'B cell'='#f2bd80',
                    'DC precursor'='#1d6d1f',
                    'DC1'='#8c3ba0',
                    'DC2'='#6533ed',
                    'Early Erythroid'='#83e3f0',
                    'Early lymphoid/T'='#fd5917',
                    'Endothelial cell'='#4f8c9d',
                    'Fibroblast'='#eb1fcb',
                    'Hepatocyte'='#f5cdaf',
                    'HSC_MPP'='#9698dc',
                    'ILC precursor'='#20f53d',
                    'Kupffer Cell'='#f283e3',
                    'Late Erythroid'='#ffb2be',
                    'Mast cell'='#f3d426',
                    'Megakaryocyte'='#5ebf72',
                    'MEMP'='#a67649',
                    'Mid Erythroid'='#2f5bb1',
                    'Mono-Mac'='#90a479',
                    'Monocyte'='#f6932e',
                    'Monocyte precursor'='#d59e9a',
                    'Neut-myeloid prog.'='#caf243',
                    'NK'='#38b5fc',
                    'pDC precursor'='#c82565',
                    'Pre pro B cell'='#d6061a',
                    'pre-B cell'='#e36f6f',
                    'pro-B cell'='#1dfee1',
                    'VCAM1+ EI macro.'='#506356',
                    'centroid' ='black')

# Custom ordering to match original author publication ordering of states
group.ordering = c("HSC_MPP", "Pre pro B cell", 'pro-B cell', 'pre-B cell', 'B cell',
            'ILC precursor', 'Early lymphoid/T', 'NK', 'Neut-myeloid prog.',
            'pDC precursor','DC precursor', 'DC1', 'DC2', 'Monocyte precursor', 'Monocyte', 
            'Mono-Mac', 'Kupffer Cell', 'VCAM1+ EI macro.', 'MEMP', 'Mast cell',
            'Megakaryocyte', 'Early Erythroid', 'Mid Erythroid', 'Late Erythroid',
            'Endothelial cell', 'Fibroblast', 'Hepatocyte') 


#' Basic function to plot cells, colored and faceted by metadata variables
#' 
#' @param metadata metadata, with UMAP labels in UMAP1 and UMAP2 slots
#' @param title Plot title
#' @param color.by metadata column name for phenotype labels
#' @param facet.by metadata column name for faceting
#' @param color.mapping custom color mapping
#' @param show.legend Show cell type legend

plotBasic = function(umap_labels,                # metadata, with UMAP labels in UMAP1 and UMAP2 slots
                        title = 'Query',         # Plot title
                        color.by = 'cell_type',  # metadata column name for coloring
                        facet.by = NULL,         # (optional) metadata column name for faceting
                        color.mapping = NULL,    # custom color mapping
                        legend.position = 'right') {  # Show cell type legend
    
    p = umap_labels %>%
            dplyr::sample_frac(1L) %>% # permute rows randomly
            ggplot(aes(x = UMAP1, y = UMAP2)) + 
            geom_point_rast(aes(col = get(color.by)), size = 0.3, stroke = 0.2, shape = 16)
        if (!is.null(color.mapping)) { p = p + scale_color_manual(values = color.mapping) }
    
    # Default formatting
    p = p + theme_bw() +
            labs(title = title, color = color.by) + 
            theme(plot.title = element_text(hjust = 0.5)) +
            theme(legend.position=legend.position) +
            theme(legend.text = element_text(size=8), legend.title=element_text(size=12)) + 
            guides(colour = guide_legend(override.aes = list(size = 4))) + guides(alpha = FALSE)

    if(!is.null(facet.by)) {
        p = p + facet_wrap(~get(facet.by)) +
                theme(strip.text.x = element_text(size = 12)) }    
    return(p)
}
