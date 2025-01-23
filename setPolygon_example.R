### this script is for
### running default pipeline for visiumHD
### conda activate visiumHD
### based on https://satijalab.org/seurat/articles/visiumhd_analysis_vignette

library(Seurat)
library(ggplot2)
library(dplyr)
library(hdf5r)
library(magrittr)
library(data.table)
library(Matrix)

### Set Parameters
## save_path : where to save the results (please include last /)
## localdir_normal : 10x normal data directory
## localdir_tumor : 10x tumor data directory
## please set polygon parameter below manually!!!!!!!!

save_path = '/home/sjcho/projects/ATLAS_visiumHD/lung_HN00227429_20241003/4.FICTURE/outs/4.0.QC_visiumHD/'
localdir_normal <- "/home/sjcho/projects/ATLAS_visiumHD/lung_HN00227429_20241003/1.preprocessing/Lung_Normal_HN00227429/outs"
localdir_tumor <- '/home/sjcho/projects/ATLAS_visiumHD/lung_HN00227429_20241003/1.preprocessing/Lung_tumor_HN00227429/outs'

### 1. load Visium HD data
parquet_path_list = list(
    'Normal_02' = paste0(localdir_normal, '/binned_outputs/square_002um/spatial/tissue_positions.parquet'),
    'Tumor_02' = paste0(localdir_tumor, '/binned_outputs/square_002um/spatial/tissue_positions.parquet')
    )

# default is load filtered_feature_bc_matrix
# And I recommend use this one, because this matrix only contain genes from the 10x probe gene set
so_list <- list(
    'Normal_02' = Load10X_Spatial(data.dir = localdir_normal, bin.size = c(2)),
    'Tumor_02' = Load10X_Spatial(data.dir = localdir_tumor, bin.size = c(2))
)

lapply(names(so_list), function(sample_name) {
    so <- so_list[[sample_name]]
    p <- ggplot(GetTissueCoordinates(so), aes(x = x, y = y)) + 
        geom_point(size = 0.1) + 
        theme_classic()
  ggsave(p, file = paste0(save_path, sample_name, '_tissue_coordinates.png'), width = 12, height = 12)
})

### delete tissue debris and almost empty spots

# gene X spot # normal@assays$Spatial.02um$counts
# spot X coordinates # GetTissueCoordinates(normal) # Read10X_Coordinates(normal)

# Setting default assay changes between 2um and 8um binning

# spatial plot
for (sample_name in names(so_list)) {
    print(sample_name)
    so <- so_list[[sample_name]]
    DefaultAssay(so) <- "Spatial.002um"

    p <- SpatialPlot(so, cells.highlight = Cells(so)[so[['log_nCount']][[1]] == log10(1)], image.alpha = 1, facet.highlight = TRUE, pt.size.factor = 2.5) + theme(legend.position = "right")
    ggsave(p, file = paste0(save_path, sample_name, '_nCount_spatial_under_log10_1.png'))
    # so_list[[sample_name]] <- so
}

umi_cutoff = list(
    'Normal_02' = 0,
    'Tumor_02' = 0
)

### YOU CAN ONLY SET POLYGON BY MANUALLY IN THIS SCRIPT
# set polygon to filter out tissue debris, systemic artifacts

polygon_cutoff = list(
    'Normal_02' = NULL,
    'Tumor_02' = matrix(c(2150, 250, 2150, 1500, 1600, 1600, 1600, 2000, 2300, 2000, 2300, 250), ncol = 2, dimnames = list(NULL, c('x', 'y')), byrow = TRUE))


# visualize tissue morphology with polygon cutoff
for (sample_name in names(so_list)) {
    so <- so_list[[sample_name]]
    p <- ggplot(GetTissueCoordinates(so), aes(x = x, y = y)) + geom_point(size = 0.1) + theme_classic()
    if (length(polygon_cutoff[[sample_name]]) > 1) {
        polygon <- polygon_cutoff[[sample_name]]
        p <- p + geom_polygon(data = polygon, aes(x = x, y = y), 
                                    fill = NA, color = 'red')
    }
    ggsave(p, file = paste0(save_path, sample_name, '_tissue_coordinates_polygon.png'))
}
