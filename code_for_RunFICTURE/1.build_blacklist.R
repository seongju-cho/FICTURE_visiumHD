### this script is for
### running default pipeline for visiumHD
### conda activate visiumHD
### based on https://satijalab.org/seurat/articles/visiumhd_analysis_vignette

library(Seurat)
library(ggplot2)
library(dplyr)
library(arrow)
library(hdf5r)
library(magrittr)
library(data.table)
library(Matrix)

### Set Parameters
## save_path : where to save the results (please include last /)
## localdir_normal : 10x normal data directory
## localdir_tumor : 10x tumor data directory
## please set polygon parameter below manually!!!!!!!!

args <- commandArgs(trailingOnly = TRUE)
spaceRanger_outs <- args[1 ]
sample_name <- args[2]
polygon_path <- args[3]
save_path <- args[4]; save_path = paste0(save_path, '/')

########
parquet_path = paste0(spaceRanger_outs, '/binned_outputs/square_002um/spatial/tissue_positions.parquet')
so <- Load10X_Spatial(data.dir = spaceRanger_outs, bin.size = c(2))
DefaultAssay(so) <- "Spatial.002um"

# plot empty Spots
p <- SpatialPlot(so, cells.highlight = Cells(so)[so[['nCount_Spatial.002um']][[1]] == 0], image.alpha = 1, facet.highlight = TRUE, pt.size.factor = 2.5) + theme(legend.position = "right")
ggsave(p, file = paste0(save_path, sample_name, '_emptySpot.png'))

### delete tissue debris and empty spots

# gene X spot # normal@assays$Spatial.02um$counts
# spot X coordinates # GetTissueCoordinates(normal) # Read10X_Coordinates(normal)

### YOU CAN ONLY SET POLYGON BY MANUALLY IN THIS SCRIPT
# set polygon to filter out tissue debris, systemic artifacts
# polygon = matrix(c(2150, 250, 2150, 1500, 1600, 1600, 1600, 2000, 2300, 2000, 2300, 250), ncol = 2, dimnames = list(NULL, c('x', 'y')), byrow = TRUE)

if(polygon_path == 'NULL') {
    polygon = NULL
} else {
    polygon <- local({load(polygon_path)
                      get(ls())})
}
 
# visualize tissue morphology with polygon cutoff
p <- ggplot(GetTissueCoordinates(so), aes(x = x, y = y)) + geom_point(size = 0.1) + theme_classic()
if (!is.null(polygon)) {
    p <- p + geom_polygon(data = polygon, aes(x = x, y = y), 
                                fill = NA, color = 'red')
}
ggsave(p, file = paste0(save_path, sample_name, '_emptySpot_with_polygon.png'))

##### Build Blacklist

# define emptySpot
emptySpot = Cells(so)[so$nCount_Spatial.002um == 0]

# define spots in polygon
polygon_filtered = c()
if (!is.null(polygon)) {
    tissue_coordinates <- GetTissueCoordinates(so)
    polygon_filtered = Cells(so)[(point.in.polygon(tissue_coordinates$x, tissue_coordinates$y, polygon[, 'x'], polygon[, 'y']))]
}

blacklist = c(emptySpot, polygon_filtered)

# save black list
df = as.data.frame(blacklist)
colnames(df) = 'barcode'
fwrite(df, file = paste0(save_path, sample_name, '_blacklist.csv'), col.names = F, sep = ',')
rm(df)
print('blacklist saved')

# save whitelist
whitelist = setdiff(Cells(so), blacklist)
df = as.data.frame(whitelist)
colnames(df) = 'barcode'
fwrite(df, file = paste0(save_path, sample_name, '_whitelist.csv'), col.names = F, sep = ',')
# system(paste0('gzip ', paste0(save_path, sample_name, '_whitelist.tsv')))
print('whitelist saved')

# save sparse matrix
writeMM(so@assays$Spatial.002um$counts[, whitelist], file = paste0(save_path, sample_name, '_raw_counts_filter.mtx'))
# we have to add %metadata_json: {"software_version": "spaceranger-3.0.1", "format_version": 2} to this mtx file
metadata_line <- '%metadata_json: {"software_version": "spaceranger-3.0.1", "format_version": 2}'
sed_command <- sprintf("sed -i '1a\\%s' %s", metadata_line, paste0(save_path, sample_name, '_raw_counts_filter.mtx'))
system(sed_command)
system(paste0('gzip ', paste0(save_path, sample_name, '_raw_counts_filter.mtx')))
print('whitelist only count matrix saved')

# save modified parquet file
# first run parquet-tools to convert parquet to csv
commands = c(
    paste0("brc_parq=", parquet_path),
    paste0("opath=", paste0(save_path, sample_name)),
    paste0("brc_raw=${opath}_tissue_positions.raw.csv"),
    "parquet-tools csv ${brc_parq} > ${brc_raw}"
)
system(paste0(commands, collapse = ' && '))
# read csv
tissue_positions = fread(paste0(save_path, sample_name, '_tissue_positions.raw.csv'))
tissue_positions[, 'in_tissue'] = 0
tissue_positions[tissue_positions[['barcode']] %in% whitelist, 'in_tissue'] = 1
fwrite(tissue_positions, file = paste0(save_path, sample_name, '_tissue_positions_FilterOut_blacklist.csv'), col.names = T, sep = ',')
print('whitelist only tissue positions saved')
