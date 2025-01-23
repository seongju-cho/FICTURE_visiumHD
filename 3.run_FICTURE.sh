#!/bin/bash
set -e
set -u
set -o pipefail

# conda activate FICTURE

spaceRanger_output=$1
formatChange_path=$2
train_width=$3
n_factor=$4
n_jobs=$5
threads=$6
bgzip_path=$7
tabix_path=$8
output_path=$9

### Normal
mtx_path=${spaceRanger_output}/binned_outputs/square_002um/filtered_feature_bc_matrix # or raw_feature_bc_matrix

ffile=${mtx_path}/features_with_header.tsv.gz
#original_pixel_to_micron=$(jq '.microns_per_pixel' ${spaceRanger_output}/binned_outputs/square_002um/spatial/scalefactors_json.json) # read from json
MU_SCALE=1
# 4.652140847631875 is the microns per pixel read from the json file

ficture run_together \
    --in-tsv ${formatChange_path}/transcripts_in_micron.tsv.gz \
    --in-minmax ${formatChange_path}/coordinate_minmax.tsv \
    --in-feature ${ffile} \
    --out-dir ${output_path} \
    --train-width ${train_width} \
    --n-factor ${n_factor} \
    --major-axis Y \
    --mu-scale ${MU_SCALE} \
    --n-jobs ${n_jobs} \
    --plot-each-factor \
    --threads ${threads} \
    --all \
    --bgzip ${bgzip_path} \
    --tabix ${tabix_path}
    # --lda-plot-um-per-pixel 10000 \
    # --decode-plot-um-per-pixel 10000 \
    # --decode-sub-um-per-pixel 10000

echo "Finish FICTURE now manually draw the pixel-level figure"

# Remove BARCODE column which cause error (But I have no idea why this is happen)
FICTURE_output_inside=${output_path}/analysis/nF${n_factor}.d_${train_width}
color_table=${FICTURE_output_inside}/figure/nF${n_factor}.d_${train_width}.rgb.tsv
decode_prefix=nF${n_factor}.d_${train_width}.decode.prj_12.r_4_5 # nF : number of Factors, d_ : size of anchor
zcat ${FICTURE_output_inside}/${decode_prefix}.pixel.sorted.tsv.gz | sed '4s/\tBARCODE//' | gzip > ${FICTURE_output_inside}/${decode_prefix}.pixel.sorted_remove_BARCODE.tsv.gz

# plot all the pixels
ficture plot_pixel_full --input ${FICTURE_output_inside}/${decode_prefix}.pixel.sorted_remove_BARCODE.tsv.gz \
                        --color_table ${color_table} \
                        --output ${FICTURE_output_inside}/figure/${decode_prefix}.pixel.png \
                        --plot_um_per_pixel 2 --full

# get factor report
ficture factor_report --path ${FICTURE_output_inside} \
                     --pref ${decode_prefix} \
                     --color_table ${color_table} \

# plot all the pixels separately (like cells.highlight of DimPlot in Seurat R package)
ficture plot_pixel_single --input ${FICTURE_output_inside}/${decode_prefix}.pixel.sorted_remove_BARCODE.tsv.gz \
                        --output ${FICTURE_output_inside}/figure/sub/${decode_prefix}.pixel \
                        --plot_um_per_pixel 2 --full --all

echo "All the FICTURE process is done"
