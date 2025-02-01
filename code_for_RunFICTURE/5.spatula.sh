#!/bin/bash

# conda activate FICTURE

outpath=$1
sample_name=$2
n_factor=$3
train_width=$4
spatula_path=$5
scale_factor=$6

FICTURE_output_path=${outpath}/3.FICTURE/${sample_name}/analysis/nF${n_factor}.d_${train_width}
FICTURE_input_path=${outpath}/2.FormatChange/${sample_name}/transcripts_in_micron.tsv.gz
spatula_output_path=${outpath}/5.spatula/${sample_name}
prefix=nF${n_factor}.d_${train_width}.decode.prj_12.r_4_5

# original barcode in pixel-level coordinates
input_file=${outpath}/1.blacklist/${sample_name}_tissue_positions_FilterOut_blacklist.csv

awk -v scale=$scale_factor '
BEGIN {FS=","; OFS="\t"}
NR==1 {print "#barcode_idx\tX\tY"; next}
{print $1, $5*scale, $6*scale}
' $input_file | sort -t$'\t' -k3,3n | gzip > "${input_file%/*}/barcodes_in_micrometer.tsv.gz"
# !!! Important !!!
### if you want to sort by X-axis, change -k3,3n to -k2,2n

## assume that the output file names are transcripts_ficture_joined.tsv.gz
${spatula_path} join-pixel-tsv \
    --mol-tsv ${input_file%/*}/barcodes_in_micrometer.tsv.gz \
    --pix-prefix-tsv nF12__,${FICTURE_output_path}/${prefix}.pixel.sorted_by_major_axis.tsv.gz \
    --out-prefix ${spatula_output_path}/original_barcode_with_factor \
    --sort-axis Y \
    --max-dist-um 2

# --in-tsv : generic input of FICTURE. It should sorted by a major-axis
# --max-dist-um : The maximum theshold (in um) of distance between the spatial coordinates in transcripts and FICTURE's pixel-level decoding output to be considered as a match.
