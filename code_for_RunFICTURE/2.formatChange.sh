#!/bin/bash
set -e
set -u
set -o pipefail

spaceRanger_output=$1
blacklist_filtered_datas=$2
sample_name=$3
opath=$4

brc_path=${blacklist_filtered_datas}/${sample_name}_tissue_positions_FilterOut_blacklist.csv
whitelist=${blacklist_filtered_datas}/${sample_name}_whitelist.csv # is it whitelist?
mfile=${blacklist_filtered_datas}/${sample_name}_raw_counts_filter.mtx.gz

mtx_path=${spaceRanger_output}/binned_outputs/square_002um/filtered_feature_bc_matrix # or raw_feature_bc_matrix
ffile=${mtx_path}/features.tsv.gz

# Check coordinate range (for future record)
microns_per_pixel=$(jq '.microns_per_pixel' ${spaceRanger_output}/binned_outputs/square_002um/spatial/scalefactors_json.json) # read from json
coor=${opath}/coordinate_minmax.tsv

IFS=' ' read -r xmin xmax ymin ymax <<<$(cut -d',' -f 5-6 ${brc_path} | tail -n +2 | awk -v mu=${microns_per_pixel} -v FS=',' 'NR == 1 { xmin = $1; xmax = $1; ymin = $2; ymax = $2 } {\
if ($1 > xmax) { xmax = $1 }\
if ($1 < xmin) { xmin = $1 }\
if ($2 > ymax) { ymax = $2 }\
if ($2 < ymin) { ymin = $2 }\
} END { print xmin*mu, xmax*mu, ymin*mu, ymax*mu }')
echo -e "xmin\t${xmin}\nxmax\t${xmax}\nymin\t${ymin}\nymax\t${ymax}" > ${coor}
echo "Create Min Max Coordinate File"

# Add header to the features file
raw_features=${mtx_path}/features.tsv.gz
edited_features=${mtx_path}/features_with_header.tsv.gz
zcat ${raw_features} | awk 'BEGIN {print "ensemble\tgene\tgene_expression"} {print $0}' | gzip > ${edited_features}

echo "Successfully Create features with header"

# First match barcode index (in the matrix file) with their spatial coordinates (in the tissue_positions file)
awk -v FS=',' -v OFS='\t' 'NR==FNR{bcd[$1]=NR; next} ($1 in bcd){ printf "%d\t%.2f\t%.2f\n", bcd[$1], $2, $3 } ' <(cat $whitelist) <(cut -d',' -f 1,5,6 ${brc_path}) | sort -k1,1n > ${opath}/barcodes.tsv
echo "Successfully Create Barcode File"

# Then annotate the coordinates and gene ID to the matrix file (assume matrix.mtx.gz is sorted by the pixel indices, which seems to be always true)
awk 'BEGIN{FS=OFS="\t"} NR==FNR{ft[NR]=$1 FS $2; next} ($4 in ft) {print $1, $2, $3, ft[$4], $5 }' \
<(zcat $ffile) \
<(\
join -t $'\t' -1 1 -2 2 ${opath}/barcodes.tsv <(zcat $mfile | tail -n +4 | sed 's/ /\t/g') \
) | sort -k3,3n -k2,2n | sed '1 s/^/#barcode_idx\tX\tY\tgene_id\tgene\tCount\n/' | gzip -c > ${opath}/transcripts.tsv.gz

zcat ${opath}/transcripts.tsv.gz | awk -v mu=${microns_per_pixel} 'NR==1{print;next} {$2=$2*mu; $3=$3*mu; print}' OFS='\t' | gzip > ${opath}/transcripts_in_micron.tsv.gz
zcat ${opath}/transcripts_in_micron.tsv.gz | head

echo "Successfully Create Transcript File and convert pixel to micron"
