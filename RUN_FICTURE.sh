#!/bin/bash
set -e
set -u
set -o pipefail

# conda activate FICTURE
# requirements : jq, bgzip, tabix, spatula

### Set Parameters

spaceRanger_output=/home/sjcho/projects/ATLAS_visiumHD/lung_HN00227429_20241003/1.preprocessing/Lung_tumor_HN00227429/outs
sample_name=Lung_Tumor
polygon_path=/home/sjcho/projects/ATLAS_visiumHD/lung_HN00227429_20241003/4.FICTURE/scripts_runTogether/Tumor_polygon.Rdata
outpath=/home/sjcho/projects/ATLAS_visiumHD/lung_HN00227429_20241003/4.FICTURE/outs_tumor
spatula_path=/home/sjcho/yard/apps/spatula/spatula/bin/spatula
scale_factor=4.620840292370345 # change pixel to micrometer in original barcode file

## FICTURE parameter
train_width=12 # size of anchor
n_factor=24 # number of factors
n_jobs=64 # number of jobs
threads=64 # number of threads
bgzip_path=/home/sjcho/.conda/envs/FICTURE/bin/bgzip
tabix_path=/home/sjcho/.conda/envs/FICTURE/bin/tabix

# build output directory
mkdir -p ${outpath}/1.blacklist
mkdir -p ${outpath}/2.FormatChange
mkdir -p ${outpath}/3.FICTURE/${sample_name}
mkdir -p ${outpath}/5.spatula/${sample_name}

# 1. Build a blacklist
/home/sjcho/.conda/envs/visiumHD/bin/Rscript /home/sjcho/projects/ATLAS_visiumHD/lung_HN00227429_20241003/4.FICTURE/scripts_runTogether/1.build_blacklist.R ${spaceRanger_output} ${sample_name} ${polygon_path} ${outpath}/1.blacklist

# 2. Filter out the blacklist from 
bash /home/sjcho/projects/ATLAS_visiumHD/lung_HN00227429_20241003/4.FICTURE/scripts_runTogether/2.formatChange.sh ${spaceRanger_output} ${outpath}/1.blacklist ${sample_name} ${outpath}/2.FormatChange

# 3. Run FICTURE
bash /home/sjcho/projects/ATLAS_visiumHD/lung_HN00227429_20241003/4.FICTURE/scripts_runTogether/3.run_FICTURE.sh ${spaceRanger_output} ${outpath}/2.FormatChange ${train_width} ${n_factor} ${n_jobs} ${threads} ${bgzip_path} ${tabix_path} ${outpath}/3.FICTURE/${sample_name}

# 4. sort FICTURE output by one axis / it has no separate storage path
bash /home/sjcho/projects/ATLAS_visiumHD/lung_HN00227429_20241003/4.FICTURE/scripts_runTogether/4.sort_ficture.sh ${outpath}/3.FICTURE/${sample_name} ${n_factor} ${train_width}

# 5. recover original coordinates information from FICTURE output by spatula
bash /home/sjcho/projects/ATLAS_visiumHD/lung_HN00227429_20241003/4.FICTURE/scripts_runTogether/5.spatula.sh ${outpath} ${sample_name} ${n_factor} ${train_width} ${spatula_path} ${scale_factor}
