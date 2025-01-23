#!/bin/bash
set -e
set -u
set -o pipefail

# conda activate FICTURE
# requirements : 
# 1. jq / conda install conda-forge::jq
# 2. FICTURE / pip install ficture
# 3. Tissue debris or Systemic error region information (polygon) / If you don't see any tissue debris or systemic error, please set NULL

### Set Parameters

spaceRanger_output=/home/sjcho/projects/ATLAS_visiumHD/lung_HN00227429_20241003/1.preprocessing/Lung_Normal_HN00227429/outs
sample_name=Lung_Normal
polygon_path=NULL
outpath=/home/sjcho/projects/ATLAS_visiumHD/lung_HN00227429_20241003/4.FICTURE/outs_normal_2

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

# 1. Build a blacklist
/home/sjcho/.conda/envs/visiumHD/bin/Rscript /home/sjcho/projects/ATLAS_visiumHD/lung_HN00227429_20241003/4.FICTURE/scripts_runTogether/1.build_blacklist.R ${spaceRanger_output} ${sample_name} ${polygon_path} ${outpath}/1.blacklist

# 2. Filter out the blacklist from 
bash /home/sjcho/projects/ATLAS_visiumHD/lung_HN00227429_20241003/4.FICTURE/scripts_runTogether/2.formatChange.sh ${spaceRanger_output} ${outpath}/1.blacklist ${sample_name} ${outpath}/2.FormatChange

# 3. Run FICTURE
bash /home/sjcho/projects/ATLAS_visiumHD/lung_HN00227429_20241003/4.FICTURE/scripts_runTogether/3.run_FICTURE.sh ${spaceRanger_output} ${outpath}/2.FormatChange ${train_width} ${n_factor} ${n_jobs} ${threads} ${bgzip_path} ${tabix_path} ${outpath}/3.FICTURE/${sample_name}
