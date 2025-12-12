# FICTURE_visiumHD
### Contact : seongju6787@postech.ac.kr

### Manually fix some errors with FICTURE. 
### Especially errors related to coordinates matching in visiumHD (and their visulization)

## requirements : 
### 1. jq / conda install conda-forge::jq
### 2. FICTURE / pip install ficture
### 3. spatula / check https://github.com/seqscope/spatula
### 4. Tissue debris or Systemic error region information (polygon) 
If you don't see any tissue debris or systemic error, please set NULL
Example code for set polygon is in the 1.build_blacklist.R file

### 5. Set path of bgzip, tabix
ex) bgzip_path=/home/sjcho/.conda/envs/FICTURE/bin/bgzip

### After annoying parameter set up and polygon set up.... just run RUN_FICTURE.sh
