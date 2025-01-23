# FICTURE_visiumHD
### Manually fix some errors with FICTURE. 
### Especially errors related to coordinates matching in visiumHD (and their visulization)

## requirements : 
### 1. jq / conda install conda-forge::jq
### 2. FICTURE / pip install ficture
### 3. Tissue debris or Systemic error region information (polygon) 
If you don't see any tissue debris or systemic error, please set NULL
### 4. Set path of bgzip, tabix
ex) bgzip_path=/home/sjcho/.conda/envs/FICTURE/bin/bgzip

### After annoying parameter set up and polygon set up.... just run RUN_FICTURE.sh
