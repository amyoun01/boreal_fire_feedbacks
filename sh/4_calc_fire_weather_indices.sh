#!/bin/bash

# -----------------------------------------------------------------------------
# Initialize workspace 
# -----------------------------------------------------------------------------

# Activate conda environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate pygeo38

# Specify working parent directory and add python scripts to current PATH
export WDIR=~/projects/boreal_fire_feedbacks
export PATH="$WDIR/code/finished/py:$PATH"

# List of GCMs
GCM_LIST=("ACCESS-CM2" "MPI-ESM1-2-HR" "MRI-ESM2-0" "MIROC6" "EC-Earth3")

ERAYR=(1980 2019)
SIMYR=(2020 2099)

# -----------------------------------------------------------------------------
# 1. Calculate CFFDRS for ERA5
# -----------------------------------------------------------------------------

SRC_DIR=$WDIR/processed_data/climate/era5
DEST_DIR=$WDIR/processed_data/cffdrs/era5

cffdrs_calculate.py -s $SRC_DIR -d $DEST_DIR -t ${ERAYR[@]} --verbose

# -----------------------------------------------------------------------------
# 2. Calculate CFFDRS for bias-corrected CMIP6 output
# -----------------------------------------------------------------------------

for GCM in ${GCM_LIST[@]}; do
  
  SRC_DIR=$WDIR/processed_data/climate/cmip6/$GCM/bias_corrected
  DEST_DIR=$WDIR/processed_data/cffdrs/cmip6/$GCM/bias_corrected
  
  echo # Blank line
  echo "Calculating fire weather indices for ${GCM} ..."

  cffdrs_calculate.py -s $SRC_DIR -d $DEST_DIR -t ${SIMYR[@]} --verbose

  echo "Finished for ${GCM}!"

done

# Deactivate conda environment
conda deactivate