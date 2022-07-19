#!/bin/bash

# 2_organize_era5_data.sh
# -----------------------------------------------------------------------------
# Organize hourly ERA5 data by finding 12:00 local standard time for each
#   grid cell and calculating needed daily statistics for CFFDRS calculations:
#     -Maximum daily air temperature, total 24-hr precip, average surface wind
#      speed, and minium daily relative humidity.
# -----------------------------------------------------------------------------

# Activate conda environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate pygeo38

# Specify working parent directory and add python scripts to current PATH
export WDIR=~/projects/boreal_fire_feedbacks
export PATH="$WDIR/code/finished/py:$PATH"

# -----------------------------------------------------------------------------
# Organize and calculate daily statistics for ERA5 data
# -----------------------------------------------------------------------------

# Geographic extent of study area using corner coords in lat/long.
GEOLIMS=(-177.1875 35.0 -35.0 79.375) # (left bottom right top)

YRLIMS=(1996 2019) # Time period and years to process for ERA5.

# Directory where hourly ERA5 data are stored in NETCDF4 format.
SRC_DIR=$WDIR/data/climate/era5/netcdf4
# Directory where to store processed ERA5 data
DEST_DIR=$WDIR/processed_data/climate/era5

# Organize and summarize ERA5 data to get relevant daily statistics for CFFDRS
# calculations.
extract_era5_metvars.py -s $SRC_DIR -d $DEST_DIR -g ${GEOLIMS[@]} \
  -t ${YRLIMS[@]} --verbose

# Deactivate conda environment
conda deactivate

# End of script ---------------------------------------------------------------