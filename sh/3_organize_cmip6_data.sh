
#!/bin/bash

# 3_organize_cmip6_data.sh
# -----------------------------------------------------------------------------
# 
# Read, organize, and do initial analysis for ERA5 and GCM climate data.
# The following steps:
# 1. Clip spatial extent of GCM data to N. America, covering study area. Also
#    keep only years needed for analysis (1980-2099).
# 2. Re-grid ERA5 to spatial resolution of each GCM using linear interpolation.
# 3. Use Quantile Delta Mapping to bias correct each meteorological variable by
#    calling MBC R package. 
#
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Initialize workspace 
# -----------------------------------------------------------------------------

# Activate conda environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate pygeo38

# Specify working parent directory and add python scripts to current PATH
export WDIR=~/projects/boreal_fire_feedbacks
export PATH="$WDIR/code/finished/py:$PATH"

# List of GCMs to process and geographic extent of study area using corner
# coords.
GCM_LIST=("MPI-ESM1-2-HR" "MRI-ESM2-0" "MIROC6" "EC-Earth3")
GEOLIMS=(-177.1875 35.0 -35.0 79.375) # (left bottom right top)

# -----------------------------------------------------------------------------
# 1. Organize and "clip" GCM data spatially and temporally
# -----------------------------------------------------------------------------

DATELIMS=("1980-01-01" "2099-12-31") # Range of dates to process for GCMs

for GCM in ${GCM_LIST[@]}; do

  # Directory where raw GCM data are stored
  SRC_DIR=$WDIR/data/climate/cmip6/$GCM 

  # Directory where clipped GCM data should be stored.
  DEST_DIR=$WDIR/processed_data/climate/cmip6/$GCM

  # Clip CMIP5 GCMs using python script. Use --unit-convert flag to convert
  # from K to C, m/s to km/hr, and kg m-2 s-1 to mm/day for temperature, 
  # wind speed, and precip, respectively. For Relative humidity it caps max 
  # daily relhum at 100%. --verbose flag produces progress bar.
  clip_cmip.py -s $SRC_DIR -d $DEST_DIR -g ${GEOLIMS[@]} -t ${DATELIMS[@]} \
    --unit-convert --verbose

done

# -----------------------------------------------------------------------------
# 2. Re-grid through interpolation the ERA5 data to the resolution of each GCM
# -----------------------------------------------------------------------------

# Directory where processed ERA5 data are stored at original spatial scale. 
SRC_DIR=$WDIR/processed_data/climate/era5
# Parent directory of all processed GCM data, used to access netcdf files and 
# import grid resolution.
GCM_DIR=$WDIR/processed_data/climate/cmip6
# Parent directory of where to send regridded data. python script below accesses
# subfolders for each specific GCM (i.e., ./era5_regridded_to_[GCM_NAME])
DEST_DIR=$SRC_DIR/regridded_data

# --gcm-list flag is a array listing all GCMs to grid using this script.
regrid_era5_metvars.py -s $SRC_DIR -g $GCM_DIR -d $DEST_DIR \
  --gcm-list ${GCM_LIST[@]} --verbose

# -----------------------------------------------------------------------------
# 3. Bias correcting using quantile delta mapping for each GCM
# -----------------------------------------------------------------------------

HSTYR=(1980 1999) # Years used for calibration of bias correcting

for (( G=0; G<${#GCM_LIST[@]}; G++ )); do # For each GCM in GCM_LIST variable ...

  # Set and export current GCM name as a environmental variable.
  export GCM=${GCM_LIST[$G]} 
  
  # Directory where regridded ERA5 data for specific GCM are located.
  REF_DIR=$WDIR/processed_data/climate/era5/regridded_data/era5_regridded_to_$GCM
  # Raw GCM output.
  GCM_DIR=$WDIR/processed_data/climate/cmip6/$GCM
  # Where to write bias corrected GCM output.
  DEST_DIR=$WDIR/processed_data/climate/cmip6/$GCM/bias_corrected

  for S in $(seq 2020 20 2080); do # For each 20 year period in the future ...

    SIMYR=($S $(($S + 19))) # Get an array that list the start and end years

    echo # Blank line
    echo "Bias correcting $GCM for ${SIMYR[0]}-${SIMYR[1]} ..."

    # Python script to do bias correcting. It calls "MBC" R package. --keep-hst
    # flag indicates that bias corrected historical data is exported as well. 
    # --verbose flag produces a progress bar.
    MBCpkg_bias_correcting.py -r $REF_DIR -g $GCM_DIR -d $DEST_DIR -t ${HSTYR[@]} \
      -s ${SIMYR[@]} --keep-hst --verbose
    
    echo "... done!"

  done

done

# Deactivate conda environment
conda deactivate

# End of script ---------------------------------------------------------------