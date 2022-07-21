#!/usr/bin/env Python3

# Calculate following statistics:

# Max for entire year (FWI_max)
# Max 90-day running mean per year (FWI_fs)
# Fire weather season length (FWI_fwsl)
# Number of days >90th percentile of fire weather (FWI_95d)

import os,sys
import glob
import numpy as np
import xarray as xr

os.chdir("/Users/adam/projects/boreal_fire_feedbacks/code/finished/py/z_functions")
import z_utils as h

src_dir = "/Users/adam/projects/boreal_fire_feedbacks/processed_data/cffdrs/cmip6/ACCESS-CM2/bias_corrected"
dst_dir = "/Users/adam/projects/boreal_fire_feedbacks/results/cffdrs_stats/cmip6/ACCESS-CM2"
yrlims = [2020,2099]
hstyr = [1980,1999]

filelist = glob.glob(os.path.join(src_dir,"*nc"))
filelist_parts = os.path.basename(filelist[-1]).split("_")
mdl_name = filelist_parts[-2]

ds = xr.open_mfdataset(filelist).load()
ds = ds[["isi","bui","fwi"]]

min_20yr = ds.sel(time=slice("1980-01-01","1999-12-31")).min(dim="time")
max_20yr = ds.sel(time=slice("1980-01-01","1999-12-31")).max(dim="time")

prctile_95th = ds.sel(time=slice("1980-01-01","1999-12-31")).quantile(0.95,dim="time")

# --------------------------------------------------------------------------
# Cacluate annual max value of all indices for given year
# --------------------------------------------------------------------------
cffdrs_max = ds.groupby("time.year").max().load()

yr = np.arange(yrlims[0],yrlims[1]+1)

for i in range(len(yr)):

    ds_i = ds.isel(time=(ds.time.dt.year == yr[i])).load()

    # --------------------------------------------------------------------------
    # Cacluate average of 90-day moving window for each year
    # --------------------------------------------------------------------------
    cffdrs_fs_i = ds_i.rolling(time=90,center=True).mean()
    cffdrs_fs_i = cffdrs_fs_i.max(dim="time")

    # --------------------------------------------------------------------------
    # Cacluate fire weather season length
    # --------------------------------------------------------------------------
    cffdrs_fwsl_i = 100 * (ds_i - min_20yr) / (max_20yr - min_20yr)
    cffdrs_fwsl_i = (cffdrs_fwsl_i > 50).sum(dim="time")

    # --------------------------------------------------------------------------
    # Cacluate number of days > 95th percentile
    # --------------------------------------------------------------------------
    cffdrs_95d_i = (ds_i > prctile_95th).sum(dim="time")

    cffdrs_max_i = cffdrs_max.sel(year=1980)
    cffdrs_max_i = cffdrs_max_i.drop_vars("year")
    
    cffdrs_max_i = cffdrs_max_i.assign_coords(coords={"stat": "max"})
    cffdrs_fs_i = cffdrs_fs_i.assign_coords(coords={"stat": "fs"})
    cffdrs_95d_i = cffdrs_fs_i.assign_coords(coords={"stat": "95d"})
    cffdrs_fwsl_i = cffdrs_fs_i.assign_coords(coords={"stat": "fwsl"})

    cffdrs_fw_stats = A = xr.combine_nested([cffdrs_max_i,
        cffdrs_fs_i,cffdrs_95d_i,cffdrs_fwsl_i],concat_dim="stat")

    fn = os.path.join(dst_dir,"cffdrs_stats_%s_%d.nc" % (mdl_name,yr[i]))
    cffdrs_fw_stats.to_netcdf(fn,engine="h5netcdf")
