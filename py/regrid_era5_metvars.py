#!/usr/bin/env python3

import os
import argparse
import glob
import xarray as xr
from tqdm import tqdm

# Arguments to parse from command line
parser = argparse.ArgumentParser(description="Regrid ERA5 data to spatial \
    resolution of each GCM.")
parser.add_argument("--src-dir","-s",type=str)
parser.add_argument("--gcm-dir","-g",type=str)
parser.add_argument("--dest-dir","-d",type=str)
parser.add_argument("--gcm-list","-m",nargs="+")
parser.add_argument('--verbose','-v', action='store_true')              

args = parser.parse_args()

src_dir = args.src_dir
gcm_dir = args.gcm_dir
dst_dir = args.dest_dir
gcms = args.gcm_list
verbose = args.verbose

filelist = glob.glob(os.path.join(src_dir,"*nc"))
filelist.sort()

N = len(filelist) * len(gcms)

with tqdm(total=N,disable=not verbose) as pbar:

    for m in range(len(gcms)):

        gcm_m = gcms[m]

        gcm_file = glob.glob(os.path.join(gcm_dir,gcm_m,"*nc"))[0]
        gcm_ds = xr.load_dataset(gcm_file)

        gcm_lat = gcm_ds["lat"].values
        gcm_lon = gcm_ds["lon"].values

        del gcm_ds

        for f in range(len(filelist)):

            era5_ds = xr.load_dataset(filelist[f],engine='h5netcdf',
                decode_cf=False)

            era5_interp = era5_ds.interp(lat=gcm_lat,lon=gcm_lon)
            era5_interp = era5_interp.fillna(-9999.)        

            era5_interp.to_netcdf(os.path.join(src_dir,"regridded_data",
                    "era5_regridded_to_" + gcm_m,
                    os.path.basename(filelist[f])))

            pbar.update()