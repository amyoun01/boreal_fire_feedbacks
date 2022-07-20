#!/usr/bin/env python3

# Import required libraries
import os, sys
import pickle
import argparse
import glob
import numpy as np
import xarray as xr
from tqdm import tqdm

abs_path = os.path.dirname(os.path.abspath(__file__))

# Import custom modules 
sys.path.append(os.path.join(abs_path,"z_functions"))
import cffdrs
import z_utils as h

standard_dict_file = os.path.join(abs_path,"z_ancillary_data",
    "standard_variable_names.pickle")
with open(standard_dict_file,'rb') as f:
    standard_variable_dict = pickle.load(f)

ffmc_attrs = standard_variable_dict["ffmc"]
dmc_attrs = standard_variable_dict["dmc"]
dc_attrs = standard_variable_dict["dc"]
isi_attrs = standard_variable_dict["isi"]
bui_attrs = standard_variable_dict["bui"]
fwi_attrs = standard_variable_dict["fwi"]

# Arguments to parse from command line
parser = argparse.ArgumentParser(description="Calculate CFFDRS indices.")
parser.add_argument("--src-dir","-s",type=str)
parser.add_argument("--dest-dir","-d",type=str)
parser.add_argument("--yrlims","-t",nargs=2,type=int)
parser.add_argument("--verbose","-v",action="store_true")

args = parser.parse_args()

src_dir = args.src_dir
dst_dir = args.dest_dir
yr = args.yrlims
verbose = args.verbose

yr = np.arange(yr[0],yr[1]+1)
N = len(yr)

with tqdm(total=N,disable=not verbose) as pbar:

    for i in range(N):

        glob_pattern = os.path.join(src_dir,"*%d*.nc" % yr[i])
        file_list = glob.glob(glob_pattern)

        if len(file_list) == 4:

            ds = xr.open_mfdataset(file_list,engine="h5netcdf")

        else:

            ds = xr.load_dataset(file_list,engine="h5netcdf")

        vars = list(ds.keys())
        stand_var_list = list()

        for v in range(len(vars)):

            stand_name = ds[vars[v]].attrs["standard_name"]
            stand_var_list.append(stand_name)

            if stand_name == ("air_temperature"):           
                tas = ds[vars[v]].values

            elif stand_name == ("relative_humidity"):
                hur = ds[vars[v]].values

            elif stand_name == ("wind_speed"):
                sfcWind = ds[vars[v]].values

            elif stand_name == ("precipitation"):
                pr = ds[vars[v]].values

        c = np.count_nonzero(np.isin(stand_var_list,["air_temperature",
            "relative_humidity","wind_speed","precipitation"]))

        if c != 4:
            print("!Not all variables present for CFFDRS calculations.")
            sys.exit(0)

        metvars = np.stack((tas,hur,sfcWind,pr),axis=-1)
        mon = ds.time.dt.month.values

        fwi_vals = cffdrs.cffdrs_calc(metvars,mon)    

        cffdrs_ds = xr.Dataset(
                data_vars=
                {"ffmc": (["time","lat","lon"],fwi_vals["ffmc"],ffmc_attrs),
                "dmc": (["time","lat","lon"],fwi_vals["dmc"],dmc_attrs),
                "dc": (["time","lat","lon"],fwi_vals["dc"],dc_attrs),
                "isi": (["time","lat","lon"],fwi_vals["isi"],isi_attrs),
                "bui": (["time","lat","lon"],fwi_vals["bui"],bui_attrs),
                "fwi": (["time","lat","lon"],fwi_vals["fwi"],fwi_attrs),
                },coords=
                {"time": ("time",ds["time"].values),
                "lat": ("lat",ds["lat"].values,ds["lat"].attrs),
                "lon": ("lon",ds["lon"].values,ds["lon"].attrs),
                })

        cffdrs_ds_dt = h.uconvert_time(cffdrs_ds["time"].values,
            "days since 1850-01-01 12:00:00")

        cffdrs_ds = cffdrs_ds.merge({"time": cffdrs_ds_dt},
            overwrite_vars="time")
        cffdrs_ds["time"].attrs["units"] = "days since 1850-01-01 12:00:00"
        cffdrs_ds["time"].attrs["calendar"] = "noleap"

        fn = os.path.basename(file_list[-1])
        fn = fn.split("_")
        fn[0] = "cffdrs"
        fn = "_".join(fn)

        fn = os.path.join(dst_dir,fn)

        cffdrs_ds.to_netcdf(fn,
            engine="h5netcdf",
            encoding={
                "time": {"dtype": "int","_FillValue": None},
                "lat":  {"dtype": "single","_FillValue": None},
                "lon": {"dtype": "single","_FillValue": None},
                "ffmc": {"dtype": "single","_FillValue": -9999},
                "dmc": {"dtype": "single","_FillValue": -9999},
                "dc": {"dtype": "single","_FillValue": -9999},
                "isi": {"dtype": "single","_FillValue": -9999},
                "bui": {"dtype": "single","_FillValue": -9999},
                "fwi": {"dtype": "single","_FillValue": -9999},
                })
        
        pbar.update() # Update progress bar
