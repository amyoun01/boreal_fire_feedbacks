#!/usr/bin/env python3

import os, sys
import glob
import numpy as np
import xarray as xr
import argparse
from tqdm import tqdm

abs_path = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(abs_path,'z_functions'))

import z_utils as h
from MBC_Rscript import MBC_Rscript

path_to_Rscript = os.path.join(os.path.split(abs_path)[0],"R/MBC_Rscript.R")

WDIR = os.getenv("WDIR")
GCM = os.getenv("GCM")
tmpdir = os.path.join(WDIR,"tmp")

# Arguments to parse from command line
parser = argparse.ArgumentParser(description="Quantile delta mapping by calling \
     MBC R package.")
parser.add_argument("--ref-dir","-r",type=str)
parser.add_argument("--gcm-dir","-g",type=str)
parser.add_argument("--dest-dir","-d",type=str)
parser.add_argument("--hstyr","-t",nargs=2,type=int)
parser.add_argument("--simyr","-s",nargs=2,type=int)
parser.add_argument("--keep-hst",action="store_true")
parser.add_argument("--verbose","-v", action="store_true")

args = parser.parse_args()

ref_dir = args.ref_dir
gcm_dir = args.gcm_dir
dst_dir = args.dest_dir
hstyr = args.hstyr
simyr = args.simyr
keephst = args.keep_hst
verbose = args.verbose

hstyr = np.arange(hstyr[0],hstyr[1]+1)
simyr = np.arange(simyr[0],simyr[1]+1)

ref_filelist = glob.glob(os.path.join(ref_dir,"*.nc"))
ref_filelist.sort()

gcm_filelist = glob.glob(os.path.join(gcm_dir,"*.nc"))
gcm_filelist.sort()

# -----------------------------------------------------------------------------
# Import reference ERA5 data for calibration period
# -----------------------------------------------------------------------------

ref_filename_parts = h.extract_str_parts(ref_filelist,extension='nc',sep='_')
ref_file_yr = ref_filename_parts[:,-1].astype('int')

ref_yr_id = np.flatnonzero(np.isin(ref_file_yr,hstyr))
ref_hst_import_files = list(np.array(ref_filelist)[ref_yr_id])

ref = xr.open_mfdataset(ref_hst_import_files,combine='nested',parallel=False,\
    engine='h5netcdf')

# -----------------------------------------------------------------------------
# Import reference GCM data for calibration period
# -----------------------------------------------------------------------------

gcm_filename_parts = h.extract_str_parts(gcm_filelist,extension='nc',sep='_')
gcm_file_yr = gcm_filename_parts[:,-1].astype('int')

hst_yr_id = np.flatnonzero(np.isin(gcm_file_yr,hstyr))
gcm_hst_import_files = list(np.array(gcm_filelist)[hst_yr_id])

hst = xr.open_mfdataset(gcm_hst_import_files,combine='nested',parallel=False,\
    engine='h5netcdf')

# -----------------------------------------------------------------------------
# Import reference GCM data for calibration period
# -----------------------------------------------------------------------------

sim_yr_id = np.flatnonzero(np.isin(gcm_file_yr,simyr))
gcm_sim_import_files = list(np.array(gcm_filelist)[sim_yr_id])

sim = xr.open_mfdataset(gcm_sim_import_files,combine='nested',parallel=False,\
    engine='h5netcdf')

# -----------------------------------------------------------------------------
# Convert reference and GCM data to numpy arrays
# -----------------------------------------------------------------------------

ref_arr = np.stack((ref['tasmax'].values,ref['hursmin'].values,\
    ref['sfcWind'].values,ref['pr'].values),axis=3)

hst_arr = np.stack((hst['tasmax'].values,hst['hursmin'].values,\
    hst['sfcWind'].values,hst['pr'].values),axis=3)

sim_arr = np.stack((sim['tasmax'].values,sim['hursmin'].values,\
    sim['sfcWind'].values,sim['pr'].values),axis=3)

rfile = os.path.join(tmpdir,"ref.npy")
hfile = os.path.join(tmpdir,"hst.npy")
sfile = os.path.join(tmpdir,"sim.npy")
qfile = os.path.join(tmpdir,"qdm.npy")

# -----------------------------------------------------------------------------
# Get time unit values (here months) from reference dataset. This should exactly
# match the hst and sim datasets as well.
# -----------------------------------------------------------------------------

time = ref.time.dt.month
time = time.values

tfile = os.path.join(tmpdir,"time.npy")
np.save(file=tfile,arr=time,allow_pickle=False)

hst_qdm = np.zeros(hst_arr.shape,dtype="single")
sim_qdm = np.zeros(sim_arr.shape,dtype="single")

_, n, m, _ = ref_arr.shape
N = n * m # Number of rows and columns in grid

with tqdm(total=N,disable=not verbose) as pbar:

    for i in range(N):

        r, c = np.unravel_index(i,(n,m))

        ref_arr_i = ref_arr[:,r,c,:]
        hst_arr_i = hst_arr[:,r,c,:]
        sim_arr_i = sim_arr[:,r,c,:]

        np.save(file=rfile,arr=ref_arr_i,allow_pickle=False)
        np.save(file=hfile,arr=hst_arr_i,allow_pickle=False)
        np.save(file=sfile,arr=sim_arr_i,allow_pickle=False)

        qdm = MBC_Rscript(rfile,hfile,sfile,qfile,tfile,path_to_Rscript)

        hst_qdm[:,r,c,:] = qdm["hst_bc"]
        sim_qdm[:,r,c,:] = qdm["sim_bc"]

        pbar.update()

hst_qdm = \
    xr.Dataset(
        data_vars={
            'tasmax': (['time','lat','lon'],hst_qdm[...,0],hst['tasmax'].attrs),
            'hursmin': (['time','lat','lon'],hst_qdm[...,1],hst['hursmin'].attrs),
            'sfcWind': (['time','lat','lon'],hst_qdm[...,2],hst['sfcWind'].attrs),
            'pr': (['time','lat','lon'],hst_qdm[...,3],hst['pr'].attrs),
            },
        coords={
            'time': ('time',hst['time'].values),
            'lat': ('lat',hst['lat'].values,hst['lat'].attrs),
            'lon': ('lon',hst['lon'].values,hst['lon'].attrs),
            },
        )

hst_tasmax = hst_qdm["tasmax"]
hst_hursmin = hst_qdm["hursmin"]
hst_sfcWind = hst_qdm["sfcWind"]
hst_pr = hst_qdm["pr"]                                

sim_qdm = \
    xr.Dataset(
        data_vars={
            'tasmax': (['time','lat','lon'],sim_qdm[...,0],sim['tasmax'].attrs),
            'hursmin': (['time','lat','lon'],sim_qdm[...,1],sim['hursmin'].attrs),
            'sfcWind': (['time','lat','lon'],sim_qdm[...,2],sim['sfcWind'].attrs),
            'pr': (['time','lat','lon'],sim_qdm[...,3],sim['pr'].attrs),
            },
        coords={
            'time': ('time',sim['time'].values),
            'lat': ('lat',sim['lat'].values,sim['lat'].attrs),
            'lon': ('lon',sim['lon'].values,sim['lon'].attrs),
            },
        )

sim_tasmax = sim_qdm["tasmax"]
sim_hursmin = sim_qdm["hursmin"]
sim_sfcWind = sim_qdm["sfcWind"]
sim_pr = sim_qdm["pr"]     

for y in range(len(simyr)):

    sim_tasmax.isel(time=sim_tasmax.time.dt.year==simyr[y]).\
        to_netcdf(os.path.join(dst_dir,"tasmax_%s_%d.nc" % (GCM,simyr[y])),
            engine="h5netcdf")

    sim_hursmin.isel(time=sim_hursmin.time.dt.year==simyr[y]).\
        to_netcdf(os.path.join(dst_dir,"hursmin_%s_%d.nc" % (GCM,simyr[y])),
            engine="h5netcdf")

    sim_sfcWind.isel(time=sim_sfcWind.time.dt.year==simyr[y]).\
        to_netcdf(os.path.join(dst_dir,"sfcWind_%s_%d.nc" % (GCM,simyr[y])),
            engine="h5netcdf")

    sim_pr.isel(time=sim_pr.time.dt.year==simyr[y]).\
        to_netcdf(os.path.join(dst_dir,"pr_%s_%d.nc" % (GCM,simyr[y])),
            engine="h5netcdf")

if keephst:

    for y in range(len(hstyr)):

        hst_tasmax.isel(time=hst_tasmax.time.dt.year==hstyr[y]).\
            to_netcdf(os.path.join(dst_dir,"tasmax_%s_%d.nc" % (GCM,hstyr[y])),
                engine="h5netcdf")

        hst_hursmin.isel(time=hst_hursmin.time.dt.year==hstyr[y]).\
            to_netcdf(os.path.join(dst_dir,"hursmin_%s_%d.nc" % (GCM,hstyr[y])),
                engine="h5netcdf")

        hst_sfcWind.isel(time=hst_sfcWind.time.dt.year==hstyr[y]).\
            to_netcdf(os.path.join(dst_dir,"sfcWind_%s_%d.nc" % (GCM,hstyr[y])),
                engine="h5netcdf")

        hst_pr.isel(time=hst_pr.time.dt.year==hstyr[y]).\
            to_netcdf(os.path.join(dst_dir,"pr_%s_%d.nc" % (GCM,hstyr[y])),
                engine="h5netcdf")

# End of script ---------------------------------------------------------------