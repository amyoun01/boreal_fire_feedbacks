#!/usr/bin/env python3

# ------------------------------------------------------------------------------
# clip_cmip.py
# ------------------------------------------------------------------------------
#
# This script re-organizes the downloaded netcdf files into a standardized 
# format for use in ultimately calculating CFFDRS indices. The primary 
# library/tool for this reorganization was 'xarray'.
#
# Standardization includes:
#     - Converting longitude from 0-thru-360 to -180-thru-180
#     - Subsetting datasets spatiall and temporally
#     - Converting to 365_day calendar for each GCM
#     - Setting time units to 'days since 1850-01-01'    
#     - Converting units on data variables ('e.g., K to C)
#     - Export these newly organized datasets at nc files
#
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Import required modules
# ------------------------------------------------------------------------------
import os, sys
import pickle
import argparse
import glob
import datetime as dt
import numpy as np 
import xarray as xr
from tqdm import tqdm

# ------------------------------------------------------------------------------
# Get directory name of current clip_cmip.py file
# Then use this directory name to import custom utils functions
# ------------------------------------------------------------------------------
abs_path = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(abs_path,"z_functions"))
import z_utils as h

# ------------------------------------------------------------------------------
# Load standardized dictionary for variables used ifn this project
# ------------------------------------------------------------------------------
standard_dict_file = os.path.join(abs_path,
                                "z_ancillary_data",
                                "standard_variable_names.pickle")
with open(standard_dict_file,"rb") as f:
    standard_variable_dict = pickle.load(f)

# ------------------------------------------------------------------------------
# Arguments to parse from command line
# ------------------------------------------------------------------------------
parser = argparse.ArgumentParser(description='Clip global CMIP6 output to user\
defined spatial and temporal limits.')
parser.add_argument('--src-dir','-s',type=str,
                default=os.getcwd())
parser.add_argument('--dest-dir','-d',type=str,
                default=os.getcwd())
parser.add_argument('--geo-lims','-g',type=float,nargs=4,
                default=[-180,-90,180,90])
parser.add_argument('--date-lims','-t',nargs=2,
                default=['1850-01-01','2100-01-01'])
parser.add_argument('--time-units',type=str,
                    default='days since 1850-01-01 12:00:00')
parser.add_argument('--unit-convert','-u',
                action='store_true')
parser.add_argument('--verbose','-v',
                action='store_true')

args = parser.parse_args()

geo_lims = args.geo_lims
date_lims = args.date_lims
time_units = args.time_units
src_dir = args.src_dir
dest_dir = args.dest_dir
unit_convert = args.unit_convert
verbose = args.verbose

# Convert str format date lims to datetime type
date_lims = [dt.datetime.strptime(date_lims[0],'%Y-%m-%d'), \
             dt.datetime.strptime(date_lims[1],'%Y-%m-%d')]

file_list = glob.glob(os.path.join(src_dir,'*nc')) # Get list of all netcdf 
                                                   # files 

# Get dates that each data file span, first splitting up filenames into 
# individual components.
filename_parts = h.extract_str_parts(file_list,extension='nc')

# Get last part of filenames where date info is stored and separate again 
# by '-'.
gcm_file_dates = list(filename_parts[:,6])
start_and_end_dates = h.extract_str_parts(gcm_file_dates,sep='-')

# Use overlap function to determine which files to import and process
file_id = h.overlap(start_and_end_dates,
                    date_lims,
                    format='%Y%m%d',
                    index_return=True)

# For each file that overlaps with time period of interest ...
for f in tqdm(range(len(file_id)),disable=not verbose): 

    i = file_id[f]

    # Use xarray to load in dataset
    ds = xr.load_dataset(file_list[i],
                        mode='r',
                        mask_and_scale=True)

    # Remove un-needed dimensions (height) and variable bounds
    if np.any(np.array(list(ds.coords)) == 'height'):
        ds = ds.drop('height')
    ds = ds.drop_vars(['time_bnds','lat_bnds','lon_bnds'])
   
    # From new xarray dataset assign key data descriptors to 
    # different variable names. Used throughout code below, 
    # including file writing.
    var_name = filename_parts[i,0] # Variable name
    gcm = filename_parts[i,2] # GCM
    experiment = filename_parts[i,3] # Experiment
    ensemble = filename_parts[i,4] # Ensemble member    

    # If cmip, convert longitude from 0 to 360 to -180 to 180
    ds = ds.assign_coords(lon=(((ds.lon + 180) % 360) - 180))
    ds = ds.sortby('lon')

    # Get parsed lat and lon lims for data subsetting
    latlim = [geo_lims[1],geo_lims[3]]
    lonlim = [geo_lims[0],geo_lims[2]]

    # Subset entire xarray dataset by time bounds and spatial limits
    ds = ds.sel(time=slice(date_lims[0].strftime('%Y-%m'),
                           date_lims[1].strftime('%Y-%m')),
                lat=slice(latlim[0],latlim[1]),
                lon=slice(lonlim[0],lonlim[1]))

    # Convert variables to units needed in CFFDRS calculations
    if unit_convert:
        
        ds = h.uconvert(ds)

    # Reassign attributes based on predetermined standardized
    # variable dictionary that was previously loaded.
    ds['lat'].attrs = standard_variable_dict['lat']
    ds['lon'].attrs = standard_variable_dict['lon']
    ds[var_name].attrs = standard_variable_dict[var_name]

    ds.attrs = h.get_global_attrs(ds,['variable_id',
                                      'mip_era',
                                      'parent_source_id',
                                      'experiment_id',
                                      'variant_label',
                                      'tracking_id'])

    # Go through each year and write new dataset to file
    yr = ds.time.dt.year
    unique_yr = np.unique(yr)

    for y in range(0,len(unique_yr)): # For each year ...

        # Subset dataset for current year
        export_ds = ds.isel(time=(ds.time.dt.year == unique_yr[y]))

        # Convert calendar to 'noleap' (i.e. 365_day). Easier to work
        # with in the long run.
        export_ds = h.convert_calendar(export_ds,'noleap')

        # Reassign time dimension to integers that indicate number of
        # days since 1850-01-01. Also update attributes for time dimension.
        dt_arr = export_ds['time'].values
        dt_new = h.uconvert_time(dt_arr,time_units)

        export_ds = export_ds.merge({'time': dt_new},overwrite_vars='time')
        export_ds['time'].attrs['units'] = time_units
        export_ds['time'].attrs['calendar'] = 'noleap'        

        # Construct name for file to write to
        export_fn_parts = [var_name,gcm,experiment,ensemble,str(unique_yr[y])]
        export_filename = '_'.join(export_fn_parts) + '.nc'

        export_filename = os.path.join(dest_dir,export_filename)

        # Write subsetted dataset to hdf5 file
        export_ds.to_netcdf(export_filename,
                            mode='w',
                            format='NETCDF4',
                            engine='h5netcdf',
                            encoding={'time': {'dtype': 'int',
                                               '_FillValue': None},
                                      'lat':  {'dtype': 'single',
                                               '_FillValue': None},
                                       'lon': {'dtype': 'single',
                                               '_FillValue': None},
                                        var_name: {'dtype': 'single',
                                                   '_FillValue': -9999.}})

# End of script ----------------------------------------------------------------