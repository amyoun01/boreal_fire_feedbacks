#!/usr/bin/env python3

# extract_era5_metvars.py
# ------------------------------------------------------------------------------
#
# This script reads in hourly ERA5 datasets and calculates relevant summary
# statistics at a daily scale for calculating CFFDRS indices. To accomplish this,
# a time zone map with UTC differences is loaded in and noon local standard time
# (12:00 LST) is identified for each grid cell, needed for CFFDRS calculations.
# The relevant summary statistics (e.g., 24 hr precip, minimum daily relative
# humidity) are calculated from these hourly data.

# This script uses the h5py to read in each day at a time and then all daily 
# summaries for each individual year are organized and exported using xarray.

# A few more things this script does:
#     - Subsets the geographic limits of the ERA5 data to a predefined bounding
#     box.
#     - Converts the calendar to "noleap" instead of "propleptic gregorian", 
#     simply removing Feb 29th from leap years. 
#     - Uses dewpoint and near-surface air temperature to calculate relative 
#     humidity.
#     -Converts:
#         * temperature to Celcius from Kelvin
#         * windspedd from m/s to km/h
#         * total precip from m to mm
#         * Contrains relative humidity to 0-100%
#
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Import required libraries
# ------------------------------------------------------------------------------
import os, sys
import pickle
import argparse
import numpy as np
import datetime as dt
import h5py
import xarray as xr
from tqdm import tqdm

# ------------------------------------------------------------------------------
# Load standardized dictionary for variables used in this project
# ------------------------------------------------------------------------------
abs_path = os.path.dirname(os.path.abspath(__file__))

# ------------------------------------------------------------------------------
# Import custom utils package
# ------------------------------------------------------------------------------
sys.path.append(os.path.join(abs_path,"z_functions"))
import z_utils as h

# ------------------------------------------------------------------------------
# Import dictionary for variable names and attributes
# ------------------------------------------------------------------------------
standard_dict_file = os.path.join(abs_path,
    "z_ancillary_data",
    "standard_variable_names.pickle")
with open(standard_dict_file,"rb") as f:
    standard_variable_dict = pickle.load(f)

# ------------------------------------------------------------------------------
# Arguments to parse from command line
# ------------------------------------------------------------------------------
parser = argparse.ArgumentParser(description="Clip global CMIP6 output to user \
    defined spatial and temporal limits.")
parser.add_argument("--src-dir","-s",type=str,
                default=os.getcwd())
parser.add_argument("--dest-dir","-d",type=str,
                default=os.getcwd())
parser.add_argument("--geo-lims","-g",type=float,nargs=4,
                default=[-180,-90,180,90])
parser.add_argument("--date-lims","-t",type=int,nargs=2,
                default=["1980","2019"])                
parser.add_argument("--verbose","-v",
                action="store_true")
                
args = parser.parse_args()
start_time = dt.datetime.strptime("1980-01-01","%Y-%m-%d")
src_dir = args.src_dir
dest_dir = args.dest_dir
geo_lims = args.geo_lims
yr_lims = args.date_lims
verbose = args.verbose

# ------------------------------------------------------------------------------
# Finish initializing variables needed for processing
# ------------------------------------------------------------------------------

# Get date range
yr = range(yr_lims[0],yr_lims[1]+1)
t0 = np.array(dt.datetime(1900,1,1,0,0),dtype="datetime64[s]")
fmt = "%Y-%m-%dT%H:%M:%S"

# Specify lat/long limits
lonlim = [geo_lims[0],geo_lims[2]]
latlim = [geo_lims[1],geo_lims[3]]

# Import UTC offset map to identify 12:00 LST for each grid cell. Based off of
# time zone map.
utc_offset_ds = xr.load_dataset(\
    os.path.join(abs_path,"z_ancillary_data","utc_offset_grid.nc"))
utc_offset_grid = utc_offset_ds["utc_offset"].values
lat = utc_offset_ds["lat"].values
lon = utc_offset_ds["lon"].values
unique_offset = np.unique(utc_offset_grid)

# Get range of values that will be processed, number of total days between first
# and last date specified by yr_lims argument
if verbose:
    
    print("Starting to process ERA5 data ...")

    tdiff = dt.datetime(yr_lims[1],12,31,12,0) - \
        dt.datetime(yr_lims[0],1,1,12,0)

    N = tdiff.days

with tqdm(total=N,disable=not verbose) as pbar: # for progress bar

    for y in range(len(yr)): # For each year ...

        os.chdir(src_dir) # Change to data source directory

        # --------------------------------------------------------------
        # First import a single file (here temperature) just to get the datetime
        # info. This will be the same for all variables.
        # --------------------------------------------------------------
        with h5py.File("%d_era5_reanalysis_2m_temperature.nc" % yr[y],"r") as f:

            # Import time data (# hours since 1900-01-01 00:00) and convert to 
            # float
            tdelta = f["time"][:]
            tdelta = tdelta.astype("timedelta64[h]")

            # Convert to timedelta values and add to 1900-01-01 00:00
            dates = t0 + tdelta

            hr_of_day = np.array(
                [dt.datetime.strptime(dates[i].astype("str"),fmt).hour 
                    for i in range(0,len(dates))],
                )
            hr_of_day = 12 - hr_of_day

            # Convert to ordinal date (i.e. "datenum") or number of days since
            datenum = \
                np.array(dates - np.array(dt.datetime(1850,1,1,0,0,0),
                    dtype="datetime64[s]"),
                    dtype="timedelta64[D]")
            datenum = datenum.astype("int")
            unique_datenum = np.unique(datenum)

            ndays = len(unique_datenum)

            time_attrs = {"units": "days since 1850-01-01 12:00:00",
                        "calendar": "proleptic_gregorian"}

        # --------------------------------------------------------------
        # Create empty arrays to store climate output
        # --------------------------------------------------------------
        tasmax = np.nan * np.ones((ndays,len(lat),len(lon)),dtype="single")
        hursmin = np.nan * np.ones((ndays,len(lat),len(lon)),dtype="single")
        sfcWind = np.nan * np.ones((ndays,len(lat),len(lon)),dtype="single")
        pr = np.nan * np.ones((ndays,len(lat),len(lon)),dtype="single")

        for d in range(len(unique_datenum)): # For each unique day of year

            # If it is the first or last day of year then import hourly data
            # from either December 31st of previous year or January 1st of 
            # following year and use these to compute daily statistics (e.g.
            # 24hr precip).
            if d == 0:

                with h5py.File("%d_era5_reanalysis_2m_temperature.nc" % \
                    (yr[y]-1),"r") as f:
                    tmp_prev = f["t2m"][-24:,:,:]
                with h5py.File("%d_era5_reanalysis_2m_dewpoint_temperature.nc" % \
                    (yr[y]-1),"r") as f:
                    dew_prev = f["d2m"][-24:,:,:]
                with h5py.File("%d_era5_reanalysis_10m_u_component_of_wind.nc" % 
                (yr[y]-1),"r") as f:
                    u10_prev = f["u10"][-24:,:,:]
                with h5py.File("%d_era5_reanalysis_10m_v_component_of_wind.nc" % \
                    (yr[y]-1),"r") as f:
                    v10_prev = f["v10"][-24:,:,:]
                with h5py.File("%d_era5_reanalysis_total_precipitation.nc" % \
                    (yr[y]-1),"r") as f:
                    ppt_prev = f["tp"][-24:,:,:]

            if d == len(unique_datenum) - 1:

                with h5py.File("%d_era5_reanalysis_2m_temperature.nc" % \
                    (yr[y]+1),"r") as f:
                    tmp_post = f["t2m"][:24,:,:]
                with h5py.File("%d_era5_reanalysis_2m_dewpoint_temperature.nc" % \
                    (yr[y]+1),"r") as f:
                    dew_post = f["d2m"][:24,:,:]
                with h5py.File("%d_era5_reanalysis_10m_u_component_of_wind.nc" % \
                    (yr[y]+1),"r") as f:
                    u10_post = f["u10"][:24,:,:]
                with h5py.File("%d_era5_reanalysis_10m_v_component_of_wind.nc" % \
                    (yr[y]+1),"r") as f:
                    v10_post = f["v10"][:24,:,:]
                with h5py.File("%d_era5_reanalysis_total_precipitation.nc" % \
                    (yr[y]+1),"r") as f:
                    ppt_post = f["tp"][:24,:,:]
                
            for u in unique_offset: # For each time zone ...

                id_day = unique_datenum[d] == datenum # Find the current day
                id_hr_of_day = u == hr_of_day # Find which hour of day 
                                              # corresponds to noon local time
                                              # based on UTC time zone info.
                id_utc_offset = u == utc_offset_grid # Find which grid cells
                                                     # correspond to current 
                                                     # time zone for 12:00.

                # Convert bool to indices to figure our which specific values are
                # for current hour and day. Else, read in 12 hour window around
                # current noon time.
                id_day_hr = np.flatnonzero(id_day & id_hr_of_day)

                # --------------------------------------------------------------
                # Read in air temperature data using h5py
                # --------------------------------------------------------------
                with h5py.File("%d_era5_reanalysis_2m_temperature.nc" % yr[y],
                    mode="r") as f:

                    # If it is the first or last day of year include  hourly data
                    # for preceding or post years.
                    if (d == 0) & (id_day_hr-6 < 0):

                        tmp = f["t2m"][np.arange(0,id_day_hr+6),:,:].astype("f4")
                        tmp = np.concatenate((tmp_prev[int(id_day_hr-6):,:,:],
                            tmp),axis=0)

                    elif (d == len(unique_datenum)-1) & (id_day_hr+6 > len(datenum)):

                        tmp = f["t2m"][np.arange(id_day_hr-6,len(datenum)),:,:].astype("f4")
                        tmp = np.concatenate((tmp,tmp_post[:int(id_day_hr+6 - len(datenum)),:,:]),
                            axis=0)

                    else:
                            
                        tmp = f["t2m"][np.arange(id_day_hr-6,id_day_hr+6),:,:].\
                            astype("single")

                    # Get scale and offset factor to convert values to proper 
                    # units.
                    missing_value = f["t2m"].attrs["missing_value"]
                    scale_factor = f["t2m"].attrs["scale_factor"]
                    add_offset = f["t2m"].attrs["add_offset"]
                    
                    tmp = add_offset + scale_factor * tmp

                    # Identify and classify nan values
                    nan_id = tmp == missing_value
                    tmp[nan_id] = np.nan

                # --------------------------------------------------------------
                # Do the same as for air temperature but with dewpoint 
                # temperature.
                # --------------------------------------------------------------
                with h5py.File("%d_era5_reanalysis_2m_dewpoint_temperature.nc" % yr[y],
                    mode="r") as f:

                    if (d == 0) & (id_day_hr-6 < 0):

                        dew = f["d2m"][np.arange(0,id_day_hr+6),:,:].astype("f4")
                        dew = np.concatenate((dew_prev[int(id_day_hr-6):,:,:],dew),
                                axis=0)

                    elif (d == len(unique_datenum)-1) & (id_day_hr+6 > len(datenum)):

                        dew = f["d2m"][np.arange(id_day_hr-6,len(datenum)),:,:].astype("f4")
                        dew = np.concatenate((dew,dew_post[:int(id_day_hr+6 - len(datenum)),:,:]),
                                axis=0)

                    else:
                            
                        dew = f["d2m"][np.arange(id_day_hr-6,id_day_hr+6),:,:].astype("f4")

                    missing_value = f["d2m"].attrs["missing_value"]
                    scale_factor = f["d2m"].attrs["scale_factor"]
                    add_offset = f["d2m"].attrs["add_offset"]
                    
                    dew = add_offset + scale_factor * dew

                    nan_id = dew == missing_value
                    dew[nan_id] = np.nan

                # --------------------------------------------------------------
                # Read in 10m wind speed for u-component
                # --------------------------------------------------------------
                with h5py.File("%d_era5_reanalysis_10m_u_component_of_wind.nc" % yr[y],
                    mode="r") as f:

                    if (d == 0) & (id_day_hr-6 < 0):

                        u10 = f["u10"][np.arange(0,id_day_hr+6),:,:].astype("f4")
                        u10 = np.concatenate((u10_prev[int(id_day_hr-6):,:,:],u10),
                            axis=0)

                    elif (d == len(unique_datenum)-1) & (id_day_hr+6 > len(datenum)):

                        u10 = f["u10"][np.arange(id_day_hr-6,len(datenum)),:,:].astype("f4")
                        u10 = np.concatenate((u10,u10_post[:int(id_day_hr+6 - len(datenum)),:,:]),
                            axis=0)

                    else:
                            
                        u10 = f["u10"][np.arange(id_day_hr-6,id_day_hr+6),:,:].astype("f4")

                    missing_value = f["u10"].attrs["missing_value"]
                    scale_factor = f["u10"].attrs["scale_factor"]
                    add_offset = f["u10"].attrs["add_offset"]
                    
                    u10 = add_offset + scale_factor * u10

                    nan_id = u10 == missing_value
                    u10[nan_id] = np.nan

                # --------------------------------------------------------------
                # Read in 10m wind speed for v-component
                # --------------------------------------------------------------
                with h5py.File("%d_era5_reanalysis_10m_v_component_of_wind.nc" % yr[y],
                    mode="r") as f:

                    if (d == 0) & (id_day_hr-6 < 0):

                        v10 = f["v10"][np.arange(0,id_day_hr+6),:,:].astype("f4")
                        v10 = np.concatenate((v10_prev[int(id_day_hr-6):,:,:],v10),
                            axis=0)

                    elif (d == len(unique_datenum)-1) & (id_day_hr+6 > len(datenum)):

                        v10 = f["v10"][np.arange(id_day_hr-6,len(datenum)),:,:].astype("f4")
                        v10 = np.concatenate((v10,v10_post[:int(id_day_hr+6 - len(datenum)),:,:]),
                            axis=0)

                    else:
                            
                        v10 = f["v10"][np.arange(id_day_hr-6,id_day_hr+6),:,:].\
                            astype("single")

                    missing_value = f["v10"].attrs["missing_value"]
                    scale_factor = f["v10"].attrs["scale_factor"]
                    add_offset = f["v10"].attrs["add_offset"]
                    
                    v10 = add_offset + scale_factor * v10

                    nan_id = v10 == missing_value
                    v10[nan_id] = np.nan
                
                # --------------------------------------------------------------
                # Read in 10m wind speed for total precipitation
                # --------------------------------------------------------------
                with h5py.File("%d_era5_reanalysis_total_precipitation.nc" % yr[y],
                    mode="r") as f:

                    if (d == 0) & (id_day_hr-23 < 0):

                        ppt = f["tp"][np.arange(0,id_day_hr+1),:,:].stype("f4")
                        ppt = np.concatenate((ppt_prev[int(id_day_hr-23):,:,:],ppt),
                            axis=0)

                    else:
    
                        ppt = f["tp"][np.arange(id_day_hr-23,id_day_hr+1,1),:,:].astype("f4")

                    missing_value = f["tp"].attrs["missing_value"]
                    scale_factor = f["tp"].attrs["scale_factor"]
                    add_offset = f["tp"].attrs["add_offset"]
                    
                    ppt = add_offset + scale_factor * ppt

                    nan_id = ppt == missing_value
                    ppt[nan_id] = np.nan

                # --------------------------------------------------------------
                # Convert meteorological data to units needed for analysis and 
                # then add them to empty arrays previously created.
                # --------------------------------------------------------------

                # Convert temperature data to proper units
                tmp = tmp - 273.15
                dew = dew - 273.15

                # Use air and dewpoint temperature to calculate relative 
                # humidity. Constrain relhum to max 100%.
                relhum = 100.0 * h.vp_calc(dew) / h.vp_calc(tmp) 
                relhum = np.where(relhum>100.0,99.9,relhum)

                # Convert wind speed from m/s to km/h
                ws = 3.6 * np.sqrt(u10**2 + v10**2)

                # Convert from m to mm
                ppt = 1000.0 * ppt                

                tasmax[d,id_utc_offset] = tmp.max(axis=0)[id_utc_offset]
                hursmin[d,id_utc_offset] = relhum.min(axis=0)[id_utc_offset]
                sfcWind[d,id_utc_offset] = ws.mean(axis=0)[id_utc_offset]
                pr[d,id_utc_offset] = ppt.sum(axis=0)[id_utc_offset]

            pbar.update() # update progress bar

        # ----------------------------------------------------------------------
        # Set global attributes and for individual variables.
        # ----------------------------------------------------------------------
        tasmax_attrs = standard_variable_dict["tasmax"]
        hursmin_attrs = standard_variable_dict["hursmin"]
        sfcWind_attrs = standard_variable_dict["sfcWind"]
        pr_attrs = standard_variable_dict["pr"]

        tasmax_attrs["Description"] = """Daily maximum temperature calculated 
        from 06:00-18:00 hours local standard time."""
        
        hursmin_attrs["Description"] = """Daily minimum relative humidity 
        calculated from 06:00-18:00 hours local standard time. Relative humidity 
        values for each hour were calculated from vapour pressure estimates 
        based on the equation: vapor_pressure = 0.611 * exp((17.502 * 
        temperature) / (243.97 + temperature). Relative humidity was 
        subsequently calculated for each hour as RH = 
        vapour_pressure(air_temperature) / vapor_pressure(dewpoint_temperature)."""

        sfcWind_attrs["Description"] = """Daily average wind speed calculated 
        from 06:00-18:00 hours local standard time from downloaded 10m u10 and 
        v10 wind speeds. Converted to km hr/1."""

        pr_attrs["Description"] = """24-hour total precipitation at 12:00 local 
        standard time."""

        lat_attrs = standard_variable_dict["lat"]
        lon_attrs = standard_variable_dict["lon"]

        global_attrs = {"source_id": "ERA5",
        "url": "https://www.ecmwf.int/en/forecasts/datasets/reanalysis-datasets/era5",
        "history": "Original data downloaded Jan 7-13,2022 in netcdf format. \
            Converted to netcdf4 using netCDF library version 4.8.1 on Mar 30, 2022."}

        # ----------------------------------------------------------------------
        # Organize all data for current year into xarray dataset and clip to 
        # predefined geographic limits.
        # ----------------------------------------------------------------------        
        export_ds = xr.Dataset(
            data_vars = {
                "tasmax": (["time","lat","lon"],tasmax,tasmax_attrs),
                "hursmin": (["time","lat","lon"],hursmin,hursmin_attrs),
                "sfcWind": (["time","lat","lon"],sfcWind,sfcWind_attrs),
                "pr": (["time","lat","lon"],pr,pr_attrs)
                },
            coords = {
                "time":("time",unique_datenum,time_attrs),
                "lat": ("lat",lat,lat_attrs),
                "lon": ("lon",lon,lon_attrs)
                },
            attrs = global_attrs)

        export_ds = export_ds.sel(lat=slice(latlim[1],latlim[0]),
            lon=slice(lonlim[0],lonlim[1]))
        export_ds = export_ds.sortby("lat")

        # ----------------------------------------------------------------------
        # Convert calendar to "noleap"
        # ----------------------------------------------------------------------
        export_ds = xr.decode_cf(export_ds)
        export_ds = h.convert_calendar(export_ds,"noleap")
        dt_new = h.uconvert_time(export_ds["time"].values,
            "days since 1850-01-01 12:00:00")

        export_ds = export_ds.merge({"time": dt_new},overwrite_vars="time")
        export_ds["time"].attrs["units"] = "days since 1850-01-01 12:00:00"
        export_ds["time"].attrs["calendar"] = "noleap"

        # ----------------------------------------------------------------------
        # Write tasmax to netcdf file.
        # ----------------------------------------------------------------------
        export_ds["tasmax"].to_netcdf(
            os.path.join(dest_dir,
                "tasmax_era5_%d.nc" % yr[y]),
            mode="w",
            format="NETCDF4",
            engine="h5netcdf",
            encoding={
                "time": {"dtype": "int","_FillValue": None},
                "lat":  {"dtype": "single","_FillValue": None},
                "lon": {"dtype": "single","_FillValue": None},
                "tasmax": {"dtype": "single","_FillValue": -9999},
                })

        # ----------------------------------------------------------------------
        # Write hursmin to netcdf file.
        # ----------------------------------------------------------------------
        export_ds["hursmin"].to_netcdf(
            os.path.join(dest_dir,
                "hursmin_era5_%d.nc" % yr[y]),
            mode="w",
            format="NETCDF4",
            engine="h5netcdf",
            encoding={
                "time": {"dtype": "int","_FillValue": None},
                "lat":  {"dtype": "single","_FillValue": None},
                "lon": {"dtype": "single","_FillValue": None},
                "hursmin": {"dtype": "single","_FillValue": -9999},
                })

        # ----------------------------------------------------------------------
        # Write sfcWind to netcdf file.
        # ----------------------------------------------------------------------
        export_ds["sfcWind"].to_netcdf(
            os.path.join(dest_dir,
                "sfcWind_era5_%d.nc" % yr[y]),
            mode="w",
            format="NETCDF4",
            engine="h5netcdf",
            encoding={
                "time": {"dtype": "int","_FillValue": None},
                "lat":  {"dtype": "single","_FillValue": None},
                "lon": {"dtype": "single","_FillValue": None},
                "sfcWind": {"dtype": "single","_FillValue": -9999},
                })

        # ----------------------------------------------------------------------
        # Write pr to netcdf file.
        # ----------------------------------------------------------------------
        export_ds["pr"].to_netcdf(
            os.path.join(dest_dir,
                "pr_era5_%d.nc" % yr[y]),
            mode="w",
            format="NETCDF4",
            engine="h5netcdf",
            encoding={
                "time": {"dtype": "int","_FillValue": None},
                "lat":  {"dtype": "single","_FillValue": None},
                "lon": {"dtype": "single","_FillValue": None},
                "pr": {"dtype": "single","_FillValue": -9999},
                })

if verbose:
    print("Finished processing ERA5 data!\n")