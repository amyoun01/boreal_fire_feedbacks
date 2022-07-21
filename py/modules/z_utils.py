"""
Set of helper functions for various tasks in boreal_fire_feedbacks
project.

"""

import os
import glob
import datetime as dt
import xarray as xr
import numpy as np
from calendar import isleap
from cftime import DatetimeNoLeap

def get_global_attrs(Dataset,attrs_list):

    for i in range(0,len(attrs_list)):

        if i == 0:

            global_attrs = {attrs_list[i]: Dataset.attrs[attrs_list[i]]}

        else:

            global_attrs[attrs_list[i]] = Dataset.attrs[attrs_list[i]]

    return global_attrs

def uconvert(ds):

    """
    For use in CFFDRS calculations, convert to the following units for GCM output:

    tasmax: convert K to C
    sfcWind: convert m/s to km/h
    pr: convert kg m-2 s-2 to mm/d 
    """

    varname = list(ds.keys())
    
    if (varname == 'tasmax'):

        ds = ds - 273.15
        # units = 'degC'

    elif (varname == 'sfcWind'):

        ds = 3.6 * ds
        # units = 'km/h'

    elif (varname == 'pr'):

        ds = ds * 60 * 60 * 24
        units = 'mm/d'

    elif (varname == 'hursmin'):

        da = ds[varname]
        da = xr.where(da > 100.0,99.99,da)
        ds[varname] = da
        units = '%'
        
    return ds

def uconvert_time(dt_arr,time_units):

    unit_str_parts = time_units.split()
    unit = unit_str_parts[0]

    if unit == 'hours':

        unit = 'H'

    elif unit == 'days':

        unit = 'D'

    if len(unit_str_parts) == 3:

        base_date = unit_str_parts[2]
        base_date = dt.datetime.strptime(base_date,'%Y-%m-%d')
        base_date = DatetimeNoLeap(base_date.year,
                                   base_date.month,
                                   base_date.day)

    elif len(unit_str_parts) == 4:

        base_date = ' '.join([unit_str_parts[2],unit_str_parts[3]])
        base_date = dt.datetime.strptime(base_date,'%Y-%m-%d %H:%M:%S')
        base_date = DatetimeNoLeap(base_date.year,
                                   base_date.month,
                                   base_date.day,
                                   base_date.hour,
                                   base_date.minute,
                                   base_date.second)

    td = np.array(dt_arr - base_date,dtype='timedelta64[%s]' % unit)
    td = td.astype('int')

    return td

def extract_str_parts(filelist,extension='',sep='_'):

    n_files = len(filelist)
    n_filename_parts = len(os.path.split(filelist[0])[1].split(sep))

    filename_parts = np.zeros((n_files,n_filename_parts),dtype=object)

    for i in range(n_files):

        file_i = os.path.split(filelist[i])[1].split('.' + extension)[0]
        filename_parts_i = file_i.split(sep)
        
        filename_parts[i,] = filename_parts_i
    
    return filename_parts

def overlap(x,lims,format='%Y-%m-%d',index_return = False):

    if (type(x) is list):

        x = np.array(x)

    if x.size == 2:

        if (type(x) is str):

            x = [dt.datetime.strptime(x[0],format), \
                 dt.datetime.strptime(x[1],format)]

        bool_1 = (x[0] < lims[0]) & (x[1] > lims[1])
        bool_2 = (x[0] > lims[1]) & (x[0] < lims[1])
        bool_3 = (x[1] > lims[1]) & (x[1] < lims[1])

        overlap_bool = bool_1 | bool_2 | bool_3

    else:

        overlap_bool = np.zeros((x.shape[0],),dtype=bool)

        for i in range(0,len(x)):

            if (type(x[i,0]) is str):
                
                x_a = [dt.datetime.strptime(x[i,0],format), \
                       dt.datetime.strptime(x[i,1],format)]

            else:

                x_a = x[i,:]

            bool_1 = (x_a[0] < lims[0]) & (x_a[1] > lims[1])
            bool_2 = (x_a[0] > lims[0]) & (x_a[0] < lims[1])
            bool_3 = (x_a[1] > lims[0]) & (x_a[1] < lims[1])

            overlap_bool[i] = bool_1 | bool_2 | bool_3

    if index_return:

        overlap_bool = np.flatnonzero(overlap_bool)

    return overlap_bool


def convert_360_to_noleap(Dataset):

    target_calendar = 'noleap'

    varnames = list(Dataset.keys())
    _, n, m = list(Dataset.sizes.values())

    missing_doy_365 = np.array([36,109,183,255,329]) + 1

        
    Dataset = Dataset.convert_calendar(target_calendar,
                                    align_on='year')

    min_yr = np.min(Dataset.time.dt.year)
    max_yr = np.max(Dataset.time.dt.year)

    ref_dt = xr.cftime_range(start = '%d-01-01' % min_yr,
                            end    = '%d-12-31' % max_yr,
                            calendar=target_calendar)

    for v in varnames:
        
        id = np.isin(ref_dt.dayofyear,missing_doy_365,invert=True)

        new_arr = np.nan * np.ones((len(ref_dt),n,m),dtype='single')
        new_arr[id,...] = Dataset[v].copy()

        id_v = np.flatnonzero(~id)

        new_arr[id_v,...] = (new_arr[id_v-1,...] + new_arr[id_v+1,...]) / 2
        
        new_arr = new_arr[id_v,...]
        ref_dt = ref_dt[id_v]

        new_ds = xr.Dataset({v: (['time','lat','lon'],new_arr)},
                            coords = {'time': ('time',ref_dt),
                                    'lat': ('lat',Dataset['lat'].values),
                                    'lon': ('lon',Dataset['lon'].values)})        

        if v is varnames[0]:

            export_ds = Dataset.merge(new_ds) 

        else:
         
            export_ds = export_ds.merge(new_ds)

    return export_ds            


def convert_360_to_standard(Dataset):

    target_calendar = 'standard'

    varnames = list(Dataset.keys())
    _, n, m = list(Dataset.sizes.values())

    missing_doy_365 = np.array([36,109,183,255,329]) + 1
    missing_doy_366 = np.array([31,91,153,213,275,335]) + 1
        
    Dataset = Dataset.convert_calendar(target_calendar,
                                    align_on='year')

    min_yr = np.min(Dataset.time.dt.year)
    max_yr = np.max(Dataset.time.dt.year)

    ref_dt = xr.cftime_range(start = '%d-01-01' % min_yr,
                            end    = '%d-12-31' % max_yr,
                            calendar=target_calendar)

    yr = ref_dt.year
    leap = np.array([isleap(x) for x in yr],dtype='bool')

    for v in varnames:
        
        id0 = ~leap & np.isin(ref_dt.dayofyear,missing_doy_365,invert=True)
        id1 = leap & np.isin(ref_dt.dayofyear,missing_doy_366,invert=True)

        id = id0 | id1

        new_arr = np.nan * np.ones((len(ref_dt),n,m),dtype='single')
        new_arr[id,...] = Dataset[v].copy()

        id_v = np.flatnonzero(~id)

        new_arr[id_v,...] = (new_arr[id_v - 1,...] + new_arr[id_v+1,...]) / 2
        
        new_arr = new_arr[id_v,...]
        ref_dt = ref_dt[id_v]

        new_ds = xr.Dataset({v: (['time','lat','lon'],new_arr)},
                            coords = {'time': ('time',ref_dt),
                                    'lat': ('lat',Dataset['lat'].values),
                                    'lon': ('lon',Dataset['lon'].values)})        

        if v is varnames[0]:

            export_ds = Dataset.merge(new_ds) 

        else:
         
            export_ds = export_ds.merge(new_ds)

    return export_ds

def convert_calendar(Dataset,target_calendar):

    calendar_types = {'standard': 'standard',
                    'gregorian': 'standard',
                    'proleptic_gregorian': 'standard',
                    'noleap': 'noleap',
                    '365_day': 'noleap',
                    '360_day': '360_day'}

    source_calendar = Dataset.time.dt.calendar
    
    source_calendar = calendar_types[source_calendar]
    target_calendar = calendar_types[target_calendar]

    try: 
        source_calendar is target_calendar

    except ValueError:
        raise(ValueError('Source and target calendar are the same,\
            no need for conversion.'))

    if (source_calendar == '360_day') and (target_calendar == 'noleap'):

        export_ds = convert_360_to_noleap(Dataset)

    elif (source_calendar == '360_day') and (target_calendar == 'standard'):

        export_ds = convert_360_to_standard(Dataset)

    elif (source_calendar == 'standard') and (target_calendar == 'noleap'):

        export_ds = Dataset.convert_calendar('noleap')

    return export_ds

def vp_calc(tmp):

    # tmp value needs to be degK, not degC
    # Calculate vp in Pa. Formula = Sonntag 1990
    vp = 611.2 * np.exp((17.62 * tmp) / (243.12 + tmp))
    vp = 0.001 * vp # Convert from Pa to kPa

    return vp    

def concat_annual_to_dataset(wdir,yr,glob_pattern="*",recursive=False,extension="nc",sep="_"):

    glob_pattern = os.path.join(wdir,glob_pattern + extension)
    filelist = glob.glob(glob_pattern,recursive=recursive)
    filelist.sort()
    
    filename_parts = extract_str_parts(filelist,extension=extension,sep=sep)

    yr_col = [i for i in range(0,len(filename_parts[0,:])) if filename_parts[0,i].isnumeric()]
    yr_col = yr_col[0]

    file_yr = filename_parts[:,yr_col].astype('int')
    
    yr_id = np.flatnonzero(np.isin(file_yr,yr))
    
    files_to_import = list(np.array(filelist)[yr_id])

    ds = xr.open_mfdataset(files_to_import,
                        combine='nested',
                        parallel=True)

    ds = ds.load()

    return ds
