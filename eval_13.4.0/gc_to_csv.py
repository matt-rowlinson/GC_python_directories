import os
import xarray as xr
import glob
from netCDF4 import Dataset
import numpy as np
import pandas as pd
import sys
import datetime
import re
import matplotlib.pyplot as plt

sys.path.append('/users/mjr583/python_lib')
import GC_tools as GC
import RowPy as rp

def find_file_list(path, substrs):
    file_list =[]
    for root, directory, files in os.walk(path):
        for f in files:
            for s in substrs:
                if s in f:
                    file_list.append(os.path.join(root, f))
    file_list.sort()
    return file_list

def get_data_as_xr(rundir, year='', variable=False, version='13.1.2'):
    path=f'/mnt/lustre/users/mjr583/GC/{version}/rundirs/{rundir}'
    ds = xr.open_mfdataset( find_file_list(path, [f'SpeciesConc.{year}'])[:-1], combine='by_coords' )
    return ds

def site_data(ds, lat=32.35, lon=-64.77, lev=0):
    x = rp.find_nearest(ds.lon, lon)-1
    y = rp.find_nearest(ds.lat, lat)+1
    if type(lev)==int:
        data = ds.isel( lon=x, lat=y, lev=lev)
    else:
        data = ds.isel( lon=x, lat=y ).isel( lev=slice(lev)).mean(dim='lev', keep_attrs=True)
    return data

def find_model_output_for_site(ds, rundirs, versions, v, years='2017'):
    timeseries=[] ; add_lowes=[]
    for n, rundir in enumerate(rundirs):
        ds0 = get_data_as_xr(rundir, year=years, variable=v, version=versions[n])

        ds0 = site_data( ds0, ds )
        ds0 = pd.DataFrame({site['save_name']:ds0.values}, index=ds0['time'].values)#.resample('M').mean()
        
        if site['save_name']=='cvao':
            ds0=ds0['2006-10-01':]
        #ds0=ds0['2017':'2018']
        timeseries.append( ds0 )
    return timeseries

def main():
    variables=['O3','HNO2','HNO3','NO','NO2','NIT','NITs','NITD1','NITD2','NITD3','NITD4']
    units = [ 1e9,1e12, 1e12, 1e12, 1e12, 1e12,1e12,1e12,1e12,1e12,1e12,1e12,1e12]
    rundir = 'TEST'
    ds = get_data_as_xr(rundir, year='', version='13.4.0')
    holder=[]
    for n, v in enumerate(variables):
        print( v )
        v0 = ds[f'SpeciesConc_{v}']
        v0 = site_data( v0, lat=32.35, lon=-64.77 )
        holder.append( v0.values*units[n] )

    df = pd.DataFrame( holder ).T#, columns=variables, index=ds['time'].values)
    df.columns=variables
    df.index=ds['time'].values
    df.index.name='datetime'
    print( df )
    df=df['2018':]
    df.to_csv('./v13-4-0_bermuda.csv' )
    
if __name__=="__main__":
    main()
