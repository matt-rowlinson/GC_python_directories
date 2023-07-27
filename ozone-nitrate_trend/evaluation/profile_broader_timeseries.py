#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#SBATCH --job-name=North_Atl
#SBATCH --ntasks=1
#SBATCH --mem=2Gb
#SBATCH --partition=interactive
#SBATCH --time=00:15:00
#SBATCH --output=Logs/NorthAtl_timeseries_%A.log
import sys
import os
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
import pandas as pd
matplotlib.use('agg')
sys.path.append('/users/mjr583/python_lib')
from sites_dicts import GAW_dict as sites
from CVAO_dict import CVAO_dict as d
import GC_tools as GC
import RowPy as rp
import argparse
from scipy.optimize import curve_fit
from sklearn.metrics import mean_squared_error, r2_score
plt.style.use('seaborn-darkgrid')

def find_file_list(path, substrs):
    file_list =[]
    for root, directory, files in os.walk(path):
        for f in files:
            for s in substrs:
                if s in f:
                    file_list.append(os.path.join(root, f))
    file_list.sort()
    return file_list

def get_data_as_xr(rundir, year='', years=False):
    path=f'/mnt/lustre/users/mjr583/GC/13.1.2/rundirs/{rundir}'
    if years:
        f_list=[]
        for y in years:
            f_list.append( find_file_list(path, [f'SpeciesConc.{y}']) )
            ds = xr.open_mfdataset( f_list, combine='by_coords' )
    else:
        ds = xr.open_mfdataset( find_file_list(path, [f'SpeciesConc.{year}']), combine='by_coords' )
    return ds

def site_data(data, ds, lat=16.9, lon=-24.9, lev=0, box=False, llon=False, ulon=False, llat=False, ulat=False):

    if box:
        ll_x = rp.find_nearest(ds.lon, llon ).values
        ur_x = rp.find_nearest(ds.lon, ulon ).values
        ll_y = rp.find_nearest(ds.lat, llat ).values
        ur_y = rp.find_nearest(ds.lat, ulat ).values

        data = data.isel( lev=0)
        data = data.isel( lon=slice(ll_x,ur_x)).mean(dim='lon', keep_attrs=True)
        data = data.isel( lat=slice(ll_y,ur_y)).mean(dim='lat', keep_attrs=True)
    else:
        data = data.isel( lon=lon).mean(dim='lon', keep_attrs=True)
        data = data.isel( lat=lat).mean(dim='lat', keep_attrs=True)
        data = data.isel( lev=slice(lev)).mean(dim='lev', keep_attrs=True)
    return data

def load_observations( site, variable ):
    if site['site_name'] == 'CVAO':
        df = pd.read_csv( site['filepath'], index_col=0, dtype={"Airmass":str, "New_Airmass":str})
        df = df[d[variable]['merge_name']]
    else:
        df = pd.read_csv( site['filepath'],index_col=0)
        df = df[variable]
    df.index = pd.to_datetime( df.index, format="%Y-%m-%d %H:%M:%S")
    return df

box_dict = {
        "North_Atlantic" : {
                            "llon" : -52,
                            'ulon' : -10,
                            'llat' : 30,
                            'ulat' : 60,
                            "c"    : 'purple' },
        "CV_North" : {
                            "llon" : -31.9,
                            'ulon' : -17.9,
                            'llat' : 16.9,
                            'ulat' : 24.9,
                            "c"    : 'red' },
        "CV_Broad" : {
                            "llon" : -34.9,
                            'ulon' : -31.9,
                            'llat' : 11.9,
                            'ulat' : 38.9,
                            "c"    : 'y' },
        "CV_NE" : {
                            "llon" : -24.9,
                            'ulon' : -17.9,
                            'llat' : 16.9,
                            'ulat' : 36.9,
                            "c"    : 'g' },
        "MH_East" : {   
                            "llon" : -34.9,
                            'ulon' : -9.9,
                            'llat' : 43.3,
                            'ulat' : 63.3,
                            "c"    : 'skyblue' },
        "MH_East_small" : {
                            "llon" : -24.9,
                            'ulon' : -9.9,
                            'llat' : 48.3,
                            'ulat' : 58.3,
                            "c"    : 'darkblue' },
        }


##----------------------------------------MAIN SCRIPT----------------------------------------##
def main():
    from CVAO_dict import CVAO_dict as d
    ds = get_data_as_xr(f'trimmed_nitrate_photol_control', year='198')
    box_name='CV_NE'
    
    llon=int(box_dict[box_name]["llon"]) ; ulon=int(box_dict[box_name]["ulon"])
    llat=int(box_dict[box_name]["llat"]) ; ulat=int(box_dict[box_name]["ulat"])

    variables=['O3','NIT','NITs','SO2','NO','NO2','CO','SO4','HNO3','NH3']
    for v in variables:
        data = ds[f'SpeciesConc_{v}']
        data = site_data( data, ds, box=True, llon=llon, ulon=ulon, llat=llat, ulat=ulat ) * float(d[v]['scale'])
        date = data.resample(time='1MS').mean(dim='time')
        data = pd.DataFrame({'XValue':np.arange(len(data.values)), box_name : data.values}, index=data['time'].values)
        data = data.resample( 'M' ).mean()
        
        f, ax = plt.subplots(1,1,figsize=(10,4))
        ax.plot( data.index, data[box_name], c=box_dict[box_name]["c"], label='v13.1.2', alpha=.5,zorder=2)
        #data = data.resample( 'M' ).mean()
        
        sys.path.append( '/users/mjr583/cvao' )
        from lowess_smoother import loess
        eval_DF = loess("XValue", box_name, data = data, alpha=.35)
        ax.plot( data.index, eval_DF['g'][1:], c=box_dict[box_name]["c"] )
        
        plt.legend(loc=0)
        #plt.savefig( f"plots/{box_name}_{v}.png" )
        plt.close()
        
if __name__ == "__main__":
    main()
