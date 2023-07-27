#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#SBATCH --job-name=jscale_time
#SBATCH --ntasks=1
#SBATCH --mem=18Gb
#SBATCH --partition=interactive
#SBATCH --time=00:30:00
#SBATCH --output=Logs/timeseries_%A.log
#import warnings
#warnings.simplefilter(action='ignore', category=FutureWarning)
import sys
import os
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import matplotlib
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
import pandas as pd
matplotlib.use('agg')
sys.path.append('/users/mjr583/python_lib')
from sites_dicts import GAW_dict as sites
from CVAO_dict import CVAO_dict as d
import RowPy as rp
import argparse
plt.style.use('seaborn-darkgrid')

def find_file_list(path, substrs):
    file_list =[]
    for root, directory, files in os.walk(path+'/OutputDir/'):
        for f in files:
            for s in substrs:
                if s in f:
                    file_list.append(os.path.join(root, f))
    file_list.sort()
    return file_list

def get_data_as_xr(rundir, year='', collection='SpeciesConc', archive=False, archive_dir='',variable=False, Chem=False, version='13.1.2'):
    if archive:
        path=f'/mnt/lustre/users/mjr583/GC/{version}/rundirs/{rundir}/{archive_dir}'
    else:
        path=f'/mnt/lustre/users/mjr583/GC/{version}/rundirs/{rundir}'

    if variable=='OH' or variable=='HO2':
        collection='ConcAfterChem'
        Chem=True
    try:
        if year=='201001':
            ds = xr.open_mfdataset( find_file_list(path, [f'{collection}.{year}']), combine='by_coords' )
        else:
            ds = xr.open_mfdataset( find_file_list(path, [f'{collection}.{year}']), combine='by_coords' )
    except:
        try:
            ds = xr.open_mfdataset( find_file_list(path, [f'{collection}.{year}']), combine='by_coords' )
        except:
            ds = xr.open_mfdataset( find_file_list(path, [f'{collection}.{year}'])[:-1], combine='by_coords' )
    return ds

def site_data(ds, lat=16.9, lon=-24.9, lev=0):
    x = rp.find_nearest(ds.lon, lon)
    y = rp.find_nearest(ds.lat, lat)
    if type(lev)==int:
        data = ds.isel( lon=x, lat=y, lev=lev)
    else:
        data = ds.isel( lon=x, lat=y ).isel( lev=slice(lev)).mean(dim='lev', keep_attrs=True)
    return data

def map_plot(species, r, rundir, n, levels=None):
    plt.figure(figsize=(10,6))
    cbar_kwargs = {'orientation':'horizontal', 'label':f'[{species}]: dev / base','spacing': 'proportional'}#, 'ticks':levels[::10]}
    ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=0))
    im = r.plot.imshow(x='lon', y='lat', ax=ax, cmap='bwr', levels=levels, 
                              transform=ccrs.PlateCarree(), cbar_kwargs=cbar_kwargs
                              )
    ax.coastlines()
    plt.title(f'v13 {rundir} / Base') 
    plt.savefig( f'plots/{rundir}_{species}_pc.png' )
    plt.close()
    return

##----------------------------------------MAIN SCRIPT----------------------------------------##
def main():
    rundirs  = ['dev_new_base','dev_both','dev_lim_HO2nNO']
    #archive_dirs  = ['dev.Base.Test_1month','dev.NIThv.Test_1month',
    #                 'dev.NO-HO2.Test_1month','dev.NIThv_HO2-NO.Test_1month']
    labels=['Base', 'dev.NoLimit','dev.Limited']#,'dev.NIThv+NO-HO2rxn']
    levels=None 
    for n, rundir in enumerate(rundirs[1:]):
        print( rundir )
        n=n+1
        ds0 = get_data_as_xr(rundirs[0], year='201801')#, archive=True, archive_dir=archive_dirs[0])
        ds  = get_data_as_xr(rundir, year='201801')#, archive=True, archive_dir=archive_dirs[n])

        v0 = ds0['SpeciesConc_O3'].mean('time').isel(lev=0) * 1e9
        v  = ds['SpeciesConc_O3'].mean('time').isel(lev=0) * 1e9
        print( v0.values.mean(), v.values.mean() ) 
        #levels = [-3., -2.5, -2., -1.5, -1., -.5, 0., .5, 1., 1.5, 2., 2.5, 3.]
        levels = np.arange( -2.8, 2.85, .05 )
        map_plot( 'O3', v / v0 , rundir, n)#, levels=levels)
        
        levels = np.arange( -22.5, 22.75, .25 )
        v0 = ( ds0['SpeciesConc_NO'] +ds0['SpeciesConc_NO2'] ).mean('time').isel(lev=0) *1e12
        v  = ( ds['SpeciesConc_NO'] +ds['SpeciesConc_NO2'] ).mean('time').isel(lev=0) * 1e12
        map_plot( 'NOx', v / v0 , rundir, n, levels=levels)

if __name__ == "__main__":
    main()
