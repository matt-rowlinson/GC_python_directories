#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#SBATCH --job-name=o3bars
#SBATCH --ntasks=1
#SBATCH --mem=156Gb
#SBATCH --partition=himem
#SBATCH --time=00:30:00
#SBATCH --output=Logs/o3_bars_%A.log
import warnings
warnings.simplefilter(action='ignore', category=RuntimeWarning)
import sys
import os
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('agg')
plt.style.use('seaborn-darkgrid')
import numpy as np
import pandas as pd

def find_file_list(path, substrs):
    file_list =[]
    for root, directory, files in os.walk(path):
        for f in files:
            for s in substrs:
                if s in f:
                    file_list.append(os.path.join(root, f))
    file_list.sort()
    return file_list

def get_data_as_xr(rundir, year='', lev=0, version='13.1.2', variable=False):
    path=f'/mnt/lustre/users/mjr583/GC/{version}/rundirs/{rundir}'
    if variable=='OH':
        ds = xr.open_mfdataset( find_file_list(path, [f'ConcAfterChem.{year}'])[:], combine='by_coords' )
        ds = ds[f'OHconcAfterChem'].isel(lev=lev) 
    else:
        try:
            ds = xr.open_mfdataset( find_file_list(path, [f'SpeciesConc.{year}'])[:], combine='by_coords' )
        except:
            ds = xr.open_mfdataset( find_file_list(path, [f'SpeciesConc.{year}'])[:-1], combine='by_coords' )

    lat = ds['lat']
    lon = ds['lon']
    return ds, lat, lon

def mean_surface_data(rundirs, versions=['13.1.2']):
    add_=[]
    for n, rundir in enumerate(rundirs):
        print( rundir )
        ds0, lat, lon = get_data_as_xr(rundir, year='', version=versions[0])
        ds0 = ds0.mean('lon').mean('lat')
        ds0 = ds0.isel( lev=0, drop=True )
        add_.append( ds0['SpeciesConc_O3'] * 1e9 )
    return add_

def ts_plot(rundirs, base_, sname=''):
    fig, ax = plt.subplots(figsize=(12,6))
    c = ['#edf8b1','#c7e9b4','#7fcdbb','#41b6c4','#1d91c0','#225ea8','#0c2c84']
    c = ['#ece2f0','#d0d1e6','#a6bddb','#67a9cf','#3690c0','#02818a','#016450']
    c = ['#d0d1e6','#a6bddb','#67a9cf','#3690c0','#02818a','#016c59','#014636']
    for n, run in enumerate(rundirs):
        # drop leap year day to plot all years alongside each other
        base_[n] = base_[n].sel(time=~((base_[n].time.dt.month == 2) & (base_[n].time.dt.day == 29)))
        
        z = j100_[n] #/ base_[n] * 100) - 100
        m = str(np.round( z.mean().values, 2))
        ax.plot( base_[0]['time'], z, c=c[n], label=f'{run[-4:]} ({m}%)' )

    import matplotlib.dates as mdates
    myFmt = mdates.DateFormatter('%b')
    ax.xaxis.set_major_formatter(myFmt)

    ax.set_ylabel( '$O_3$ (ppbv)' )
    plt.legend()
    plt.savefig(f'plots/TEST.{sname}.png')
    plt.close()
    return

##----------------------------------------MAIN SCRIPT----------------------------------------##
def main():
    global d, site, rundirs, v
    versions = ['13.1.2']
    interval = [10, 10, 5, 5, 5, 5, 5]
    v='$O_3$'

    rundirs  = ['base_run_1750','base_run_1980', 'base_run_1990', 'base_run_2000', 'base_run_2005',
                                 'base_run_2010', 'base_run_2015', 'base_run_2020']   
    rundirs = ['nitrate_photol_2015_all_scale-25']
    base_ = mean_surface_data(rundirs, versions)

    first=True
    for ds in base_:
        if first:
            df = pd.DataFrame({"O3" : ds}, index=ds['time'])
            first=False
        else:
            df2 = pd.DataFrame({"O3" : ds}, index=ds['time'])
            df = pd.concat([df, df2]) 
    savedf =  df.resample( 'Y' ).mean()
    print( savedf )
    savedf.to_csv('csv/2015_j25_mean_surface_o3.csv')
    
        
if __name__ == "__main__":
    main()
