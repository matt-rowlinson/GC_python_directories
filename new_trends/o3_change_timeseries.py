#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#SBATCH --job-name=jscale_maps
#SBATCH --ntasks=1
#SBATCH --mem=16Gb
#SBATCH --partition=test
#SBATCH --time=00:30:00
#SBATCH --output=Logs/jscale_map_plots_%A.log
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
    year = rundir[-4:]
    if year=='e-50':
        year='2019'
    print( year )
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

def ts_plot(rundirs, base_, j100_, sname=''):
    fig, ax = plt.subplots(figsize=(12,6))
    c = ['#edf8b1','#c7e9b4','#7fcdbb','#41b6c4','#1d91c0','#225ea8','#0c2c84']
    c = ['#ece2f0','#d0d1e6','#a6bddb','#67a9cf','#3690c0','#02818a','#016450']
    c = ['#d0d1e6','#a6bddb','#67a9cf','#3690c0','#02818a','#016c59','#014636']
    c = ['#a65628','#c7e9b4','#7fcdbb','#41b6c4','#1d91c0','#225ea8','#253494','#081d58','#e41a1c','#984ea3','darkgrey']
    for n, run in enumerate(rundirs):
        print( n )
        # drop leap year day to plot all years alongside each other
        base_[n] = base_[n].sel(time=~((base_[n].time.dt.month == 2) & (base_[n].time.dt.day == 29)))
        j100_[n] = j100_[n].sel(time=~((j100_[n].time.dt.month == 2) & (j100_[n].time.dt.day == 29)))
        
        z = (j100_[n].values / base_[n].values * 100) - 100
        print( z )
        m = str(np.round( z.mean(), 2))
        if "viral" in run:
            ax.plot( base_[0]['time'], z, c=c[n], label=f'Shah et al. ({m}%)' )
        elif "j25" in run:
            ax.plot( base_[0]['time'], z, c=c[n], label=f'J25 ({m}%)' )
        elif "scale-50" in run:
            ax.plot( base_[0]['time'], z, c=c[n], label=f'J50 ({m}%)' )

        else:
            ax.plot( base_[0]['time'], z, c=c[n], label=f'{run[-4:]} ({m}%)' )
    
    import matplotlib.dates as mdates
    myFmt = mdates.DateFormatter('%b')
    ax.xaxis.set_major_formatter(myFmt)
    ax.set_ylim( bottom=0. )
    ax.set_ylabel( '$\Delta$ O$_3$ due to NIThv (%)' )

    plt.tight_layout()
    plt.subplots_adjust(right=.82)
    plt.legend(loc="center left", bbox_to_anchor=[1., .7], fancybox=True, frameon=True, facecolor='w')
    plt.savefig(f'plots/{sname}.png')
    plt.close()
    return

##----------------------------------------MAIN SCRIPT----------------------------------------##
def main():
    global d, site, rundirs, v
    versions = ['13.1.2']
    interval = [10, 10, 5, 5, 5, 5, 5]
    v='$O_3$'

    rundirs  = ['base_run_1751', 'base_run_1980', 'base_run_1990', 'base_run_2000', 'base_run_2005',
                                 'base_run_2010', 'base_run_2015', 'base_run_2020', 
                                 'base_run_2015', 'base_run_2015', 'base_run_2015']
    base_ = mean_surface_data(rundirs, versions)
    
    rundirs  = ['j100_run_1751', 'j100_run_1980', 'j100_run_1990', 'j100_run_2000', 'j100_run_2005',
                                 'j100_run_2010', 'j100_run_2015', 'j100_run_2020', 
                                 'viral_run_2015', 'j25_run_2015', 'nitrate_photol_2015_all_scale-50']
    j100_ = mean_surface_data(rundirs, versions)
    print( len( j100_) )
    print( j100_[10] )
    print( base_[10] )
    #print( j100_[-2].values )
    #print( j100_[-1].values )
    print( base_[10].shape )
    print( j100_[10].shape )

    #sys.exit()
    ts_plot(rundirs, base_, j100_, sname='TEST.delta-o3_timeseries')
        
if __name__ == "__main__":
    main()
