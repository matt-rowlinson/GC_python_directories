#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#SBATCH --job-name=jscale_time
#SBATCH --ntasks=1
#SBATCH --mem=18Gb
#SBATCH --partition=interactive
#SBATCH --time=00:30:00
#SBATCH --output=Logs/jscale_time_plots_%A.log
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
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
import RowPy as rp
import argparse
from scipy.optimize import curve_fit
from sklearn.metrics import mean_squared_error, r2_score
plt.style.use('seaborn-darkgrid')

def get_cosweights(lat):
    cosweights = np.cos(np.deg2rad(lat))
    cosweights = np.swapaxes((np.tile(cosweights,72)).reshape(72,46),0,1)
    return cosweights

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
    path=f'/mnt/lustre/users/mjr583/GC/{version}/rundirs/{rundir}/OutputDir/'
    ds = xr.open_mfdataset( find_file_list(path, [f'SpeciesConc.{year}']), combine='by_coords' )
    cosw = get_cosweights( ds.lat )
    ds = ds[f'SpeciesConc_{variable}'] * 1e12 * cosw
    return ds

def site_data(ds, lat=16.9, lon=-24.9, lev=0):
    x = rp.find_nearest(ds.lon, lon)-1
    y = rp.find_nearest(ds.lat, lat)+1
    if type(lev)==int:
        data = ds.isel( lon=x, lat=y, lev=lev)
    elif type(lev) == bool:
        data = ds.isel( lon=x, lat=y )
    else:
        data = ds.isel( lon=x, lat=y ).isel( lev=slice(lev)).mean(dim='lev', keep_attrs=True)
    return data

def find_model_output_for_site(rundirs, versions, v, years='2017', site='CVO', diurnal=False):
    timeseries=[] ; add_lowes=[]
    for n, rundir in enumerate(rundirs):
        print( rundir )
        ds0 = get_data_as_xr(rundir, year=years, variable=v, version=versions[0])
        ds0 = site_data( ds0, lon=site['longitude'], lat=site['latitude'], lev=False )
        ds0 = ds0.where(ds0['time.hour'] == 13, drop=True)
        ds_ = ds0.mean(dim='time')
        d25 = ds0.quantile(.25, dim='time')
        d75 = ds0.quantile(.75, dim='time')
        
        ds0 = pd.DataFrame({site['save_name']:ds_.values}, index=ds_['lev'].values)
        ds0["q25"] = d25
        ds0["q75"] = d75

        timeseries.append( ds0 )
    return timeseries


##----------------------------------------MAIN SCRIPT----------------------------------------##
def main():
    
    rundirs  = ['base_run_2015','j25_run_2015','viral_run_2015','Ander22b_run_2015']
    variable='HNO2'
    years='201502'

    lev=pd.read_csv('/users/mjr583/GC/info_files/GC_72_vertical_levels.csv')['Altitude (km)']
    times = find_model_output_for_site(rundirs, ['13.1.2'], variable, years=years, site=sites['CVO'], diurnal=False)

    f, ax = plt.subplots(1,1, figsize=(5,5))
    cs = ['#377eb8','#984ea3','#ff7f00','#e41a1c']
    labels = ['Base','Kasibhatla et al. 2018','Shah et al. 2022','Andersen et al. 2022']
    # Delete when j25 is done
    #cs = ['#377eb8','#ff7f00','#e41a1c']
    #labels = ['Base','Shah et al. 2022','Andersen et al. 2022']
    l=32
    for n in range(len(rundirs)):
        ax.plot( times[n].cvao.values[:l], lev[:l], c=cs[n], label=labels[n] )
        ax.fill_betweenx( lev[:l], times[n].q25.values[:l], times[n].q75.values[:l],
                                                    color=cs[n], alpha=.2, zorder=3 )
    ax.set_ylabel( 'Altitude (km)')
    ax.set_xlabel( 'HONO (ppt)')
    ax.set_ylim( bottom=1, top=5 )
    ax.legend()
    plt.savefig(f'plots/HNO2_vertical_Feb.png', dpi=200)
    plt.close()

if __name__ == "__main__":
    main()
