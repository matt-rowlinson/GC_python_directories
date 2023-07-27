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
#import cartopy.crs as ccrs
import matplotlib
#from matplotlib.colors import LogNorm
#from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
import pandas as pd
matplotlib.use('agg')
#sys.path.append('/users/mjr583/python_lib')
#from sites_dicts import GAW_dict as sites
#from CVAO_dict import CVAO_dict as d
#import RowPy as rp
#import argparse
import yaml
from matplotlib.backends.backend_pdf import PdfPages
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
            ds = xr.open_mfdataset( find_file_list(path, [f'{collection}.{year}']), combine='by_coords' )
    return ds


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
    spec_db = yaml.load(open("./species_database.yml"), Loader=yaml.FullLoader)

    Ycore_spec=['NIT']#,'CO','NO','NO2']
    rundirs  = ['base_run_1980','base_run_1990','base_run_2000','base_run_2005',
                'base_run_2010','base_run_2015','base_run_2020']
    
    years=[]
    for n, rundir in enumerate(rundirs):
        print( rundir )
        year = rundir[-4:]
        ds  = get_data_as_xr(rundir, year=year)
        CaC = get_data_as_xr(rundir, collection='ConcAfterChem', year=year)
        Met = get_data_as_xr(rundir, collection='StateMet', year=year)
        trop = Met['Met_TropP'].mean( dim='time' )
        AirMass = Met['Met_AD'].mean( dim='time' )
        pmid_press = Met['Met_PMID'].mean( dim='time' )
        MASK =  ( pmid_press > trop )
        
        a=[]
        for nn, spec in enumerate(Ycore_spec):
            var = ds[f'SpeciesConc_{spec}'] + ds[f'SpeciesConc_{spec}s'] + ds[f'SpeciesConc_{spec}D1'] \
                    + ds[f'SpeciesConc_{spec}D2'] + ds[f'SpeciesConc_{spec}D3'] + ds[f'SpeciesConc_{spec}D4']
            var = var.mean(dim='time')

            var = var.where(MASK)
            RMM_air = (.78*(2.*14.)+.22*(2.*16.))
            conversion_factor = (AirMass*1E3 / RMM_air)
            var = var * conversion_factor
            var = var * float(spec_db[Ycore_spec[0]]['MW_g'] * 1e-3 ) / 1E9

            a.append( var.sum().values )
            years.append( year )
            print( var.sum().values ) 
            #a.append( ds[f'SpeciesConc_{spec}'].mean(dim='time').isel(lev=0).mean().values * 1e9 )
    print( a )
    a =  np.array( a )
    print( a )


    fig, ax =plt.subplots(figsize=(5,5))
    ax.bar( years, a )
    plt.savefig(f'plots/TEST.NIT-burden_timeseries.png')

if __name__ == "__main__":
    main()
