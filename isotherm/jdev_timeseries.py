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

##----------------------------------------MAIN SCRIPT----------------------------------------##
def main():
    site = sites['CVO']
    rundirs = ['dev_new_base','dev_onlyNIThv','dev_onlyNO-HO2','dev_both']
    labels  = ['Base', 'dev.NIThv','dev.NO-HO2rxn','dev.NIThv+NO-HO2rxn']
    cs=['#1b9e77','#e7298a','#d95f02','#7570b3','#e7298a','#66a61e','#e6ab02'][::-1]
    
    x=-24.9 ; y=16.9 ; z=0
    x = site['longitude'] ; y=site['latitude'] ; z = 0
    style=['-','--','-','--']*3
    a=[]

    for n, rundir in enumerate(rundirs):
        print( rundir )
        SpeciesConc = get_data_as_xr(rundir, year='201801')
        StateMet    = get_data_as_xr(rundir, collection='StateMet', year='201801')
        Jvalues     = get_data_as_xr(rundir, collection='JValues', year='201801')

        NO3 = site_data(  SpeciesConc['SpeciesConc_NIT']   + SpeciesConc['SpeciesConc_NITs'] 
                        + SpeciesConc['SpeciesConc_NITD1'] + SpeciesConc['SpeciesConc_NITD2']
                        + SpeciesConc['SpeciesConc_NITD3'] + SpeciesConc['SpeciesConc_NITD4'],                        
                        lon=x, lat=y, lev=z) * 1e12
        HNO2 = site_data(  SpeciesConc['SpeciesConc_HNO2'], lon=x, lat=y, lev=z) * 1e12

        jHNO3 = site_data( Jvalues['Jval_HNO3'], lon=x, lat=y, lev=z)
        jHNO2 = site_data( Jvalues['Jval_HNO2'], lon=x, lat=y, lev=z)
        a.append( jHNO3 * NO3  )
        a.append( jHNO2 * HNO2 )

    f, ax = plt.subplots( figsize=(10,6) )
 
    ax.plot( a[0]['time'], a[0], label='jHNO3 * NO3', ls=style[0])
    ax.plot( a[1]['time'], a[1], label='jHONO * HONO', ls=style[1])

    plt.ylabel( f"1e12 * {SpeciesConc['SpeciesConc_NIT'].units} * {Jvalues['Jval_HNO3'].units}" )
    plt.legend()
    plt.savefig( 'plots/jdev_TEST.png' )
    plt.close()

if __name__ == "__main__":
    main()
