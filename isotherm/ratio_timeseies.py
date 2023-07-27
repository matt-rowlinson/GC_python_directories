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
    rundirs  = ['dev_new_base','dev_onlyNIThv','dev_onlyNO-HO2','dev_both']
    archive_dirs  = ['dev.Base.Test_1month','dev.NIThv.Test_1month',
                     'dev.NO-HO2.Test_1month','dev.NIThv_HO2-NO.Test_1month']

    labels=['Base', 'dev.NIThv','dev.NO-HO2rxn','dev.NIThv+NO-HO2rxn']
    cs=['#1b9e77','#e7298a','#d95f02','#7570b3','#e7298a','#66a61e','#e6ab02'][::-1]
    
    x=-24.9 ; y=16.9 ; z=0
    x = site['longitude'] ; y=site['latitude'] ; z = 0
    
    f, ax = plt.subplots( figsize=(10,6) )
    style=['-','--','-','--']

    for n, rundir in enumerate(rundirs):
        print( rundir )
        ds  = get_data_as_xr(rundir, year='2018012', archive=True, archive_dir=archive_dirs[n])
        Met = get_data_as_xr(rundir, collection='StateMet', year='', archive=True, archive_dir=archive_dirs[n])

        tn16 = site_data( ds['SpeciesConc_TN016'], lon=x, lat=y, lev=z)
        tn17 = site_data( ds['SpeciesConc_TN017'], lon=x, lat=y, lev=z)
        T    = site_data( Met['Met_T'], lon=x, lat=y, lev=z)
        p    = site_data( Met['Met_PMIDDRY'], lon=x, lat=y, lev=z)
        H    = site_data( Met['Met_SPHU'], lon=x, lat=y, lev=z)
        
        print( T.values.mean(), p.values.mean(), H.values.mean() )
        r = tn17 / tn16 * 100
    
        ax.plot( ds.time, r, label=labels[n], ls=style[n])
    ax.text(.1, .9, "TEXT", fontsize=50) 
    plt.ylabel( 'Prod_NOnHO2_Channel2 (%)' )
    plt.legend()
    plt.savefig( 'plots/ratio_CVO.png' )
    plt.close()

if __name__ == "__main__":
    main()
