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
    print( file_list )
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

def site_data(ds, lat=16.9, lon=-24.9, lev=None):
    x = rp.find_nearest(ds.lon, lon)
    y = rp.find_nearest(ds.lat, lat)
    if type(lev)==int:
        data = ds.isel( lon=x, lat=y, lev=lev)
    elif lev==None:
        data = ds.isel( lon=x, lat=y )
    else:
        data = ds.isel( lon=x, lat=y ).isel( lev=slice(lev)).mean(dim='lev', keep_attrs=True)
    return data

##----------------------------------------MAIN SCRIPT----------------------------------------##
def main():
    site = sites['TUD']
    rundirs = ['new_base_4x5','andersen_4x5']#,'dev_onlyNIThv','dev_onlyNO-HO2','dev_both']
    labels  = ['Capped','Uncapped']#', 'dev.NIThv','dev.NO-HO2rxn','dev.NIThv+NO-HO2rxn']
    cs=['#1b9e77','#e7298a','#d95f02','#7570b3','#e7298a','#66a61e','#e6ab02'][::-1]
    
    x=-24.9 ; y=16.9 ; z=0
    x = site['longitude'] ; y=site['latitude'] ; z = 0
    style=['-','--','-','--']*3
    a=[]

    for n, rundir in enumerate(rundirs):
        print( rundir )
        SpeciesConc = get_data_as_xr(rundir, year='201401')
        StateMet    = get_data_as_xr(rundir, collection='StateMet', year='201401')
        Jvalues     = get_data_as_xr(rundir, collection='JValues', year='201401')

        NO3 = site_data(  SpeciesConc['SpeciesConc_NIT']   + SpeciesConc['SpeciesConc_NITs'] 
                        + SpeciesConc['SpeciesConc_NITD1'] + SpeciesConc['SpeciesConc_NITD2']
                        + SpeciesConc['SpeciesConc_NITD3'] + SpeciesConc['SpeciesConc_NITD4'],                        
                        lon=x, lat=y, lev=z) * 1e12
        HNO2 = site_data(  SpeciesConc['SpeciesConc_HNO2'], lon=x, lat=y, lev=z) * 1e12

        jHNO3 = site_data( Jvalues['Jval_HNO3'], lon=x, lat=y )
        jNIT = site_data( Jvalues['Jval_NIT'], lon=x, lat=y )
        jNITs = site_data( Jvalues['Jval_NITs'], lon=x, lat=y )
        #print( jNIT.shape )
        #sys.exit()

        #print( (jNIT / jHNO3).values )
        print( jNIT.shape)
        jNIT = jNIT.mean(dim='time')
        jHNO3 = jHNO3.mean(dim='time')

        print( jNIT.shape )
        ## Read alts from info_file and plot vertical against that (avg time dim)
        lev=pd.read_csv('/users/mjr583/GC/info_files/GC_72_vertical_levels.csv')['Altitude (km)'] #* 1e3
        aa = (jNIT / jHNO3).values
        print( aa )
        print( aa.shape )
        a.append( aa )#(jNIT / jHNO3).values  )
        #a.append( jHNO2 * HNO2 )

    f, ax = plt.subplots( figsize=(6,6) )
    print( lev ) 
    ax.plot( a[0][:35], lev[:35], label=labels[0])#'J-scale at CVAO')
    ax.plot( a[1][:35], lev[:35], label=labels[1])#'J-scale at CVAO')
    
    ax.set_ylabel( 'Altitude (km)')
    ax.set_xlabel( 'J-scale ')
    ax.set_ylim(bottom=0. )
    ax.set_xlim(left=0. )

    plt.legend()
    plt.suptitle(f'{site["save_name"]} J-values')
    plt.savefig( f'plots/J-value_{site["save_name"]}.png' )
    plt.close()

if __name__ == "__main__":
    main()
