#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#SBATCH --job-name=o3_rf_calc
#SBATCH --ntasks=1
#SBATCH --mem=16Gb
#SBATCH --partition=test
#SBATCH --time=00:30:00
#SBATCH --output=Logs/o3_rf_calc_%A.log
import sys
import os
import xarray as xr
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
import pandas as pd
matplotlib.use('agg')
sys.path.append('/users/mjr583/python_lib')
from CVAO_dict import CVAO_dict as d
import RowPy as rp
import argparse
import matplotlib.colors as mcolors
plt.style.use('seaborn-darkgrid')
import yaml
 
current_dir = os.path.dirname(__file__)
rgb_WhGrYlRd = np.genfromtxt('/users/mjr583/python_lib/colormaps/WhGrYlRd.txt',
                                     delimiter=' ')
WhGrYlRd = mcolors.ListedColormap(rgb_WhGrYlRd/255.0)
cmap=WhGrYlRd

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

def get_data_as_xr(rundir, year='', lev=0, collection='SpeciesConc', version='13.1.2', variable=False):
    path=f'/mnt/lustre/users/mjr583/GC/{version}/rundirs/{rundir}'
    year = rundir[-4:]
    print( path, f'{collection}.{year}' ) 
    if variable=='OH':
        ds = xr.open_mfdataset( find_file_list(path, [f'ConcAfterChem.{year}'])[:], combine='by_coords' )
        ds = ds[f'OHconcAfterChem'].isel(lev=lev) 
    else:
        try:
            ds = xr.open_mfdataset( find_file_list(path, [f'{collection}.{year}'])[:], combine='by_coords' )
        except:
            ds = xr.open_mfdataset( find_file_list(path, [f'{collection}.{year}'])[:-1], combine='by_coords' )

    lat = ds['lat']
    lon = ds['lon']
    return ds, lat, lon

def data_3d(rundir, collection='SpeciesConc',variable='O3', version='13.1.2'):
    #add_=[]
    #for n, rundir in enumerate(rundirs):
        #print( rundir, variable )
    ds0, lat, lon = get_data_as_xr(rundir, collection=collection, variable=variable, year='', version=version)
    if variable=='AREA':
        return ds0[f'{variable}'].mean(dim='time') 
    elif collection=='StateMet':
        return ds0[f'Met_{variable}'].mean(dim='time') 
    elif collection=="HEMCO_diagnostics":
        return ds0[f'{variable}'].mean(dim='time') 
    else:
        return  ds0[f'SpeciesConc_{variable}'].mean(dim='time') # * 1e9 )
    return add_

def cart_plot(rundirs, base, sname=''):
    fig=plt.figure(figsize=(6,4))
    ax = fig.add_subplot(1,1,1, projection=ccrs.EqualEarth(), aspect='auto')
    im = base.plot.imshow( x='lon',y='lat',  ax=ax, vmin=0, vmax=1.2 , transform=ccrs.PlateCarree(),
            cbar_kwargs = {'orientation':'horizontal', 'extend':'neither','label':'W m$^{-2}$','spacing': 'proportional'}) 
    ax.coastlines()
    
    plt.tight_layout()
    plt.savefig( f'plots/{sname}.png' , dpi=300)
    plt.close()
    return

##----------------------------------------MAIN SCRIPT----------------------------------------##
def main():
    mm_air  = (.78*(2.*14.)+.22*(2.*16.))
    AVO     = 6.0221412e23

    #rundirs  = ['j25-SS_run_1751','j25-SS_run_2015']#,'j100_run_1751','j100_run_2000']#
    rundirs  = ['base_run_1751','base_run_2015']
    #rundirs  = ['Ander22b_run_1751','Ander22b_run_2015']

    DU_=[]
    for rundir in rundirs:
        print( rundir )
        vmr     = data_3d(rundir, version='13.1.2')
        AIRDEN  = data_3d(rundir, variable='AIRDEN'  , collection="StateMet")
        AIRVOL  = data_3d(rundir, variable='AIRVOL'  , collection="StateMet")
        AREA    = data_3d(rundir, variable='AREA'    , collection="StateMet")
        trop    = data_3d(rundir, variable='TropP'   , collection="StateMet")
        pmid    = data_3d(rundir, variable='PMID'    , collection="StateMet")
        MASK =  ( pmid > trop )

        vmr = vmr.where(MASK)  # Get tropospheric ozone only
        molecs = (AIRDEN / 1e6) / (mm_air / 1e3) * AVO * (AIRVOL) # ( kg/cm3 ) / kg/mol * molecules/cm3 * cm3
        DU = vmr * molecs               # molecules = mol/mol * molecules/gridbox
        DU = DU.sum(dim='lev')      # Sum over vertical levels to get column O3
        DU = DU / AREA              # Divide by area to get molecules/cm2
        DU = DU / 2.69e20           # Conversion to DU
        DU = DU * 1e6               ## Correcting earlier AIRVOL conversion ( leading to overflow )
        DU_.append( DU )
        print( DU.mean().values )
    
    DU1 = DU_[1] #- DU_[0]
    DU2 = DU_[0] #- DU_[2] 
    DU = DU1 - DU2
    NRF = DU * 42 * 1e-3                        # get RF in mW m-2
    weighted_NRF = np.sum( NRF * AREA) / np.sum( AREA )  # weighted by area
    print( weighted_NRF.mean().values)
    print( NRF.min().values, NRF.max().values )

    cart_plot( rundirs, NRF, sname=f'{rundirs[0][:4]}_PD-PI_NRF_2010') 

        
if __name__ == "__main__":
    main()
