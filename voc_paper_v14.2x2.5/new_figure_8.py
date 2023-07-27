#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#SBATCH --job-name=map_plots
#SBATCH --ntasks=1
#SBATCH --mem=15gb
#SBATCH --partition=interactive
#SBATCH --time=00:11:30
#SBATCH --output=Logs/map_plots.log
import os
import sys
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import matplotlib
matplotlib.use('agg')
sys.path.append('/users/mjr583/python_lib')
import GC_tools as GC
import RowPy as rp
import matplotlib.colors as mcolors
from CVAO_dict import CVAO_dict as d
from mpl_toolkits.axes_grid1 import make_axes_locatable
import cartopy.crs as ccrs

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

def get_data_as_xr(rundir, year='', version='14.0.1', lev=0, variable=False):
    path=f'/mnt/lustre/users/mjr583/GC/{version}/rundirs/{rundir}'
    ds = xr.open_mfdataset( find_file_list(path, [f'SpeciesConc.{year}']), combine='by_coords' )
    if variable:
        ds = ds[f'SpeciesConc_{d[variable]["GC_name"]}']
    return ds

def plot_max_min(plottable):
    ax_Max=[]
    for i in plottable:
        Min, Max = rp.get_abs_max_min(i)
        ax_Max.append( Max.values )
    return ax_Max

def main():
    version='14.0.1'
    variables=['O3','OH']
    rundirs = ['geo','scale_all']
    rundirs = ['geo_2x25','all_2x25']

    base = get_data_as_xr('geo_2x25', version=version, year='2016', variable='O3')
    dev0 = get_data_as_xr('all_2x25', version=version, year='2016', variable='O3')

    base_p1 = base.mean(dim='time').isel(lev=0) * 1e9
    dev0_p1 = dev0.mean(dim='time').isel(lev=0) * 1e9
    p1 =  - ( 100 - ( dev0_p1 / base_p1 * 100 ))

    base_p2 = base.mean(dim='time').mean(dim='lon') * 1e9
    dev0_p2 = dev0.mean(dim='time').mean(dim='lon') * 1e9
    p2 =  - ( 100 - ( dev0_p2 / base_p2 * 100 ))

    base = get_data_as_xr('geo_2x25', version=version, year='2016', variable='OH')
    dev0 = get_data_as_xr('all_2x25', version=version, year='2016', variable='OH')

    base_p3 = base.mean(dim='time').isel(lev=0) * 1e9
    dev0_p3 = dev0.mean(dim='time').isel(lev=0) * 1e9
    p3 =  - ( 100 - ( dev0_p3 / base_p3 * 100 ))

    base_p4 = base.mean(dim='time').mean(dim='lon') * 1e9
    dev0_p4 = dev0.mean(dim='time').mean(dim='lon') * 1e9
    p4 =  - ( 100 - ( dev0_p4 / base_p4 * 100 ))

    p = pd.read_csv( '/users/mjr583/GC/info_files/GC_72_vertical_levels.csv')['Pressure (hPa)']
    panel_labels=['(a)','(b)','(c)','(d)']
    fig=plt.figure(figsize=(10,8))
    
    ## (a) O3 surface diff
    ax=fig.add_subplot( 2,2,1, projection=ccrs.EqualEarth(), aspect='auto')
    im1 = p1.plot.imshow( x='lon',y='lat', ax=ax, transform=ccrs.PlateCarree(), 
                    center=0., cmap='bwr', 
                    cbar_kwargs={'orientation':'horizontal','pad':0.01,'shrink':0.6, 'aspect':40,'label':r'$\Delta$ O$_3$ (%)'})
    ax.coastlines()
    ax.set_title('(a) Surface [O$_3$] difference') 

    ## (c) O3 zonal diff
    ax=fig.add_subplot( 2,2,3)
    im2 = ax.pcolormesh( p2['lat'], p[:32], p2[:32,:], cmap='bwr', vmin=-.8, vmax=.8)#norm=norm )
    ax.invert_yaxis()
    ax.set_xlabel( 'Latitude' ) ; ax.set_ylabel( 'Pressure (hPa)' )
    divider = make_axes_locatable(ax)
    cax = divider.append_axes('bottom', size='5%', pad=0.5)
    cbar = fig.colorbar(im2, cax=cax, orientation='horizontal')
    cbar.ax.set_xlabel(r'$\Delta$ O$_3$ (%)',fontsize=10)
    ax.set_title('(c) Zonal mean [O$_3$] difference') 

    ## (c) OH surface diff
    ax=fig.add_subplot( 2,2,2, projection=ccrs.EqualEarth(), aspect='auto')
    im3 = p3.plot.imshow( x='lon',y='lat', ax=ax, transform=ccrs.PlateCarree(), 
                    center=0., cmap='bwr', 
                    cbar_kwargs={'orientation':'horizontal','pad':0.01,'shrink':0.6, 'aspect':40,'label':r'$\Delta$ OH (%)'})
    ax.coastlines()
    ax.set_title('(b) Surface [OH] difference') 

    ## (d) OH zonal diff
    ax=fig.add_subplot( 2,2,4)
    im4= ax.pcolormesh( p4['lat'], p[:32], p4[:32,:], cmap='bwr',vmin=-1.6, vmax=1.6)#norm=norm )
    ax.invert_yaxis()
    ax.set_xlabel( 'Latitude' ) ; ax.set_ylabel( 'Pressure (hPa)' )
    divider = make_axes_locatable(ax)
    cax = divider.append_axes('bottom', size='5%', pad=0.5)
    cbar = fig.colorbar(im4, cax=cax, orientation='horizontal')
    cbar.ax.set_xlabel(r'$\Delta$ OH (%)',fontsize=10)#, rotation=270)
    ax.set_title('(d) Zonal mean [OH] difference') 
   
    plt.tight_layout()
    plt.savefig( f"plots/figure_8.png", dpi=300 )
    plt.close()

if __name__=="__main__":
    main()
