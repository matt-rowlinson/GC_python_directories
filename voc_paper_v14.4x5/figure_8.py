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

current_dir = os.path.dirname(__file__)
rgb_WhGrYlRd = np.genfromtxt('/users/mjr583/python_lib/colormaps/WhGrYlRd.txt',
                                     delimiter=' ')
WhGrYlRd = mcolors.ListedColormap(rgb_WhGrYlRd/255.0)
cmap=WhGrYlRd

def find_file_list(path, substrs):
    file_list =[]
    for root, directory, files in os.walk(path):
        for f in files:
            for s in substrs:
                if s in f:
                    file_list.append(os.path.join(root, f))
    file_list.sort()
    print( file_list )
    return file_list

def get_data_as_xr(rundir, year='', version='14.0.1', lev=0, variable=False):
    path=f'/mnt/lustre/users/mjr583/GC/{version}/rundirs/{rundir}'
    ds = xr.open_mfdataset( find_file_list(path, [f'SpeciesConc.{year}']), combine='by_coords' )
    if variable:
        ds = ds[f'SpeciesConc_{d[variable]["GC_name"]}'].isel(lev=lev) * float(d[variable]['scale'])
    lat = ds['lat']
    lon = ds['lon']
    return ds, lat, lon


def plot_max_min(plottable):
    ax_Max=[]
    for i in plottable:
        Min, Max = rp.get_abs_max_min(i)
        ax_Max.append( Max.values )
    return ax_Max

def main():
    version='13.1.2'
    variables=['O3','OH']
    rundirs = ['geo','scale_all']
    rundirs = ['ceds_only','asia-meic_scale']
    add_plottable=[] 
    for v in variables:
        print( v )
        dss=[]
        for n, rundir in enumerate(rundirs):
            ds0, lat, lon = get_data_as_xr(rundir, version=version, year='2017', lev=0, variable=v)
            ds0 = np.mean( ds0, 0 )
            dss.append( ds0 )
        add_plottable.append( - ( 100 - ( dss[1] / dss[0] * 100 ))) 

    ax_Max = plot_max_min(add_plottable)
    X, Y = np.meshgrid( lon, lat )
    panel_labels=['(a)','(b)']
    fig=plt.figure(figsize=(10,4))
    for n in range(2):
        print( n )
        ax=fig.add_subplot( 1,2,n+1, projection=ccrs.EqualEarth(), aspect='auto')
        print( 'mean:', add_plottable[n].values.max() ) 
        im = add_plottable[n].plot.imshow( x='lon',y='lat', ax=ax, transform=ccrs.PlateCarree(), 
                    center=0., cmap='bwr', 
                    cbar_kwargs={'orientation':'horizontal','shrink':0.6, 'aspect':40,'label':'%'})
        ax.coastlines()

        ax.set_title(f'{panel_labels[n]} {d[variables[n]]["longname"]}', fontsize=12 )
        
    plt.tight_layout()
    plt.savefig( f"plots/OLD.figure_8.png", dpi=200 )
    plt.close()

if __name__=="__main__":
    main()
