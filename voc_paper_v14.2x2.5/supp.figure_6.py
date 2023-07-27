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
    return file_list

def get_data_as_xr(rundir, year='', lev=0, version='14.0.1', variable=False):
    path=f'/mnt/lustre/users/mjr583/GC/{version}/rundirs/{rundir}'
    ds = xr.open_mfdataset( find_file_list(path, [f'SpeciesConc.{year}']), combine='by_coords' )
    if variable:
        ds = ds[f'SpeciesConc_{d[variable]["GC_name"]}'] * float(d[variable]['scale'])
    lat = ds['lat']
    lon = ds['lon']
    return ds, lat, lon


def plot_max_min(plottable):
    ax_Max=[]
    for i in plottable:
        Min, Max = rp.get_abs_max_min(i)
        ax_Max.append( Max.values )
    return ax_Max

####------------------------------------------------------------------------

def main():
    version='14.0.1'
    variables=['C2H6','C2H6','EOH','CH2O','PRPE','MEK','ALK4','BENZ','TOLU','XYLE']#,'
    rundirs = ['geo_2x25','all_2x25']
    p = pd.read_csv( '/users/mjr583/GC/info_files/GC_72_vertical_levels.csv')['Pressure (hPa)']

    add_plottable=[] 
    for v in variables:
        print( v )
        dss=[]
        for n, rundir in enumerate(rundirs):
            ds0, lat, lon = get_data_as_xr(rundir, year='2016', version=version, variable=v)
            ds0 = ds0.mean(dim='time').mean(dim='lon')
            print( ds0.shape )
            dss.append( ds0 )
        #add_plottable.append( - ( 100 - ( dss[1] / dss[0] * 100 ))) 
        
        add_plottable.append( ( dss[1] - dss[0] ) / dss[0] * 100 ) 


    #ax_Max = plot_max_min(add_plottable)
    #X, Y = np.meshgrid( lon, lat )
    panel_labels=['(a)','(b)','(c)','(d)','(e)','(f)','(g)','(h)','(i)','(j)']
    fig = plt.figure(figsize=(12,12))
    for n in range(10):
        print( n )
        ax = fig.add_subplot(4,3,n+1)
        im2 = ax.pcolormesh( ds0['lat'], p[:32], add_plottable[n][:32,:], cmap='bwr', vmin=-80., vmax=80.)#norm=norm )
        ax.invert_yaxis()
        ax.set_xlabel( 'Latitude' ) ; ax.set_ylabel( 'Pressure (hPa)' )
        divider = make_axes_locatable(ax)
        cax = divider.append_axes('bottom', size='5%', pad=0.5)
        cbar = fig.colorbar(im2, cax=cax, orientation='horizontal')
        cbar.ax.set_xlabel(r'$\Delta$'+f' {variables[n]} (%)',fontsize=10)

        ax.set_title(f'{panel_labels[n]} {variables[n]}', fontsize=12 )
    
    plt.tight_layout()
    plt.savefig( f"plots/Supp.figure_6.png", dpi=200 )
    plt.close()
    sys.exit()

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

if __name__=="__main__":
    main()
