#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#SBATCH --job-name=cart_test
#SBATCH --ntasks=1
#SBATCH --mem=16Gb
#SBATCH --partition=test
#SBATCH --time=00:30:00
#SBATCH --output=Logs/cart_test_%A.log
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

current_dir = os.path.dirname(__file__)
rgb_WhGrYlRd = np.genfromtxt('/users/mjr583/python_lib/colormaps/WhGrYlRd.txt',
                                     delimiter=' ')
WhGrYlRd = mcolors.ListedColormap(rgb_WhGrYlRd/255.0)
cmap=WhGrYlRd

def get_arguments():
    parser = argparse.ArgumentParser(description="Parse arguments to pass to GC processing scripts")
    parser.add_argument("-r", "--rundir", type=str, 
                        help='Name of desired GC rundir')
    parser.add_argument("-v", "--var", type=list,
                        default=["O3"],
                        help="Name of GC variable")
    parser.add_argument("-V", "--version", type=str,
                        default='13.1.2',
                        help="Version of GEOS-Chem")
    args=parser.parse_args()
    return args.rundir, args.var, args.version

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
    if variable=='OH':
        ds = xr.open_mfdataset( find_file_list(path, [f'ConcAfterChem.{year}'])[:], combine='by_coords' )
        ds = ds[f'OHconcAfterChem'].isel(lev=lev) #* 1e9
    else:
        try:
            ds = xr.open_mfdataset( find_file_list(path, [f'SpeciesConc.{year}'])[:], combine='by_coords' )
        except:
            ds = xr.open_mfdataset( find_file_list(path, [f'SpeciesConc.{year}'])[:-1], combine='by_coords' )

    lat = ds['lat']
    lon = ds['lon']
    return ds, lat, lon

def plot_max_min(plottable):
    ax_Max=[]
    for i in plottable:
        Min, Max = rp.get_abs_max_min(i)
        ax_Max.append( Max.values )
    return ax_Max

def plot(f, axes, plottable, lat, lon, c='', labels=False,
        cbar_label=[],sname="4_panel_plot.png",
        panel_labels=['(a)','(b)']):
    plottable=plottable[2:]
    X, Y = np.meshgrid( lon, lat )
    ax_Max = plot_max_min(plottable)
    for n in range(len(plottable)):
        print( n )
        m = rp.get_basemap(freq=60, ax=axes[n])
        im = m.pcolormesh(X, Y, plottable[n], cmap='bwr', vmax=ax_Max[n], vmin=-(ax_Max[n]),)
        axes[n].set_title(f'{panel_labels[n]} {labels[n]}', fontsize=14 )
        
        divider = make_axes_locatable(axes[n])
        cax = divider.append_axes('bottom', size='5%', pad=0.25)
        cbar = f.colorbar(im, cax=cax, orientation='horizontal')
        cbar.ax.set_xlabel(f'{cbar_label[n]}',fontsize=12)

    plt.tight_layout()
    plt.savefig( f"plots/{sname}" )
    plt.close()
    return


##----------------------------------------MAIN SCRIPT----------------------------------------##
def main():
    global d, site, rundirs, v
    rundir, variables, version = get_arguments()
    rundirs  = ['base_run_1980', 'base_run_1990','base_run_2000', 'base_run_2005',
                                 'base_run_2010', 'base_run_2015', 'base_run_2020']
    versions = ['13.1.2']

    variables=['O3']
    units = ['ppbv']
    for nn, v in enumerate(variables):
        print( v )
        add_=[]
        for n, rundir in enumerate(rundirs):
            print( rundir )
            ds0, lat, lon = get_data_as_xr(rundir, year='', version=versions[0])
            ds0 = ds0.mean('time')
            ds0 = ds0.isel( lev=0 )
            add_.append( ds0 )

        print( len( add_ ) ) 
        
        fig=plt.figure(figsize=(8,9))
        first=True
        for n in range(1,7):

            ds = (( add_[n] - add_[0] ) * 1e9)
            ds = ds['SpeciesConc_NIT'] + ds['SpeciesConc_NITs'] + ds['SpeciesConc_NITD1'] \
               + ds['SpeciesConc_NITD2'] + ds['SpeciesConc_NITD3'] + ds['SpeciesConc_NITD4']
            ax = fig.add_subplot(3,2,n, projection=ccrs.Robinson(), aspect='auto')
            im = ds.plot.imshow( x='lon',y='lat', ax=ax, transform=ccrs.PlateCarree(),
                        cbar_kwargs={'orientation':'horizontal','shrink':0.6, 'aspect':40,'label':f' {v} ppbv'})
            ax.coastlines()
            #ax.title( rundirs[n].replace("base_run_","") )

        plt.savefig( 'plots/NIT_ccrs_TEST.png' )
        plt.close()
        
if __name__ == "__main__":
    main()
