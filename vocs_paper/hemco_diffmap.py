#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#SBATCH --job-name=viral_comp
#SBATCH --ntasks=1
#SBATCH --mem=48Gb
#SBATCH --partition=himem
#SBATCH --time=03:30:00
#SBATCH --output=Logs/viral_%A.log
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
    #year = rundir[-4:]
    ds = xr.open_mfdataset( find_file_list(path, [f'HEMCO_diagnostics.{year}'])[:-1], combine='by_coords' )

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

def surface_data(species, rundirs, versions=['13.1.2'], year='2017'):
    add_=[]
    for n, rundir in enumerate(rundirs):
        #print( rundir, versions, species )
        ds0, lat, lon = get_data_as_xr(rundir, year=year, version=versions[0], variable=species)
        ds0 = ds0.mean('time')
        ds0 = ds0.isel( lev=0 )
        ds0 = ds0[species]
        add_.append( ds0 )
    return add_


def make_plottable( ds_, plot_type='abs'):
    plottable_=[]
    for n in range( 1, len(ds_) ):
        if plot_type=='abs':
            plottable_.append( ( ds_[n] - ds_[n-1] ) )
        elif plot_type=='pc':
            plottable_.append((- ( 100 - ( ds_[n] / ds_[n-1] * 100 )))) 
    return plottable_

def the_diff( a, b ):
    c_=[]
    for n, aa in enumerate( a ):
        c_.append( b[n] - a[n] )
    return c_

def cart_plot(species, rundirs, ds, plot_type='abs', sname=''):
    fig=plt.figure(figsize=(15,4))
    first=True

    all_plottable = [ ds[0], ds[1] ]#, (ds[1] / ds[0] * 100)-100 ]
    for n in range(0,3):
        #print( n, plot_type )
        if plot_type=='abs':
            labels=["Base", "Geo-only", "Geo - Base"]
            all_plottable.append( ds[1] - ds[0] ) 
            unit = 'kg m-2 s-1'
            center=0.
        elif plot_type=='pc':
            labels=["Base", "Geo-only", "Geo / Base"]
            all_plottable.append( ds[1] / ds[0] ) #np.divide( ds[1], ds[0], where=ds[0]!=0) )
            center=1.
            if first:
                unit = '%'
        levels=None
        ax = fig.add_subplot(1,3,n+1, projection=ccrs.Robinson(), aspect='auto')
        if n < 2:
            Max= np.nanmax([all_plottable[0], all_plottable[1] ] )
            im = all_plottable[n].plot.imshow( x='lon',y='lat', vmax=Max,  ax=ax, transform=ccrs.PlateCarree(), 
                   add_colorbar=False)
        else:
            im0 = all_plottable[n].plot.imshow( x='lon',y='lat', ax=ax, transform=ccrs.PlateCarree(),
                    cmap='bwr', center=center, add_colorbar=False) 

        ax.coastlines()
        ax.set_title(labels[n], fontsize=16)
    
    plt.tight_layout()
    plt.subplots_adjust( bottom=.1) 
    cbar_ax = fig.add_axes([0.02, 0.15, 0.6, 0.03])
    cbar = fig.colorbar(im, cax=cbar_ax, orientation='horizontal')
    cbar.ax.set_xlabel(f'{species} (kg m-2 s-1)', size=15)

    cbar_ax = fig.add_axes([0.68, 0.15, 0.3, 0.03])
    cbar = fig.colorbar(im0, cax=cbar_ax, orientation='horizontal')
    cbar.ax.set_xlabel(f'$\Delta$ {species} ({unit})', size=15)
    
    plt.savefig( f'plots/{sname}.png' , dpi=300)
    plt.close()
    return

##----------------------------------------MAIN SCRIPT----------------------------------------##
def main():
    global d, site, rundirs, v
    versions = ['14.0.1']
    
    species = [#'EmisNO_Anthro','EmisNO_Ship','EmisNO_BioBurn','EmisNO_Aircraft','EmisNO_Soil','EmisNO_Anthro',
              #'EmisNH3_Ship','EmisNH3_Seabirds','EmisNH3_Natural','EmisNH3_BioBurn', 'EmisNH3_Anthro','EmisNH3_Anthro']
               #'EmisC2H6_Anthro','EmisC3H8_Anthro',
               #'EmisBENZ_Anthro',
               'EmisEOH_Anthro','EmisPRPE_Anthro','EmisTOLU_Anthro']
               #'EmisCO_Anthro','EmisCO_Anthro',
               #'EmisNO_Anthro','EmisNO_Anthro']
    year='2017'
    for s in species:
        print( s )
        rundirs  = ['base_run','geo']   
        ds_ = surface_data(s, rundirs, versions, year=year)
        
        cart_plot( s, rundirs, ds_, plot_type='abs', sname=f'v14.geo-base.HEMCO_{s}_ab')
        cart_plot( s, rundirs, ds_, plot_type='pc' , sname=f'v14.geo-base.HEMCO_{s}_pc')

if __name__ == "__main__":
    main()
