#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#SBATCH --job-name=delta_2010-2015
#SBATCH --ntasks=1
#SBATCH --mem=25Gb
#SBATCH --partition=nodes
#SBATCH --time=04:30:00
#SBATCH --output=Logs/2010-2015_delta_%A.log
import sys
import os
import xarray as xr
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.colors import LogNorm
import numpy as np
import pandas as pd
matplotlib.use('agg')
sys.path.append('/users/mjr583/python_lib')
import RowPy as rp
import matplotlib.colors as mcolors
plt.style.use('seaborn-darkgrid')

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

def get_data_as_xr(rundir, year='', lev=0, version='13.1.2', variable=False):
    path=f'/mnt/lustre/users/mjr583/GC/{version}/rundirs/{rundir}'
    year = rundir[-4:]
    if variable=='OH':
        ds = xr.open_mfdataset( find_file_list(path, [f'ConcAfterChem.{year}'])[:], combine='by_coords' )
        #ds = ds[f'OHconcAfterChem'].isel(lev=lev) 
    else:
        try:
            ds = xr.open_mfdataset( find_file_list(path, [f'SpeciesConc.{year}'])[:], combine='by_coords' )
        except:
            ds = xr.open_mfdataset( find_file_list(path, [f'SpeciesConc.{year}'])[:-1], combine='by_coords' )

    lat = ds['lat']
    lon = ds['lon']
    return ds, lat, lon

def plot_max(plottable):
    ax_Max=[]
    for i in plottable:
        ax_Max.append( i.max().values )
    ax_Max = np.nanmax(ax_Max)
    return ax_Max

def surface_data(species, rundirs, versions=['13.1.2']):
    add_=[]
    for n, rundir in enumerate(rundirs):
        #print( rundir )
        ds0, lat, lon = get_data_as_xr(rundir, year='', version=versions[0], variable=species)
        ds0 = ds0.mean('time')
        ds0 = ds0.isel( lev=0 )
        if species=='OH':
            add_.append( ds0['OHconcAfterChem']  )
        elif species=='all_NIT':
            add_.append( ( ds0[f'SpeciesConc_NIT'] + ds0[f'SpeciesConc_NITs'] + ds0[f'SpeciesConc_NITD1'] + 
                           ds0[f'SpeciesConc_NITD2'] + ds0[f'SpeciesConc_NITD3'] + ds0[f'SpeciesConc_NITD4']) * 1e9 )
        elif species=='all_SO4':
            add_.append( ( ds0[f'SpeciesConc_SO4'] + ds0[f'SpeciesConc_SO4s'] + ds0[f'SpeciesConc_SO4D1'] + 
                           ds0[f'SpeciesConc_SO4D2'] + ds0[f'SpeciesConc_SO4D3'] + ds0[f'SpeciesConc_SO4D4']) * 1e9 )
        elif species=='NOx':
            add_.append( ( ds0[f'SpeciesConc_NO'] + ds0[f'SpeciesConc_NO2'] ) * 1e9 )
        elif species=='Bry':
            add_.append(( ds[f'SpeciesConc_Br'] + (ds[f'SpeciesConc_Br2']*2) + ds[f'SpeciesConc_HOBr'] \
                 + ds[f'SpeciesConc_BrO'] + ds[f'SpeciesConc_HBr'] + ds[f'SpeciesConc_BrNO2']  \
                 + ds[f'SpeciesConc_BrNO3'] + ds[f'SpeciesConc_IBr'] + ds[f'SpeciesConc_BrCl']   \
                         ) * 1e12 )
        elif species=='Cly':
            add_.append(( ds[f'SpeciesConc_Cl'] + (ds[f'SpeciesConc_Cl2']*2) + ds[f'SpeciesConc_HOCl'] \
                 + ds[f'SpeciesConc_ClO'] + ds[f'SpeciesConc_HCl'] + ds[f'SpeciesConc_ClNO2']  \
                 + ds[f'SpeciesConc_ClNO3'] + ds[f'SpeciesConc_ICl'] + ds[f'SpeciesConc_BrCl']   \
                 + ds[f'SpeciesConc_ClOO'] + ds[f'SpeciesConc_OClO'] + (ds[f'SpeciesConc_Cl2O2']*2)   \
                         ) * 1e12 )
        elif species=='Iy':
            add_.append(( ds[f'SpeciesConc_I'] + (ds[f'SpeciesConc_I2']*2) + ds[f'SpeciesConc_HOI'] \
                 + ds[f'SpeciesConc_IO'] + ds[f'SpeciesConc_OIO'] + ds[f'SpeciesConc_HI']  \
                 + ds[f'SpeciesConc_INO'] \
                 + (ds[f'SpeciesConc_I2O2']*2) + (ds[f'SpeciesConc_I2O3']*2) + (ds[f'SpeciesConc_I2O4']*2)   \
                         ) * 1e12 )
        else:
            add_.append( ds0[f'SpeciesConc_{species}'] * 1e9 )  
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

def cart_plot(species, rundirs, base, jrun, sname=''):
    fig=plt.figure(figsize=(15,4))
    all_plottable = [jrun[0] / base[0], jrun[1] / base[0], jrun[2] / base[0] ]
    all_plottable = [- (100 - (jrun[0] / base[0] * 100)), 
                     - (100 - (jrun[1] / base[0] * 100)),
                     - (100 - (jrun[2] / base[0] * 100)) ]

    labels=[ "Kasibhatla et al. / Base", "Shah et al. / Base", "Andersen et al. / Base" ]

    Max = plot_max(all_plottable)
    #print(Max)
    for n in range(len(all_plottable)):
        print( n )
        #print( all_plottable[n].min().values )
        unit = f'$\Delta$ {species}'
        unit = f'$\Delta$ {species} (%) '

        ax = fig.add_subplot(1,3,n+1, projection=ccrs.Robinson(), aspect='auto')
        
        im = all_plottable[n].plot.imshow( x='lon',y='lat', ax=ax, transform=ccrs.PlateCarree(), 
                                      vmin=1e-4, vmax=Max, norm=LogNorm(), center=0., cmap=cmap, add_colorbar=False)

        ax.coastlines()
        ax.set_title(labels[n], fontsize=16)
    plt.tight_layout() 
    plt.subplots_adjust( bottom=.08)
    cbar_ax = fig.add_axes([0.1, 0.13, 0.8, 0.04])
    cbar = fig.colorbar(im, cax=cbar_ax, orientation='horizontal')
    cbar.ax.set_xlabel(f'{unit}', size=15)
    
    plt.savefig( f'plots/{sname}.png' , dpi=300)
    plt.close()
    sys.exit()
    return

##----------------------------------------MAIN SCRIPT----------------------------------------##
def main():
    global d, site, rundirs, v
    versions = ['13.1.2']
    
    species = ['O3','OH','HNO2','NOx'][::-1]

    for s in species:
        print( s )
        rundirs  = ['base_run_2015']   
        base_ = surface_data(s, rundirs, versions)
        
        rundirs  = ['j25_run_2015','viral_run_2015','Ander22b_run_2015']
        j100_ = surface_data(s, rundirs, versions)
        
        cart_plot( s, rundirs, base_, j100_, sname=f'3panel_diff-maps_{s}_Log_pc')
        
if __name__ == "__main__":
    main()
