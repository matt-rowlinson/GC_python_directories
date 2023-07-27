#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#SBATCH --job-name=delta_2010-2015
#SBATCH --ntasks=1
#SBATCH --mem=55Gb
#SBATCH --partition=himem
#SBATCH --time=05:30:00
#SBATCH --output=Logs/2010-2015_delta_%A.log
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

def surface_data(species, rundirs, versions=['13.1.2']):
    add_=[]
    for n, rundir in enumerate(rundirs):
        #print( rundir )
        ds0, lat, lon = get_data_as_xr(rundir, year='', version=versions[0], variable=species)
        ds0 = ds0.mean('time')
        ds0 = ds0.isel( lev=0 )
        if species=='OH':
            ds0 = ds0['OHconcAfterChem']  
        elif species=='all_NIT':
            ds0 = ( ( ds0[f'SpeciesConc_NIT'] + ds0[f'SpeciesConc_NITs'] + ds0[f'SpeciesConc_NITD1'] + 
                           ds0[f'SpeciesConc_NITD2'] + ds0[f'SpeciesConc_NITD3'] + ds0[f'SpeciesConc_NITD4']) * 1e9 )
        elif species=='all_SO4':
            ds0 = ( ( ds0[f'SpeciesConc_SO4'] + ds0[f'SpeciesConc_SO4s'] + ds0[f'SpeciesConc_SO4D1'] + 
                           ds0[f'SpeciesConc_SO4D2'] + ds0[f'SpeciesConc_SO4D3'] + ds0[f'SpeciesConc_SO4D4']) * 1e9 )
        elif species=='NOx':
            ds0 = ( ( ds0[f'SpeciesConc_NO'] + ds0[f'SpeciesConc_NO2'] ) * 1e9 )
        elif species=='Bry':
            ds0 = (( ds[f'SpeciesConc_Br'] + (ds[f'SpeciesConc_Br2']*2) + ds[f'SpeciesConc_HOBr'] \
                 + ds[f'SpeciesConc_BrO'] + ds[f'SpeciesConc_HBr'] + ds[f'SpeciesConc_BrNO2']  \
                 + ds[f'SpeciesConc_BrNO3'] + ds[f'SpeciesConc_IBr'] + ds[f'SpeciesConc_BrCl']   \
                         ) * 1e12 )
        elif species=='Cly':
            ds0 = (( ds[f'SpeciesConc_Cl'] + (ds[f'SpeciesConc_Cl2']*2) + ds[f'SpeciesConc_HOCl'] \
                 + ds[f'SpeciesConc_ClO'] + ds[f'SpeciesConc_HCl'] + ds[f'SpeciesConc_ClNO2']  \
                 + ds[f'SpeciesConc_ClNO3'] + ds[f'SpeciesConc_ICl'] + ds[f'SpeciesConc_BrCl']   \
                 + ds[f'SpeciesConc_ClOO'] + ds[f'SpeciesConc_OClO'] + (ds[f'SpeciesConc_Cl2O2']*2)   \
                         ) * 1e12 )
        elif species=='Iy':
            ds0 = (( ds[f'SpeciesConc_I'] + (ds[f'SpeciesConc_I2']*2) + ds[f'SpeciesConc_HOI'] \
                 + ds[f'SpeciesConc_IO'] + ds[f'SpeciesConc_OIO'] + ds[f'SpeciesConc_HI']  \
                 + ds[f'SpeciesConc_INO'] \
                 + (ds[f'SpeciesConc_I2O2']*2) + (ds[f'SpeciesConc_I2O3']*2) + (ds[f'SpeciesConc_I2O4']*2)   \
                         ) * 1e12 )
        else:
            ds0 = ds0[f'SpeciesConc_{species}'] * 1e9

        ds0 = ds0.where( 30. < ds0.lat , drop=True )
        ds0 = ds0.where( ds0.lat < 70. , drop=True )
        ds0 = ds0.where( -16. < ds0.lon , drop=True )
        ds0 = ds0.where( ds0.lon < 45. , drop=True )
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

def cart_plot(species, rundirs, base, jrun, plot_type='abs', sname=''):
    fig=plt.figure(figsize=(15,10))
    first=True
    base_ = make_plottable( base, plot_type=plot_type )
    jrun_ = make_plottable( jrun, plot_type=plot_type )
    diff_of_diff = the_diff( base_, jrun_ )
    all_plottable = [base[0], jrun[0], jrun[0], base[1], jrun[1], jrun[1], base_[0], jrun_[0], diff_of_diff[0]]
    labels=["Base: 2010", "j100: 2010", "x",
            "Base: 2015", "j100: 2015", "x",
            "Base: 2015-2010", "j100: 2015-2010", "Diff: j100-Base"]
    for n in range(0,9):
        print( n )
        if n==2 or n==5:
            continue
        if plot_type=='abs':
            unit = f'{species} (ppbv)'
        elif plot_type=='pc':
            if first:
                sname+='_PC'
                unit = f' {species} (%)'
        levels=None
        ax = fig.add_subplot(3,3,n+1, projection=ccrs.Robinson(), aspect='auto')
        if n <= 5:
            Max = np.max( [base[0], jrun[0], base[1], jrun[1] ] )
            im0 = all_plottable[n].plot.imshow( x='lon',y='lat', ax=ax, transform=ccrs.PlateCarree(), 
                    vmax=Max, cmap=cmap, add_colorbar=False)
            first=False
        elif n==8:
            im2 = all_plottable[n].plot.imshow( x='lon',y='lat', ax=ax, transform=ccrs.PlateCarree(), 
                    levels=levels, center=0., cmap='bwr', add_colorbar=False) 
        else:
            im = all_plottable[n].plot.imshow( x='lon',y='lat', ax=ax, transform=ccrs.PlateCarree(), 
                    levels=levels, center=0., cmap='bwr', add_colorbar=False)

        ax.coastlines()
        ax.set_extent([-16, 45, 30, 70])
        ax.set_title(labels[n], fontsize=16)
    plt.tight_layout() 
    plt.subplots_adjust( bottom=.1)#, hspace=.1) 
    cbar_ax = fig.add_axes([0.67, 0.45, 0.02, 0.5])
    cbar = fig.colorbar(im0, cax=cbar_ax, orientation='vertical')
    cbar.ax.set_ylabel(f'{species} ppbv', size=15)

    cbar_ax = fig.add_axes([0.03, 0.05, 0.6, 0.02])
    cbar = fig.colorbar(im, cax=cbar_ax, orientation='horizontal')
    cbar.ax.set_xlabel(f'{unit}', size=15)

    cbar_ax = fig.add_axes([0.67, 0.05, 0.3, 0.02])
    cbar = fig.colorbar(im2, cax=cbar_ax, orientation='horizontal')
    cbar.ax.set_xlabel(f'$\Delta${unit}', size=15)

    #plt.tight_layout()
    plt.savefig( f'plots/{sname}.png' , dpi=300)
    plt.close()

    return

##----------------------------------------MAIN SCRIPT----------------------------------------##
def main():
    global d, site, rundirs, v
    versions = ['13.1.2']
    
    species = ['O3','HNO2','CO','ALK4','NO','NO2','C2H6',
               'NIT','NITs','NITD1','NITD2','NITD3','NITD4',
               'SO4','SO4','SO4s','SO4D1','SO4D2','SO4D3','SO4D4',
               'OH',
               'all_NIT','all_SO4',
               'Cly','Bry','Iy']#[::-1]
    for s in species:
        print( s )
        rundirs  = ['base_run_2010', 'base_run_2015']   
        base_ = surface_data(s, rundirs, versions)
        
        rundirs  = ['j100_run_2010', 'j100_run_2015']
        j100_ = surface_data(s, rundirs, versions)
        
        cart_plot( s, rundirs, base_, j100_, plot_type='pc' , sname=f'EUall5_maps_2010-2015_{s}')
        cart_plot( s, rundirs, base_, j100_, plot_type='abs', sname=f'EUall5_maps_2010-2015_{s}')
        
if __name__ == "__main__":
    main()
