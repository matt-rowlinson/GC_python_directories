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
    year = rundir[-4:]
    if variable=='OH':
        ds = xr.open_mfdataset( find_file_list(path, [f'ConcAfterChem.{year}'])[:], combine='by_coords' )
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
        print( rundir )
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
            ds0 = (( ds0[f'SpeciesConc_Br'] + (ds0[f'SpeciesConc_Br2']*2) + ds0[f'SpeciesConc_HOBr'] \
                 + ds0[f'SpeciesConc_BrO'] + ds0[f'SpeciesConc_HBr'] + ds0[f'SpeciesConc_BrNO2']  \
                 + ds0[f'SpeciesConc_BrNO3'] + ds0[f'SpeciesConc_IBr'] + ds0[f'SpeciesConc_BrCl']   \
                         ) * 1e12 )
        elif species=='Cly':
            ds0 = (( ds0[f'SpeciesConc_Cl'] + (ds0[f'SpeciesConc_Cl2']*2) + ds0[f'SpeciesConc_HOCl'] \
                 + ds0[f'SpeciesConc_ClO'] + ds0[f'SpeciesConc_HCl'] + ds0[f'SpeciesConc_ClNO2']  \
                 + ds0[f'SpeciesConc_ClNO3'] + ds0[f'SpeciesConc_ICl'] + ds0[f'SpeciesConc_BrCl']   \
                 + ds0[f'SpeciesConc_ClOO'] + ds0[f'SpeciesConc_OClO'] + (ds0[f'SpeciesConc_Cl2O2']*2)   \
                         ) * 1e12 )
        elif species=='Iy':
            ds0 = (( ds0[f'SpeciesConc_I'] + (ds0[f'SpeciesConc_I2']*2) + ds0[f'SpeciesConc_HOI'] \
                 + ds0[f'SpeciesConc_IO'] + ds0[f'SpeciesConc_OIO'] + ds0[f'SpeciesConc_HI']  \
                 + ds0[f'SpeciesConc_INO'] \
                 + (ds0[f'SpeciesConc_I2O2']*2) + (ds0[f'SpeciesConc_I2O3']*2) + (ds0[f'SpeciesConc_I2O4']*2)   \
                         ) * 1e12 )
        else:
            ds0 = ds0[f'SpeciesConc_{species}'] * 1e9

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
    labels=["Jscale_100", "Shah et al. 2022", "Shah - j100"]
    for n in range(0,3):
        print( n )
        if plot_type=='abs':
            all_plottable.append( ds[1] - ds[0] ) 
            unit = 'ppbv'
        elif plot_type=='pc':
            all_plottable.append( (ds[1] / ds[0] * 100)-100 ) 
            if first:
                sname=sname+'_PC'
                unit = '%'
        levels=None
        ax = fig.add_subplot(1,3,n+1, projection=ccrs.Robinson(), aspect='auto')
        if n < 2:
            print( 'to plot' )
            #Max= np.nanmax([all_plottable[0], all_plottable[1] ] )
            im = all_plottable[n].plot.imshow( x='lon',y='lat', ax=ax, transform=ccrs.PlateCarree(), 
                   add_colorbar=False)
            print( 'after plot' )
        else:
            im0 = all_plottable[n].plot.imshow( x='lon',y='lat', ax=ax, transform=ccrs.PlateCarree(),
                    cmap='bwr', center=0., add_colorbar=False) 

        ax.coastlines()
        ax.set_title(labels[n], fontsize=16)
    
    plt.tight_layout()
    plt.subplots_adjust( bottom=.1) 
    cbar_ax = fig.add_axes([0.02, 0.15, 0.6, 0.03])
    cbar = fig.colorbar(im, cax=cbar_ax, orientation='horizontal')
    cbar.ax.set_xlabel(f'{species} (ppbv)', size=15)

    cbar_ax = fig.add_axes([0.68, 0.15, 0.3, 0.03])
    cbar = fig.colorbar(im0, cax=cbar_ax, orientation='horizontal')
    cbar.ax.set_xlabel(f'$\Delta$ {species} ({unit})', size=15)
    
    plt.savefig( f'plots/{sname}.png' , dpi=300)
    plt.close()
    return

##----------------------------------------MAIN SCRIPT----------------------------------------##
def main():
    global d, site, rundirs, v
    versions = ['13.1.2']
    
    species = ['O3','NOx','O3','HNO2','CO','ALK4','NOx','NO','NO2','C2H6',
               'NIT','NITs','NITD1','NITD2','NITD3','NITD4',
               'SO4','SO4','SO4s','SO4D1','SO4D2','SO4D3','SO4D4',
               'OH',
               'all_NIT','all_SO4',
               'Cly','Bry','Iy']
    for s in species[:1]:
        print( s )
        rundirs  = ['base_run_2015','viral_run_2015']   
        ds_ = surface_data(s, rundirs, versions)

        cart_plot( s, rundirs, ds_, plot_type='pc' , sname=f'Shah.-viral_diff_{s}')
        cart_plot( s, rundirs, ds_, plot_type='abs', sname=f'Shah.-viral_diff_{s}')
        
if __name__ == "__main__":
    main()
