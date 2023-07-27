#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#SBATCH --job-name=jscale_maps
#SBATCH --ntasks=1
#SBATCH --mem=150Gb
#SBATCH --partition=himem
#SBATCH --time=05:00:00
#SBATCH --output=Logs/jscale_map_plots_%A.log
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
        ds = ds[f'OHconcAfterChem'].isel(lev=lev) 
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

def surface_data(rundirs, versions=['13.1.2']):
    add_=[]
    for n, rundir in enumerate(rundirs):
        print( rundir )
        ds0, lat, lon = get_data_as_xr(rundir, year='', version=versions[0])
        ds0 = ds0.mean('time')
        ds0 = ds0.isel( lev=0 )
        add_.append( (ds0['SpeciesConc_NO'] + ds0['SpeciesConc_NO2'] )* 1e9 ) 
    return add_

def make_plottable( ds_, interval, plot_type='abs'):
    plottable_=[]
    for n in range( 1, 7 ):
        if plot_type=='abs':
            plottable_.append( ( ds_[n] - ds_[n-1] ) / interval[n-1] )
        elif plot_type=='pc':
            plottable_.append((- ( 100 - ( ds_[n] / ds_[n-1] * 100 ))) / interval[n])
    return plottable_

def the_diff( a, b ):
    c_=[]
    for n, aa in enumerate( a ):
        c_.append( b[n] - a[n] )
    return c_

def cart_plot(rundirs, base, jrun, interval, plot_type='abs', sname=''):
    fig=plt.figure(figsize=(12,12))
    first=True
    base_ = make_plottable( base, interval, plot_type=plot_type )
    jrun_ = make_plottable( jrun, interval, plot_type=plot_type )
    diff_of_diff = the_diff( base_, jrun_ )
    xss = [base_, jrun_, diff_of_diff]
    all_plottable = [x for xs in xss for x in xs]
    
    ax_placement=[1, 4, 7, 10, 13, 16, 2, 5, 8, 11, 14, 17, 3, 6, 9, 12, 15, 18]
    labels = ['1990 - 1980', '2000 - 1990', '2005 - 2000', '2010 - 2005', '2015 - 2010','2020 - 2015',
              '1990 - 1980', '2000 - 1990', '2005 - 2000', '2010 - 2005', '2015 - 2010','2020 - 2015',
              '(b)-(a)', '(e)-(d)','(h)-(g)','(k)-(j)','(n)-(m)','(q)-(p)']
    fig_n=['(a)', '(b)', '(c)', '(d)', '(e)', '(f)', '(g)', '(h)', '(i)', '(j)', '(k)', '(l)', '(m)',
            '(n)','(o)','(p)','(q)', '(r)']
           
    for n in range(0,18):
        if plot_type=='abs':
            unit = r'$\Delta$ total NOx (ppbv/yr$^{-1}$)'
            if n < 11:
                levels=np.arange( -4, 4.1, .1 )
            else:
                levels=np.arange( -1,1.01, .01 )
        elif plot_type=='pc':
            if n < 11: 
                levels=np.arange( -18, 18.5, .5 )
            else:
                levels=np.arange( -10,10.2, .2 )
            if first:
                sname+='_PC'
                unit = r'$\Delta$ total NOx (%/yr$^{-1}$)'

        ax = fig.add_subplot(6,3,ax_placement[n], projection=ccrs.Robinson(), aspect='auto')
        if first:
            im0 = all_plottable[n].plot.imshow( x='lon',y='lat', ax=ax, transform=ccrs.PlateCarree(), 
                    levels=levels, center=0., add_colorbar=False)
            first=False
        else:
            im = all_plottable[n].plot.imshow( x='lon',y='lat', ax=ax, transform=ccrs.PlateCarree(), 
                    levels=levels, center=0., add_colorbar=False) 
        print( n, all_plottable[n].values.max() )

        #cbar_kwargs={'orientation':'horizontal','shrink':0.6, 'aspect':40,'label':f'{unit}'})
        ax.coastlines()
        plt.title( None)#"{fig_n[n]} {labels[n]}", fontsize=12)
    plt.tight_layout()
    plt.subplots_adjust( bottom=.08) 
    cbar_ax = fig.add_axes([0.05, 0.04, 0.6, 0.02])
    cbar = fig.colorbar(im0, cax=cbar_ax, orientation='horizontal')
    cbar.ax.set_xlabel(f'{unit}')

    cbar_ax = fig.add_axes([0.67, 0.04, 0.3, 0.02])
    cbar = fig.colorbar(im, cax=cbar_ax, orientation='horizontal')
    cbar.ax.set_xlabel(f'{unit}')
    
    plt.subplots_adjust( top=.95)
    fig.text( .14, .97, "Base", fontsize=23, weight='bold' )
    fig.text( .45, .97, "Jscale$_{100}$", fontsize=23, weight='bold' )
    fig.text( .72, .97, "Jscale$_{100}$ - Base", fontsize=23, weight='bold' )

    fig.text( .02, .84, "1990-1980", rotation='vertical', fontsize=16, weight='bold' )
    fig.text( .02, .68, "2000-1990", rotation='vertical', fontsize=16, weight='bold' )
    fig.text( .02, .53, "2005-2000", rotation='vertical', fontsize=16, weight='bold' )
    fig.text( .02, .39, "2010-2005", rotation='vertical', fontsize=16, weight='bold' )
    fig.text( .02, .24, "2015-2010", rotation='vertical', fontsize=16, weight='bold' )
    fig.text( .02, .1, "2020-2015", rotation='vertical', fontsize=16, weight='bold' )
    
    plt.savefig( f'plots/{sname}.png' , dpi=300)
    plt.close()
    return

##----------------------------------------MAIN SCRIPT----------------------------------------##
def main():
    global d, site, rundirs, v
    versions = ['13.1.2']
    interval = [10, 10, 5, 5, 5, 5, 5]
    v='NOx'
    print( v )

    rundirs  = ['base_run_1980', 'base_run_1990','base_run_2000', 'base_run_2005',
                                 'base_run_2010', 'base_run_2015', 'base_run_2020']   
    base_ = surface_data(rundirs, versions)
    

    rundirs  = ['j100_run_1980', 'j100_run_1990','j100_run_2000', 'j100_run_2005',
                                 'j100_run_2010', 'j100_run_2015', 'j100_run_2020']  
    j100_ = surface_data(rundirs, versions)
    
    cart_plot( rundirs, base_, j100_, interval, plot_type='pc', sname='18_panel_NOx_diffs')
    cart_plot( rundirs, base_, j100_, interval, plot_type='abs', sname='18_panel_NOx_diffs')
        
if __name__ == "__main__":
    main()
