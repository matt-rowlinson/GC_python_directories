#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#SBATCH --job-name=jscale_time
#SBATCH --ntasks=1
#SBATCH --mem=2Gb
#SBATCH --partition=interactive
#SBATCH --time=00:30:00
#SBATCH --output=Logs/jscale_time_plots_%A.log
import sys
import os
import xarray as xr
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
import constants as const
import argparse
import matplotlib.colors as mcolors
plt.style.use('seaborn-darkgrid')

import yaml
spec_db = yaml.load(open("/mnt/lustre/users/mjr583/GC/13.1.2/rundirs/nitrate_photol_2010_control/species_database.yml"), Loader=yaml.FullLoader)

current_dir = os.path.dirname(__file__)
rgb_WhGrYlRd = np.genfromtxt('/users/mjr583/python_lib/colormaps/WhGrYlRd.txt',delimiter=' ')
WhGrYlRd = mcolors.ListedColormap(rgb_WhGrYlRd/255.0)
cmap=WhGrYlRd

def new_cmap():
    # sample the colormaps that you want to use. Use 128 from each so we get 256
    colors1 = plt.cm.bwr(np.linspace(0., .5, 128))
    colors2 = plt.cm.afmhot_r(np.linspace(0., .1, 8))
    colors3 = plt.cm.RdYlBu_r(np.linspace(.5, 1, 120))

    colors = np.vstack((colors1, colors2, colors3))
    mymap = mcolors.LinearSegmentedColormap.from_list('my_colormap', colors)
    return mymap

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

def get_data_as_xr(rundir, mw_kg, year='', lev=0, variable=False):
    path=f'/mnt/lustre/users/mjr583/GC/13.1.2/rundirs/{rundir}'
    ds = xr.open_mfdataset( find_file_list(path, [f'SpeciesConc.{year}']), combine='by_coords' )
    if variable:
        if variable == 'NOx':
            ds = ( ds[f'SpeciesConc_NO'] + ds[f'SpeciesConc_NO2'] ).isel(lev=lev) * float(d[variable]['scale'])
        else:
            if lev=='sum':
                if variable=='SO4_all':
                    ds = ds[f'SpeciesConc_SO4'] + ds[f'SpeciesConc_SO4s'] + ds[f'SpeciesConc_SO4H1'] +ds[f'SpeciesConc_SO4H2']
                    ds = ds * 1e9 * ( mw_kg / 224. ) * 1e3
                    ds = ds.sum(dim='lev')
                elif variable=='NIT_all':
                    ds = ds[f'SpeciesConc_NIT'] + ds[f'SpeciesConc_NITs']
                    ds = ds * 1e9 * ( mw_kg / 224. ) * 1e3
                    ds = ds.sum(dim='lev')
                else:
                    ds = ds[f'SpeciesConc_{variable}'] * 1e9 * ( mw_kg / 224. ) * 1e3
                    ds = ds.sum(dim='lev')
            else:
                ds = ds[f'SpeciesConc_{variable}'].isel(lev=lev) * 1e9 * ( mw_kg / 224. ) #* 1e3
    lat = ds['lat']
    lon = ds['lon']
    return ds, lat, lon

def plot_max_min(plottable):
    ax_Max=[]
    for i in plottable:
        Min, Max = rp.get_abs_max_min(i)
        ax_Max.append( Max )
    return ax_Max

def plot(f, axes, plottable, lat, lon, c='', labels=False,
        cbar_label=[],sname="4_panel_plot.png",
        panel_labels=['(a)','(b)','(c)','(d)']):

    X, Y = np.meshgrid( lon, lat )
    ax_Max = plot_max_min(plottable)
    for n in range(len(plottable)):
        m = rp.get_basemap(freq=60, ax=axes[n])
        if n <2:
            c.set_bad(color='w')
            im = m.pcolormesh(X, Y, plottable[n], cmap=c, vmax=np.max( [ax_Max[0], ax_Max[1] ] ))
        else:
            cmap = new_cmap()
            im = m.pcolormesh(X, Y, plottable[n], cmap=cmap, vmax=ax_Max[n], vmin=-(ax_Max[n]),)
        if labels:
            axes[n].set_title(f'{panel_labels[n]} {labels[n]}', fontsize=14 )
        
        divider = make_axes_locatable(axes[n])
        cax = divider.append_axes('bottom', size='5%', pad=0.25)
        cbar = f.colorbar(im, cax=cax, orientation='horizontal')
        cbar.ax.set_xlabel(f'{cbar_label[n]}',fontsize=12)

    plt.tight_layout()
    plt.savefig( f"plots/{sname}" )
    plt.close()
    return


def column_aerosol():
    variables=['SO4_all','NIT_all','SO4','SO4s','SO4H1','SO4H2','NITs','NIT']
    for v in variables:
        if v=='SO4_all':
            mw_kg = spec_db['SO4']["MW_g"] * 1.0e-3
        elif v=='NIT_all':
            mw_kg = spec_db['NIT']["MW_g"] * 1.0e-3
        else:
            mw_kg = spec_db[v]["MW_g"] * 1.0e-3

        f, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2,figsize=(12,8))
        axes=[ax1,ax2,ax3,ax4]
        add_plottable=[]
            
        ds0, lat, lon = get_data_as_xr('nitrate_photol_2000_control', mw_kg, year='2000', lev='sum', variable=v)
        ds0 = np.mean( ds0, 0 ).values
        add_plottable.append( ds0 )

        ds1, lat, lon = get_data_as_xr('nitrate_photol_2015_control', mw_kg, year='2020', lev='sum', variable=v)
        ds1 = np.mean( ds1, 0 ).values
        add_plottable.append( ds1 )

        add_plottable.append( add_plottable[1] - add_plottable[0] ) 
        add_plottable.append( - ( 100 - ( add_plottable[1] / add_plottable[0] * 100 ))) 

        plot( f, axes, add_plottable, lat, lon, c=cmap, 
                cbar_label=[f'{v} mg m-2', f'{v} mg m-2', f'{v} mg m-2', '%'],
                labels=[f'2000 {v}', f'2020 {v}', f'2020-2000','2020/2000'],
                sname=f'bauer_4_panel_2020-2000_{v}')
        return

##----------------------------------------MAIN SCRIPT----------------------------------------##
def main():
    global d, site, rundirs, v
    from mean_oh import get_data_as_xr as metrics
    rundir = 'nitrate_photol_2017_control'

    variables=['SO4','NOx','SO2','HNO3','NH4','NH3','NIT']
    for v in variables:
        f, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2,figsize=(12,8))
        axes=[ax1,ax2,ax3,ax4]
        add_plottable=[]
            
        print( v )
        if v=='NOx':
            mw_kg = spec_db['NO2']["MW_g"] * 1.0e-3
        else:
            mw_kg = spec_db[v]["MW_g"] * 1.0e-3

            
        ds0, lat, lon = get_data_as_xr('nitrate_photol_2000_control', mw_kg, year='2000', lev=0, variable=v)
        ds0 = np.mean( ds0, 0 ).values
        add_plottable.append( ds0 )

        ds1, lat, lon = get_data_as_xr('nitrate_photol_2015_control', mw_kg, year='2020', lev=0, variable=v)
        ds1 = np.mean( ds1, 0 ).values
        add_plottable.append( ds1 )

        add_plottable.append( add_plottable[1] - add_plottable[0] ) 
        add_plottable.append( - ( 100 - ( add_plottable[1] / add_plottable[0] * 100 ))) 

        plot( f, axes, add_plottable, lat, lon, c=cmap, 
                cbar_label=[f'{v} mg m-2', f'{v} mg m-2', f'{v} mg m-2', '%'],
                labels=[f'2000 {v}', f'2020 {v}', f'2020-2000','2020/2000'],
                sname=f'surface_bauer_4_panel_2020-2000_{v}')
        
if __name__ == "__main__":
    main()
