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
    print( path, substrs )
    sys.exit()
    for root, directory, files in os.walk(path):
        for f in files:
            for s in substrs:
                if s in f:
                    file_list.append(os.path.join(root, f))
    file_list.sort()
    return file_list

def get_data_as_xr(rundir, year='', lev=0, version='13.1.2', variable=False):
    path=f'/mnt/lustre/users/mjr583/GC/{version}/rundirs/{rundir}'
    print( path )
    sys.exit()
    if variable=='OH':
        ds = xr.open_mfdataset( find_file_list(path, [f'ConcAfterChem.{year}']), combine='by_coords' )
        ds = ds[f'OHconcAfterChem'].isel(lev=lev) #* 1e9
    else:
        ds = xr.open_mfdataset( find_file_list(path, [f'SpeciesConc.{year}']), combine='by_coords' )
        if variable:
            if variable=='NOx':
                ds = ( ds[f'SpeciesConc_NO'] + ds[f'SpeciesConc_NO2'] ).isel(lev=lev) * float(d[variable]['scale'])
            elif variable=='NITA':
                ds = ( ds[f'SpeciesConc_NIT'] + ds[f'SpeciesConc_NITs'] ).isel(lev=lev) * 1e9
            elif variable=='SO4A':
                ds = ( ds[f'SpeciesConc_SO4'] + ds[f'SpeciesConc_SO4s'] ).isel(lev=lev) * 1e9
            elif variable=='NOy':
                ds = ( ds[f'SpeciesConc_NO'] + ds[f'SpeciesConc_NO2'] + ds[f'SpeciesConc_N2O5'] \
                     + ds[f'SpeciesConc_ClNO2'] + ds[f'SpeciesConc_HNO3'] + ds[f'SpeciesConc_NIT']  \
                     + ds[f'SpeciesConc_NITs'] + ds[f'SpeciesConc_HNO2'] + ds[f'SpeciesConc_NH3']   \
                     + ds[f'SpeciesConc_PAN'] + ds[f'SpeciesConc_HNO4'] \
                             ).isel(lev=lev) * 1e12
            elif variable=='DST':
                ds = (  ds[f'SpeciesConc_DST1'] + ds['SpeciesConc_DST2'] + ds['SpeciesConc_DST3'] + \
                        ds['SpeciesConc_DST4'] ).isel(lev=lev) * 1e12
            elif variable == 'all_nitrate':
                try:
                    ds = ( ds[f'SpeciesConc_NIT']   + ds[f'SpeciesConc_NITs']  + \
                           ds[f'SpeciesConc_NITD1'] + ds[f'SpeciesConc_NITD2'] + \
                           ds[f'SpeciesConc_NITD3'] + ds[f'SpeciesConc_NITD4']).isel(lev=0) * 1e9
                    print( 'option A')
                except:
                    ds = ( ds[f'SpeciesConc_NIT'] + ds[f'SpeciesConc_NITs'] ).isel(lev=0) * 1e9
                    print( 'option B' )

            elif variable == 'all-nitrate-plus-HNO3':
                try:
                    ds = ( ds[f'SpeciesConc_NIT']   + ds[f'SpeciesConc_NITs']  + \
                           ds[f'SpeciesConc_NITD1'] + ds[f'SpeciesConc_NITD2'] + \
                           ds[f'SpeciesConc_NITD3'] + ds[f'SpeciesConc_NITD4'] + ds[f'SpeciesConc_HNO3']).isel(lev=0) * 1e9
                    print( 'option A')
                except:
                    ds = ( ds[f'SpeciesConc_NIT'] + ds[f'SpeciesConc_NITs'] + ds[f'SpeciesConc_HNO3'] ).isel(lev=0) * 1e9
                    print( 'option B' )

            elif variable == 'all_sulphate':
                try:
                    ds = ( ds[f'SpeciesConc_SO4']   + ds[f'SpeciesConc_SO4s']  + \
                           ds[f'SpeciesConc_SO4D1'] + ds[f'SpeciesConc_SO4D2'] + \
                           ds[f'SpeciesConc_SO4D3'] + ds[f'SpeciesConc_SO4D4']).isel(lev=0) * 1e9
                except:
                    ds = ( ds[f'SpeciesConc_SO4'] + ds[f'SpeciesConc_SO4s'] ).isel(lev=0) * 1e9
            else:
                ds = ds[f'SpeciesConc_{variable}'].isel(lev=lev) * float(d[variable]['scale'])
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
    #ax_Max[3] = 200
    for n in range(len(plottable)):
        m = rp.get_basemap(freq=60, ax=axes[n])
        #if n <2:
        #    im = m.pcolormesh(X, Y, plottable[n], cmap=cmap, vmax=np.max( [ax_Max[0], ax_Max[1] ] ))
        #else:
        im = m.pcolormesh(X, Y, plottable[n], cmap='bwr', vmax=ax_Max[n], vmin=-(ax_Max[n]),)
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


##----------------------------------------MAIN SCRIPT----------------------------------------##
def main():
    global d, site, rundirs, v
    rundir, variables, version = get_arguments()
    rundirs  = ['new_base_4x5', 'langmuir_andersen_4x5']
    versions = ['14.1.0','14.1.0']

    variables=['NH4','NOx','NITA','SO4A','NIT', 'NITs','SO4','all_nitrate','all_sulphate','O3','CO'][::-1]
    units = ['ppbv','ppbv','pptv','pptv','pptv','ppbv','ppbv','pptv','ppbv','ppbv']
    for nn, v in enumerate(variables):
        print( v )
        f, ((ax1,ax2)) = plt.subplots(1,2,figsize=(14,6))
        axes=[ax1,ax2]
        add_plottable=[]
        for n, rundir in enumerate(rundirs):
            ds0, lat, lon = get_data_as_xr(rundir, year='2014', lev=0, variable=v, version=versions[n])
            ds0 = np.mean( ds0, 0 )
            add_plottable.append( ds0 )

        add_plottable.append( add_plottable[1] - add_plottable[0] ) 
        add_plottable.append( - ( 100 - ( add_plottable[1] / add_plottable[0] * 100 ))) 
        plot( f, axes, add_plottable, lat, lon, c=cmap, 
                cbar_label=[f'{units[nn]}','%'],
                labels=['Langmuir RH - Base','Langmuir RH / Base'],
                sname=f'RH_Langmuir_{v}')
        
if __name__ == "__main__":
    main()
