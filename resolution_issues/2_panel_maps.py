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
import xesmf as xe
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

def get_data_as_xr(rundir, year='', lev=0, variable=False):
    path=f'/mnt/lustre/users/mjr583/GC/13.1.2/rundirs/{rundir}'
    if variable=='OH':
        ds = xr.open_mfdataset( find_file_list(path, [f'ConcAfterChem.{year}']), combine='by_coords' )
        ds = ds[f'OHconcAfterChem'].isel(lev=lev) #* 1e9
    else:
        ds = xr.open_mfdataset( find_file_list(path, [f'SpeciesConc.{year}']), combine='by_coords' )
        if variable:
            if variable=='NOx':
                ds = ( ds[f'SpeciesConc_NO'] + ds[f'SpeciesConc_NO2'] ).isel(lev=lev) * float(d[variable]['scale'])
            elif variable=='NOy':
                ds = ( ds[f'SpeciesConc_NO'] + ds[f'SpeciesConc_NO2'] + ds[f'SpeciesConc_N2O5'] \
                     + ds[f'SpeciesConc_ClNO2'] + ds[f'SpeciesConc_HNO3'] + ds[f'SpeciesConc_NIT']  \
                     + ds[f'SpeciesConc_NITs'] + ds[f'SpeciesConc_HNO2'] + ds[f'SpeciesConc_NH3']   \
                     + ds[f'SpeciesConc_PAN'] + ds[f'SpeciesConc_HNO4'] \
                             ).isel(lev=lev) * 1e9
            elif variable=='DST':
                ds = (  ds[f'SpeciesConc_DST1'] + ds['SpeciesConc_DST2'] + ds['SpeciesConc_DST3'] + \
                        ds['SpeciesConc_DST4'] ).isel(lev=lev) * 1e9
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
        panel_labels=['(a)','(b)','(c)','(d)']):

    X, Y = np.meshgrid( lon, lat )
    ax_Max = plot_max_min(plottable)
    for n in range(len(plottable)):
        m = rp.get_basemap(freq=60, ax=axes[n])
        if n <2:
            im = m.pcolormesh(X, Y, plottable[n], cmap=cmap, vmax=np.max( [ax_Max[0], ax_Max[1] ] ))
        else:
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

def regrid(ds_in, ds_out):
    ds_out = xr.Dataset(
            {   "lat": (["lat"], ds_out['lat']),
                "lon": (["lon"], ds_out['lon']),  } )
    regridder = xe.Regridder(ds_in, ds_out, "bilinear")
    ds = regridder( ds_in )
    return ds

##----------------------------------------MAIN SCRIPT----------------------------------------##
def main():
    global d, site, rundirs, v
    rundir, variables, version = get_arguments()
    rundirs = ['ryan_low_res_CVAO', f'ryan_high_res_CVAO']

    variables=['NIT','NITs','NOx','NOy']
    for v in variables:
        print( v )
        f, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2,figsize=(12,8))
        axes=[ax1,ax2,ax3,ax4]
        add_plottable=[]

        ds1, Hlat, Hlon = get_data_as_xr(rundirs[1], year='201706', lev=0, variable=v)
        ds1 = np.mean( ds1, 0 )

        ds0, lat, lon = get_data_as_xr(rundirs[0], year='201706', lev=0, variable=v)
        ds0 = np.mean( ds0, 0 )
        ds0 = ds0.sel( lon=slice(Hlon.min()-5, Hlon.max()+5) )
        ds0 = ds0.sel( lat=slice(Hlat.min()-5, Hlat.max()+5) )
        lat = lat.sel( lat=slice(Hlat.min()-5, Hlat.max()+5) )
        lon = lon.sel( lon=slice(Hlon.min()-5, Hlon.max()+5) )

        ds0_high = regrid( ds0, ds1 )
        ds1_low = regrid( ds1, ds0 )

        Max= np.max([ds0.values.max(), ds1.values.max()])
        X, Y = np.meshgrid( lon, lat)
        HX, HY = np.meshgrid( Hlon, Hlat) 
        m1 = rp.get_basemap(lllon=Hlon.min(), urlon=Hlon.max(), lllat=Hlat.min(), urlat=Hlat.max(), freq=60, ax=ax1)
        im = m1.pcolormesh(X, Y, ds0, cmap=cmap, vmin=0., vmax=Max )
        
        m2 = rp.get_basemap(lllon=Hlon.min(), urlon=Hlon.max(), lllat=Hlat.min(), urlat=Hlat.max(), freq=60, ax=ax2)
        im = m2.pcolormesh(HX, HY, ds1, cmap=cmap, vmin=0., vmax=Max )
        
        Max= np.max([(ds0.values - ds1_low.values).max(), -(ds0.values - ds1_low.values).min()])
        m3 = rp.get_basemap(lllon=Hlon.min(), urlon=Hlon.max(), lllat=Hlat.min(), urlat=Hlat.max(), freq=60, ax=ax3)
        im = m3.pcolormesh(X, Y, ds0 - ds1_low, cmap='bwr', vmin=-Max, vmax=Max )
        
        m4 = rp.get_basemap(lllon=Hlon.min(), urlon=Hlon.max(), lllat=Hlat.min(), urlat=Hlat.max(), freq=60, ax=ax4)
        im = m4.pcolormesh(HX, HY, ds0_high - ds1, cmap='bwr', vmin=-Max, vmax=Max )

        idy = rp.find_nearest(lat, 16.9)
        idx = rp.find_nearest(lon, -24.9)
        m1.scatter( lon[idx],lat[idy], marker='x', color='k', s=100)
        m3.scatter( lon[idx],lat[idy], marker='x', color='k', s=100)

        idy = rp.find_nearest(Hlat, 16.9)
        idx = rp.find_nearest(Hlon, -24.9)
        m2.scatter( Hlon[idx],Hlat[idy], marker='x', color='k', s=120)
        m4.scatter( Hlon[idx],Hlat[idy], marker='x', color='k', s=120)

        plt.tight_layout()
        plt.savefig( f"plots/CVAO_high_res_{v}.png" )
        plt.close()
        
if __name__ == "__main__":
    main()
