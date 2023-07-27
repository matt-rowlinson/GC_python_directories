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

def get_data_as_xr(rundir, year='', lev=0, variable=False):
    path=f'/mnt/lustre/users/mjr583/GC/13.1.2/rundirs/{rundir}'
    ds = xr.open_mfdataset( find_file_list(path, [f'SpeciesConc.{year}']), combine='by_coords' )
    if variable:
        #try:
        #    ds = ds[f'SpeciesConc_{d[variable]["GC_name"]}'].isel(lev=lev) * float(d[variable]['scale'])
        #except:
        try:
            ds = ds[f'SpeciesConc_{d[variable]["GC_name"]}'].isel(lev=lev) * 1e12
        except:
            ds = ds[f'SpeciesConc_{variable}'].isel(lev=lev) * 1e12

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

def main():
    variables=['EOH']#ALK4','MEK','EOH','C2H6','C3H8','PRPE','BENZ','TOLU','XYLE','O3']
    rundirs = ['ceds_only','new_scale_all_vocs']
    for v in variables[::-1]:
        add_plottable=[] 

        print( v )
        dss=[]
        for n, rundir in enumerate(rundirs):
            ds0, lat, lon = get_data_as_xr(rundir, year='2017', lev=0, variable=v)
            ds0 = np.mean( ds0, 0 )
            dss.append( ds0 )
        add_plottable.append( dss[1] - dss[0] )
        add_plottable.append( - ( 100 - ( dss[1] / dss[0] * 100 )))

        f, ((ax1,ax2)) = plt.subplots(1,2,figsize=(10,4))
        axes=[ax1,ax2]
        ax_Max1 = plot_max_min(add_plottable[0])
        ax_Max2 = plot_max_min(add_plottable[1])
        ax_Max = [ax_Max1,ax_Max2]

        X, Y = np.meshgrid( lon, lat )
        panel_labels=['(a) Base - Scaled emissions','(b) Base / Scaled emissions']
        cbar_labels=['ppt','%']
        for n, ax in enumerate(axes):
            m = rp.get_basemap(freq=60, ax=ax, lines=False)
            Max = np.max(  (add_plottable[n].values.max(), -(add_plottable[n].values.min())) ) * 1.1
            Max=[1000,300]
            im = m.pcolormesh(X, Y, add_plottable[n], cmap='bwr', vmax= Max[n], vmin=-Max[n])
            ax.set_title(f'{panel_labels[n]}', fontsize=12 )
            
            divider = make_axes_locatable(ax)
            cax = divider.append_axes('bottom', size='5%', pad=0.25)
            cbar = f.colorbar(im, cax=cax, orientation='horizontal')
            cbar.ax.set_xlabel(cbar_labels[n],fontsize=12)
        plt.suptitle( f'{v}' )
        plt.tight_layout()
        plt.savefig( f"new_plots/Sfigs_{v}.png", dpi=200 )
        plt.close()

if __name__=="__main__":
    main()
