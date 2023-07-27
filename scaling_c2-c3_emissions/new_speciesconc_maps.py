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
'''
for variable in variables:
    control, lat, lon, lev, time = GC.get_gc_var('ceds_only', variable=variable, version='13.3.1', year=2017)
    fltr=np.where(control < 0.)
    control[fltr] = np.nan

    nohals, lat, lon, lev, time = GC.get_gc_var('scale_all_vocs', variable=variable,version='13.3.1', year='2017')
    fltr=np.where(nohals < 0.)
    nohals[fltr] = np.nan

    control=np.mean(control[:,0,:,:],0)
    nohals=np.mean(nohals[:,0,:,:],0)

    Min, Max = rp.get_abs_max_min( nohals / control )

    fig,(ax1, ax2, ax3)= plt.subplots(1,3,figsize=(12,3))
    X,Y=np.meshgrid(lon,lat)
    m1=rp.get_basemap(ax=ax1, lines=False, freq=60.)#lllon=-70, urlon=0., lllat=-15., urlat=30.,ax=ax1)
    m2=rp.get_basemap(ax=ax2, lines=False, freq=60.)#lllon=-70, urlon=0., lllat=-15., urlat=30.,ax=ax2)
    m3=rp.get_basemap(ax=ax3, lines=False, freq=60.)#lllon=-70, urlon=0., lllat=-15., urlat=30.,ax=ax3)

    MAX=np.nanmax(control)

    diff = - (100 - ( nohals / control ) * 100 )
    if variable == 'NO' or variable == 'NO2':

        im=ax1.pcolormesh(X,Y, control,  cmap=cmap, vmax=MAX, norm=LogNorm() )
        im2=ax2.pcolormesh(X,Y, nohals,  cmap=cmap, vmax=MAX, norm=LogNorm())
        im3=ax3.pcolormesh(X,Y, nohals / control,  cmap='bwr', vmin=0.97, vmax=1.03)

    else:
        im=ax1.pcolormesh(X,Y, control,  cmap=cmap, vmax=MAX)
        im2=ax2.pcolormesh(X,Y, nohals,  cmap=cmap, vmax=MAX)
        im3=ax3.pcolormesh(X,Y, diff,  cmap='bwr', vmin=-4, vmax=4)#, vmin=Min)

    ax1.set_title('Base')
    ax2.set_title('Scaled CEDS')
    ax3.set_title('Scaled CEDS / Base')

    plt.subplots_adjust(bottom=0.24)
    cbar_ax = fig.add_axes([0.126, 0.15, 0.502, 0.07])
    fig.colorbar(im, cax=cbar_ax, orientation='horizontal')
    cbar_ax.set_xlabel('%s (%s)' %(variable, d[variable]['unit']))

    cbar_ax = fig.add_axes([0.673,0.15,0.225, 0.07])
    fig.colorbar(im3, cax=cbar_ax, orientation='horizontal')
    cbar_ax.set_xlabel(f'%')#delta {variable}')

    #plt.suptitle(f'SpeciesConc_{variable}')

    plt.savefig(f'plots/SpeciesConc_map_{variable}.png')
    plt.close()
'''

def main():
    variables=['O3','OH']#,'
    rundirs = ['ceds_only','scale_all_vocs']

    for v in variables:
        print( v )
        f, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2,figsize=(12,8))
        axes=[ax1,ax2,ax3,ax4]
        add_plottable=[]
        for n, rundir in enumerate(rundirs):
            ds0, lat, lon = get_data_as_xr(rundir, year='2017', lev=0, variable=v)
            ds0 = np.mean( ds0, 0 )
            add_plottable.append( ds0 )
        add_plottable.append( add_plottable[1] - add_plottable[0] ) 
        add_plottable.append( - ( 100 - ( add_plottable[1] / add_plottable[0] * 100 ))) 

        plot( f, axes, add_plottable, lat, lon, c=cmap, 
                cbar_label=[f'{v} {d[v]["unit"]}',
                            f'{v} {d[v]["unit"]}',
                            u'Δ'+f'{v} {d[v]["unit"]}',
                            u'Δ%'],
                labels=['Base', 'Scaled VOC emissions', 'Scaled - Base','Scaled / Base'],
                sname=f'4_panel_scaled-ceds_{v}')


if __name__=="__main__":
    main()
