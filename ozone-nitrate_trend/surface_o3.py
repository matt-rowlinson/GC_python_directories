#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#SBATCH --job-name=NOxtimeseries
#SBATCH --ntasks=1
#SBATCH --mem=18Gb
#SBATCH --partition=nodes
#SBATCH --time=00:45:00
#SBATCH --output=Logs/timeseries_%A.log
import sys
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
matplotlib.use('agg')
sys.path.append('/users/mjr583/python_lib')
import GC_tools as GC
import RowPy as rp
import os

def find_file_list(path, substrs):
    file_list =[]
    for root, directory, files in os.walk(path):
        for f in files:
            for s in substrs:
                if s in f:
                    file_list.append(os.path.join(root, f))
    file_list.sort()
    return file_list

def get_inputs():
    inputs=GC.get_arguments()
    variable=inputs.var
    if variable==None:
        variable='O3'
    return variable

def get_data_as_xr(rundir, year, lonlat=False, area=False):
    path=f'/mnt/lustre/users/mjr583/GC/13.1.2/rundirs/{rundir}'
    ds = xr.open_mfdataset( find_file_list(path, [f'SpeciesConc.{year}']), combine='by_coords' )
    o3 = ds['SpeciesConc_O3'].isel(lev=0).mean(dim='time', keep_attrs=True) * 1e9
    if lonlat==True:
        if area==True:
            return o3, ds['lon'], ds['lat'], ds['AREA']
        return o3, ds['lon'], ds['lat']      
    return o3

def get_cosweights(lat):
    cosweights = np.cos(np.deg2rad(lat))
    cosweights = np.swapaxes((np.tile(cosweights,72)).reshape(72,46),0,1)
    return cosweights

def main():
    ## Get the data for each year/scale
    con_o3, lon, lat, area  = get_data_as_xr(f'nitrate_photol_1980_control', '1980', lonlat=True, area=True)
    cosweights = get_cosweights(lat)
    con_o3 = con_o3 * cosweights
    s25_o3  = get_data_as_xr(f'nitrate_photol_1980_scale-25', '1980') * cosweights
    s100_o3 = get_data_as_xr(f'nitrate_photol_1980_scale-100', '1980') * cosweights

    con_o3_2010, lon, lat, area  = get_data_as_xr(f'nitrate_photol_2010_control', '2010', lonlat=True, area=True)
    con_o3_2010 = con_o3_2010 * cosweights
    s25_o3_2010  = get_data_as_xr(f'nitrate_photol_2010_scale-25', '2010') * cosweights
    s100_o3_2010 = get_data_as_xr(f'nitrate_photol_2010_scale-100', '2010') * cosweights

    con_o3_2017, lon, lat, area  = get_data_as_xr(f'nitrate_photol_2017_control', '2017', lonlat=True, area=True)
    con_o3_2017 = con_o3_2017 * cosweights
    s25_o3_2017  = get_data_as_xr(f'nitrate_photol_2017_scale-25', '2017') * cosweights
    s100_o3_2017 = get_data_as_xr(f'nitrate_photol_2017_scale-100', '2017') * cosweights
    
    ## Do plotting (3x3 grid)
    to_plot = [con_o3, s25_o3 - con_o3, s100_o3 - con_o3,
                con_o3_2010,  s25_o3_2010 - con_o3_2010, s100_o3_2010 - con_o3_2010,
                con_o3_2017,  s25_o3_2017 - con_o3_2017, s100_o3_2017 - con_o3_2017]

    subtitles = ['Control',r'$J_{25scale}$ - Control',r' $J_{100scale}$ - Control','','','','','','',
                 '2010 control',r'2010 $J_{25scale}$ - control',r'2010 $J_{100scale}$ - control',
                 '2017 control',r'2017 $J_{25scale}$ - control',r'2017 $J_{100scale}$ - control']

    labels=['(a)','(b)','(c)','(d)','(e)','(f)','(g)','(h)','(i)'] 
    norms=[LogNorm(), LogNorm(), None, LogNorm(), None, None, LogNorm(), None, None]
    mins=[0, -15, -15, 0, -15, -15, 0, -15, -15]
    maxs=[50, 15, 15, 50, 15, 15, 50, 15, 15]

    cbar_labels=['ppbv','ppbv', 'ppbv','','','','','','']
    ylabels=['1980','','','2010','','','2017','','']
    cmaps = [plt.cm.inferno_r,plt.cm.bwr,plt.cm.bwr,plt.cm.inferno_r,plt.cm.bwr,plt.cm.bwr,plt.cm.inferno_r,plt.cm.bwr,plt.cm.bwr] 
    fig = plt.figure(figsize=(10,8))
    X, Y = np.meshgrid( lon, lat )
    #cmap = plt.cm.inferno_r
    #cmap.set_under('w')
    for n, ax in enumerate(range(len(to_plot))):
        ax = fig.add_subplot(3,3,n+1)
        m = rp.get_basemap(proj='robin', ax=ax, lines=False)
        xx,yy=m(X,Y)
        cmaps[n].set_over('k')
        im = m.pcolormesh( xx, yy, to_plot[n], cmap=cmaps[n], vmin=mins[n],vmax=maxs[n],zorder=100)

        m.drawcoastlines(zorder=101)
        #m.fillcontinents(color='grey', zorder=101)
        ax.title.set_text(subtitles[n])
        ax.text(-0.03, 1.05, labels[n], fontsize=12, transform=ax.transAxes, zorder=103)

        divider = make_axes_locatable(ax)
        cax = divider.append_axes('bottom', size='5%', pad=0.05)
        cbar=fig.colorbar(im, cax=cax, orientation='horizontal', extend='max')
        cbar.ax.set_xlabel(cbar_labels[0])
        cbar.ax.set_ylabel(ylabels[n])


    #fig.suptitle(year) 
    plt.savefig(f'plots/surface_o3.png')
    plt.close()
    print('Done')

if __name__=="__main__":
    main()
