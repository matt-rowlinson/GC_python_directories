#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#SBATCH --job-name=global_ratios
#SBATCH --ntasks=1
#SBATCH --mem=8Gb
#SBATCH --partition=interactive
#SBATCH --time=00:45:00
#SBATCH --output=Logs/global_ratios_%A.log
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
import matplotlib.colors as colors
import os
import tomas_custom_cmap

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

def get_data_as_xr(rundir, year, lonlat=False):
    path=f'/mnt/lustre/users/mjr583/GC/13.1.2/rundirs/{rundir}'
    ds = xr.open_mfdataset( find_file_list(path, [f'SpeciesConc.{year}']), combine='by_coords' )
    nox = ds['SpeciesConc_NO'] + ds['SpeciesConc_NO2'] * 1e12
    nox = nox.isel(lev=slice(0,9)).mean(dim='lev', keep_attrs=True).mean(dim='time', keep_attrs=True)
    o3 = ds['SpeciesConc_O3'].isel(lev=slice(0,9)).mean(dim='lev', keep_attrs=True).mean(dim='time', keep_attrs=True) * 1e9
    oh = xr.open_mfdataset( find_file_list( path, [f'ConcAfterChem.{year}']), 
            combine='by_coords')['OHconcAfterChem'].isel(lev=slice(0,9)).mean(dim='lev',
            keep_attrs=True).mean(dim='time',keep_attrs=True)
    if lonlat==True:
        return nox, oh, o3, ds['lon'], ds['lat']
    return nox, oh, o3

def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = colors.LinearSegmentedColormap.from_list(
            'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
            cmap(np.linspace(minval, maxval, n)))
    return new_cmap

def main():
    ## Get the data for each year/scale
    years=['1980','2000','2010','20175']
    for year in years:
        print(year)
        year2=year#+"01"
        con_nox,  con_oh,  con_o3, lon, lat  = get_data_as_xr(f'nitrate_photol_{year}_control', year2, lonlat=True)
        s25_nox,  s25_oh,  s25_o3  = get_data_as_xr(f'nitrate_photol_{year}_all_scale-25', year2)
        #s100_nox, s100_oh, s100_o3 = get_data_as_xr(f'nitrate_photol_{year}_all_scale-', year2)
        
        ## Do plotting (3x3 grid)
        to_plot = [con_nox, con_oh, con_o3, s25_nox / con_nox, s25_oh / con_oh, s25_o3 / con_o3,s100_nox / con_nox, s100_oh / con_oh, s100_o3 / con_o3 ]
        subtitles = [r'$J_{0scale} NO_x$',r'$J_{0scale} OH$',r'$J_{0scale} O_3$',
                r'$J_{25scale} NO_x$ / $J_{0scale} NO_x$',r'$J_{25scale} OH$ / $J_{0scale} OH$',r'$J_{25scale} O_3$ / $J_{0scale} O_3$',
                r'$J_{100scale} NO_x$ / $J_{0scale} NO_x$',r'$J_{100scale} OH$ / $J_{0scale} OH$',r'$J_{100scale} O_3$ / $J_{0scale} O_3$']
        labels=['(a)','(b)','(c)','(d)','(e)','(f)','(g)','(h)','(i)'] 
        norms=[LogNorm(), LogNorm(), None, LogNorm(), None, None, LogNorm(), None, None]
        mins=[1, 1e5, 0, 1, 1, 1, 1, 1, 1]
        maxs=[1000, 1e7, 50, 20, 1.6, 1.30, 20, 1.6, 1.30]

        cbar_labels=['pptv','molecule cm-3', 'ppbv','','','','','','']
        
        fig = plt.figure(figsize=(10,8))
        X, Y = np.meshgrid( lon, lat )

        cmap = tomas_custom_cmap.get_colormap('TMS_custom_colormap_CMRmap.cpt')
        #cmap = plt.cm.inferno_r
        #cmap = truncate_colormap(cmap, 0., 0.7)
        #cmap.set_under('w')
        for n, ax in enumerate(range(len(to_plot))):
            #print( np.nanmean(to_plot[n]) )
            ax = fig.add_subplot(3,3,n+1)
            m = rp.get_basemap(proj='robin', ax=ax, lines=False)
            xx,yy=m(X,Y)
            im = m.pcolormesh( xx, yy, to_plot[n], cmap=cmap, vmin=mins[n],vmax=maxs[n], norm=norms[n],zorder=100)

            m.drawcoastlines(zorder=101)
            m.fillcontinents(color='grey', zorder=101)
            ax.title.set_text(subtitles[n])
            ax.text(-0.03, 1.05, labels[n], fontsize=12, transform=ax.transAxes, zorder=103)

            divider = make_axes_locatable(ax)
            cax = divider.append_axes('bottom', size='5%', pad=0.05)
            cbar=fig.colorbar(im, cax=cax, orientation='horizontal', extend='both')
            cbar.ax.set_xlabel(cbar_labels[n])

        fig.suptitle(year2) 
        plt.savefig(f'plots/global_{year2}.png')
        plt.close()
        print('Done')

if __name__=="__main__":
    main()
