#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#SBATCH --job-name=hemco_diffmap
#SBATCH --ntasks=1
#SBATCH --mem=18gb
#SBATCH --partition=interactive
#SBATCH --time=08:00:00
#SBATCH --output=Logs/hemco_diffmap_%A.log
import xarray as xr
import numpy as np
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('agg')
from matplotlib.colors import LogNorm
import os
import sys

def find_file_list(path, substrs):
    file_list = []
    for root, directory, files in os.walk(path):
        for f in files:
            for s in substrs:
                if s in f:
                    file_list.append(os.path.join(root, f))
    file_list.sort()
    return file_list

def read_hemco(rundir, year):
    path=f'/users/mjr583/scratch/GC/13.1.2/rundirs/{rundir}/OutputDir/'
    flist=find_file_list(path, [f'HEMCO_diagnostics.{year}'])
    ds = xr.open_mfdataset(flist)
    return ds


def main():
    b1980 = read_hemco('base_run_1980', '1980').sum(dim='time')*1e9
    b1990 = read_hemco('base_run_1990', '1990').sum(dim='time')*1e9
    b2000 = read_hemco('base_run_2000', '2000').sum(dim='time')*1e9
    b2005 = read_hemco('base_run_2005', '2005').sum(dim='time')*1e9
    b2010 = read_hemco('base_run_2010', '2010').sum(dim='time')*1e9
    b2015 = read_hemco('base_run_2015', '2015').sum(dim='time')*1e9
    b2020 = read_hemco('base_run_2020', '2020').sum(dim='time')*1e9

    keys=sorted(list(set(b1980.keys())))
    for key in keys[:-4][::-1]:
        if b1980[key].max() == 0.:
            print( key, "No emissions - skipping" )
            continue
        print( key )
        
        p = [
            b1980[key], b1980[key],
            b1990[key], (b1990[key]-b1980[key])*100,
            b2000[key], (b2000[key]-b1990[key])*100,
            b2005[key], (b2005[key]-b2000[key])*100,
            b2010[key], (b2010[key]-b2005[key])*100,
            b2015[key], (b2015[key]-b2010[key])*100,
            b2020[key], (b2020[key]-b2015[key])*100
            ]

        f = plt.figure(figsize=(5,9))
        for n, pp in enumerate(p):
            try:
                pp=pp.isel(lev=0)
            except:
                pass 
            ax = f.add_subplot(7,2,n+1, projection=ccrs.Robinson(), aspect='auto')
            ax.coastlines()
            
            if n % 2 == 0:
                vmin=0 
                vmax=np.nanmax([p[0],p[2],p[4],p[6],p[8],p[10],p[12]])
                im0 = pp.plot.imshow( x='lon',y='lat', ax=ax, norm=LogNorm(), vmin=0., 
                    vmax=np.nanmax([p[0],p[2],p[4],p[6],p[8],p[10],p[12]]), 
                    transform=ccrs.PlateCarree(), center=0., add_colorbar=False)
            else:
                Min=np.nanmin([p[3],p[5],p[7],p[9],p[11],p[13]])
                Max=np.nanmax([p[3],p[5],p[7],p[9],p[11],p[13]])
                lim = np.max([ abs(Min), abs(Max) ]) * .1
                im1 = pp.plot.imshow( x='lon',y='lat', ax=ax, vmin=-lim, vmax=lim,
                    transform=ccrs.PlateCarree(), cmap='bwr', center=0., add_colorbar=False)

            ax.set_title( None)
            if n==1:
                plt.delaxes(ax)
        plt.tight_layout()
        plt.subplots_adjust( bottom=.08) 
        #cbar_ax = f.add_axes([0.05, 0.04, 0.6, 0.02])
        #cbar = f.colorbar(im0, cax=cbar_ax, orientation='horizontal')
        #cbar.ax.set_xlabel(f'kg s-1 m-2', size=15)

        cbar_ax = f.add_axes([0.6, 0.06, 0.3, 0.01])
        cbar = f.colorbar(im1, cax=cbar_ax, orientation='horizontal')
        cbar.ax.set_xlabel('%', size=15)
        
        f.text( .5, .9, key, fontsize=14, weight='bold' )

        plt.savefig(f'hemco_plots/HEMCO_diff_{key}.png' )
        plt.close()


if __name__=="__main__":
    main()
