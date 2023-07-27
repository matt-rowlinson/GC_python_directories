#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#SBATCH --job-name=jscale_maps
#SBATCH --ntasks=1
#SBATCH --mem=16Gb
#SBATCH --partition=test
#SBATCH --time=00:30:00
#SBATCH --output=Logs/jscale_map_plots_%A.log
import os
import xarray as xr
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('agg')
import sys
import pandas as pd
import numpy as np
import warnings
warnings.filterwarnings("ignore")

def find_file_list(path, substrs):
    file_list =[]
    for root, directory, files in os.walk(path):
        for f in files:
            for s in substrs:
                if s in f:
                    file_list.append(os.path.join(root, f))
    file_list.sort()
    return file_list

def get_data_as_xr(rundir, year='', collection='SpeciesConc', version='14.1.0'):
    path=f'/mnt/lustre/users/mjr583/GC/{version}/rundirs/{rundir}'
    try:
        ds = xr.open_mfdataset( find_file_list(path, [f'{collection}.{year}']), combine='by_coords' )
    except:
        ds = xr.open_mfdataset( find_file_list(path, [f'{collection}.{year}'])[:-1], combine='by_coords' )
    return ds

def surface_data(rundir, year='', version='14.0.1'):
    add_=[]
    ds0, lat, lon = get_data_as_xr(rundir, year='', version=version)
    ds0 = ds0.mean('time')
    ds0 = ds0.isel( lev=0 ) * 1e12
    return ds0


##----------------------------------------MAIN SCRIPT----------------------------------------##
def main():
    ## Read alts from info_file and plot vertical against that (avg time dim)
    lev=pd.read_csv('/users/mjr583/GC/info_files/GC_72_vertical_levels.csv')['Altitude (km)'] #* 1e3

    ## Makes simple global map plot
    base='new_base_4x5'
    dev ='shah_4x5'
    year = 2015
    
    for rundir in rundirs:
        ds0 = get_data_as_xr(base, collection='SpeciesConc', year=year).mean(dim='time').isel(lev=0)
        ds0 = ds0.where( ds0.lon < 10).where( ds0.lon > -100 )
        ds0 = ds0.where( ds0.lat < 60 ).where( ds0.lat > -30 )
        base_total_NO3 = (ds0['SpeciesConcVV_NIT']   + ds0['SpeciesConcVV_NITs']  + ds0['SpeciesConcVV_NITD1'] +\
                         ds0['SpeciesConcVV_NITD2'] + ds0['SpeciesConcVV_NITD3'] + ds0['SpeciesConcVV_NITD4']) * 1e9
        
        ds1 = get_data_as_xr(rundir, collection='SpeciesConc', year=year).mean(dim='time').isel(lev=0)
        ds1 = ds1.where( ds1.lon < 10).where( ds1.lon > -100 )
        ds1 = ds1.where( ds1.lat < 60 ).where( ds1.lat > -30 )
        dev_total_NO3 =  (ds1['SpeciesConcVV_NIT']   + ds1['SpeciesConcVV_NITs']  + ds1['SpeciesConcVV_NITD1'] +\
                         ds1['SpeciesConcVV_NITD2'] + ds1['SpeciesConcVV_NITD3'] + ds1['SpeciesConcVV_NITD4']) * 1e9
        
        ## Check values in specific gridbox
        x,y = -24.9, 16.9
        abslat = np.abs(ds0.lat-y)
        abslon = np.abs(ds0.lon-x)
        c = np.maximum(abslon, abslat)
        ([xloc], [yloc]) = np.where(c == np.min(c))
        b0 = base_total_NO3.isel(lon=xloc+1, lat=yloc)
        b1 = dev_total_NO3.isel(lon=xloc+1, lat=yloc)
        print( b0.values, b1.values )

        ### Plotting
        fig = plt.figure(figsize=(8,8))
        from matplotlib.colors import LogNorm
        
        Max = 1.01 * np.max( [base_total_NO3.max(), dev_total_NO3.max() ])
        ax1 = plt.subplot(2,2,1, projection=ccrs.PlateCarree(central_longitude=0))
        im = base_total_NO3.plot.imshow(x='lon',y='lat',ax=ax1,
                norm=LogNorm(vmax=Max), transform=ccrs.PlateCarree(), 
                cbar_kwargs={'orientation':'horizontal', 'label':'Base: Total NO3 (ppb)'})

        ax2 = plt.subplot(2,2,2, projection=ccrs.PlateCarree(central_longitude=0))
        im = base_total_NO3.plot.imshow(x='lon',y='lat',ax=ax2, 
                norm=LogNorm(vmax=Max), transform=ccrs.PlateCarree(), 
                cbar_kwargs={'orientation':'horizontal', 'label':f'{rundir.title()}: Total NO3 (ppb)'})

        ax3 = plt.subplot(2,2,3, projection=ccrs.PlateCarree(central_longitude=0))
        delta = dev_total_NO3 - base_total_NO3  
        im = delta.plot.imshow(x='lon',y='lat',ax=ax3, transform=ccrs.PlateCarree(), 
                cbar_kwargs={'orientation':'horizontal', 'label':f'{rundirs} - {base} (ppb)'})

        ax4 = plt.subplot(2,2,4, projection=ccrs.PlateCarree(central_longitude=0))
        delta = (dev_total_NO3 - base_total_NO3)  / base_total_NO3 * 100
        im = delta.plot.imshow(x='lon',y='lat',ax=ax4,transform=ccrs.PlateCarree(), 
                cbar_kwargs={'orientation':'horizontal', 'label':r'$\Delta$%'})
        
        for ax in [ax1,ax2,ax3,ax4]:
            ax.coastlines()
            ax.set_title(None)
            ax.set_extent([-100,10,-30,60])

        plt.tight_layout()
        plt.savefig( f'plots/diff_pNO3_{rundir}.png' , dpi=300)
        plt.close()
        
if __name__ == "__main__":
    main()
