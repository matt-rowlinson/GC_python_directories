#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#SBATCH --job-name=diff_maps
#SBATCH --ntasks=1
#SBATCH --mem=16Gb
#SBATCH --partition=test
#SBATCH --time=00:30:00
#SBATCH --output=Logs/simple_diff_map_plots_%A.log
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
    base='langmuir_4x5'
    rundirs = ['DD_langmuir_4x5']
    labels = ['Langmuir wDep']
    #rundirs = ['shah_4x5','andersen_4x5','langmuir_4x5','langmuir_RH_4x5',
    #        'langmuir_RH_90_4x5','langmuir_RH_95_4x5','langmuir_RH_90_10_4x5','langmuir_RH_95_10_4x5']
    #labels=['Shah et al.','Andersen et al.','New langmuir', r'New langmuir $\times$ RH',
    #            r'Langmuir $\times$ RH, UL=90th',r'Langmuir $\times$ RH, UL=95th',
    #            r'Langmuir $\times$ RH, LL=10th, UL=90th', r'Langmuir $\times$ RH, LL=10th, UL=95th'] 
    year = 2015
    
    for rundir,label in zip(rundirs,labels):
        ds0 = get_data_as_xr(base, collection='SpeciesConc', year=year).mean(dim='time').isel(lev=0)
        base_total_NO3 = (ds0['SpeciesConcVV_NIT']   + ds0['SpeciesConcVV_NITs']  + ds0['SpeciesConcVV_NITD1'] +\
                         ds0['SpeciesConcVV_NITD2'] + ds0['SpeciesConcVV_NITD3'] + ds0['SpeciesConcVV_NITD4']) * 1e9
        
        ds1 = get_data_as_xr(rundir, collection='SpeciesConc', year=year).mean(dim='time').isel(lev=0)
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
                cbar_kwargs={'orientation':'horizontal', 'label':f'{label}: Total NO3 (ppb)'})

        ax3 = plt.subplot(2,2,3, projection=ccrs.PlateCarree(central_longitude=0))
        delta = dev_total_NO3 - base_total_NO3  
        im = delta.plot.imshow(x='lon',y='lat',ax=ax3, transform=ccrs.PlateCarree(), 
                cbar_kwargs={'orientation':'horizontal', 'label':f'{label} - {base} (ppb)'})

        ax4 = plt.subplot(2,2,4, projection=ccrs.PlateCarree(central_longitude=0))
        delta = (dev_total_NO3 - base_total_NO3)  / base_total_NO3 * 100
        im = delta.plot.imshow(x='lon',y='lat',ax=ax4,transform=ccrs.PlateCarree(), 
                cbar_kwargs={'orientation':'horizontal', 'label':r'$\Delta$%'})
        
        for ax in [ax1,ax2,ax3,ax4]:
            ax.coastlines()
            ax.set_title(None)

        plt.tight_layout()
        plt.savefig( f'plots/diff_NIThvOn_pNO3_{rundir}.png' , dpi=300)
        plt.close()
        
if __name__ == "__main__":
    main()
