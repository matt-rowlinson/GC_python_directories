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
    rundirs = ['shah_4x5','andersen_4x5','langmuir_4x5','langmuir_RH_4x5','langmuir_RH_90_4x5','langmuir_RH_95_4x5','langmuir_RH_90_10_4x5']
    labels=['Shah et al.','Andersen et al.','New langmuir', r'New langmuir $\times$ RH',
                r'Langmuir $\times$ RH, max$_{90th}$',r'Langmuir $\times$ RH, max$_{95th}$',
                r'Langmuir $\times$ RH, min$_{10th}$ max$_{90th}$']
    
    for n, rundir in enumerate(rundirs):
        print( '___',rundir,'___' )
        Jvalues     = get_data_as_xr(rundir, collection='JValues', year='2015')

        jHNO3 = Jvalues['Jval_HNO3'].mean(dim='time')
        jNIT  = Jvalues['Jval_NIT'].mean(dim='time')
        EF    = jNIT / jHNO3
    
        fig = plt.figure(figsize=(10,6))
        
        ax1 = plt.subplot(1,2,1, projection=ccrs.PlateCarree(central_longitude=0))
        im = EF.isel(lev=0).plot.imshow(x='lon',y='lat',ax=ax1, transform=ccrs.PlateCarree(), 
                cbar_kwargs={'orientation':'horizontal', 'label':'Enhancement Factor'})
        ax1.coastlines()
        ax1.set_title(None) 
        ax1.text( -120,100, labels[n], fontsize=16, weight='bold' )
        
        ax2 = plt.subplot(1,2,2)
        print( 'Min: ' , EF.min( ).values )
        print( 'Max: ' , EF.max( ).values )
        print( 'Mean: ', EF.mean().values,'\n' )

        ax2.pcolormesh( EF.lat, lev[:35], EF.mean(dim='lon')[:35] )
        ax2.set_ylabel( 'Altitude (km)' )
            
        plt.tight_layout()
        plt.savefig( f'plots/EF.{rundir}.png' , dpi=300)
        plt.close()
    
if __name__ == "__main__":
    main()
