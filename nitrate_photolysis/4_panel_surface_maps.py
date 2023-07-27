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

def get_data_as_xr(rundir, year='20140101', lev=0, version='13.1.2', variable=False):
    path=f'/mnt/lustre/users/mjr583/GC/{version}/rundirs/{rundir}'
    if variable=='OH':
        ds = xr.open_mfdataset( find_file_list(path, [f'ConcAfterChem.'])[:], combine='by_coords' )
        ds = ds[f'OHconcAfterChem'].isel(lev=lev) 
    else:
        try:
            ds = xr.open_mfdataset( find_file_list(path, [f'SpeciesConc.{year}'])[:], combine='by_coords' )
        except:
            ds = xr.open_mfdataset( find_file_list(path, [f'SpeciesConc.{year}'])[:-1], combine='by_coords' )

    lat = ds['lat']
    lon = ds['lon']
    return ds, lat, lon

def surface_data(rundir, version='14.0.1', year='2014'):
    add_=[]
    ds0, lat, lon = get_data_as_xr(rundir, year=year, version=version)
    ds0 = ds0.mean('time')
    ds0 = ds0.isel( lev=0 )
    return ds0

def cart_plot(species, base, dev, unit, vmax, sname=''):
    fig=plt.figure(figsize=(10,6))

    z_ = [base, dev, dev-base, 100 - (base/dev*100)]
    for n, z in zip(range(4),z_):
        ax = fig.add_subplot(2,2,n+1, projection=ccrs.EqualEarth(), aspect='auto' )
        if n < 2:
            im = z.plot.imshow( x='lon',y='lat', ax=ax, transform=ccrs.PlateCarree(), 
                    vmax=vmax, cbar_kwargs={'orientation':'horizontal','label':f'SpeciesConc {species} ({unit})'}) 
        else:
            if n==3:
                unit='%'
            im = z.plot.imshow( x='lon',y='lat', ax=ax, transform=ccrs.PlateCarree(), cmap='bwr', 
                    center=0., cbar_kwargs={'orientation':'horizontal','label':f'SpeciesConc {species} ({unit})'}) 
        ax.coastlines()
        ax.set_title(None)
    
    plt.tight_layout()
    plt.savefig( f'plots/{sname}.png' , dpi=300)
    plt.close()
    return

##----------------------------------------MAIN SCRIPT----------------------------------------##
def main():
    ## Specify the GC version, rundir and species to get surface plots of
    rundirs = ['base_4x5', 'new_base_4x5']
    version = '14.1.0'
    species_= ['O3','CO','NOx']#SALA','SALC','Cl','SALCCL','SALACL']
    scale   = 1e9
    unit    = 'ppbv'
    vmax    = None
    ## Makes simple global map plot 
    base_ = surface_data(rundirs[0], version ,year='20140101')
    dev_  = surface_data(rundirs[1], version, year='20140101')

    for s in species_:
        print( s )
        if s == 'NOx':
            base = (base_[f'SpeciesConcVV_NO'] + base_[f'SpeciesConcVV_NO2']  )* scale
            dev  = (dev_[ f'SpeciesConcVV_NO'] + base_[f'SpeciesConcVV_NO2'] )* scale
        else:
            base = base_[f'SpeciesConcVV_{s}'] * scale
            dev  = dev_[ f'SpeciesConcVV_{s}'] * scale
        cart_plot( s, base, dev, unit, vmax, sname=f'TEST.new_base.v{version}.{s}_surfacemap')
        
if __name__ == "__main__":
    main()
