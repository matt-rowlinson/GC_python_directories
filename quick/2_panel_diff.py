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

def find_file_list(path, substrs):
    file_list =[]
    for root, directory, files in os.walk(path):
        for f in files:
            for s in substrs:
                if s in f:
                    file_list.append(os.path.join(root, f))
    file_list.sort()
    return file_list

def get_data_as_xr(rundir, year='', lev=0, version='13.1.2', variable=False):
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

def surface_data(rundir, year='', version='14.0.1'):
    add_=[]
    ds0, lat, lon = get_data_as_xr(rundir, year='', version=version)
    ds0 = ds0.mean('time')
    ds0 = ds0.isel( lev=0 ) * 1e12
    return ds0

def cart_plot(ds, sname=''):
    fig=plt.figure()#figsize=(6,6))
    titles=['Base - Langmuir', 'Base / Langmuir (%)']
    for n in range(len( ds )):
        ax = fig.add_subplot(2,1,n+1, projection=ccrs.EqualEarth(), aspect='auto')
        im = ds[n].plot.imshow( x='lon',y='lat', ax=ax, transform=ccrs.PlateCarree(), cbar_kwargs={'orientation':'vertical'}) 
        ax.coastlines()
        ax.set_title(titles[n] ) 
    plt.tight_layout()
    plt.savefig( f'plots/{sname}.png' , dpi=300)
    plt.close()
    return

##----------------------------------------MAIN SCRIPT----------------------------------------##
def main():
    ## Specify the GC version, rundir and species to get surface plots of
    version = '14.1.0'
    species_=['HNO2']#SALA','SALC','Cl','SALCCL','SALACL']
    
    ## Makes simple glboal map plot
    import sys
    base_ = surface_data( 'new_base_4x5', '2015', version )
    dev_  = surface_data( 'DD_new_base_4x5', '2015', version )
    for species in species_:
        print( species )
        base = base_[f'SpeciesConcVV_{species}']
        dev  = dev_[f'SpeciesConcVV_{species}']
        cart_plot( [dev-base,- (100 - (dev/base*100))], sname=f'DD_Base.{species}_surfacemap')
        
if __name__ == "__main__":
    main()
