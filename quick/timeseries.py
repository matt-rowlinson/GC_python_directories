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
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('agg')
import numpy as np
plt.style.use('seaborn-darkgrid') 
## temp
import sys
def find_file_list(path, substrs):
    file_list =[]
    for root, directory, files in os.walk(path):
        for f in files:
            for s in substrs:
                if s in f:
                    file_list.append(os.path.join(root, f))
    file_list.sort()
    return file_list

def get_data_as_xr(rundir, year='2017', lev=0, version='13.1.2', variable=False):
    path=f'/mnt/lustre/users/mjr583/GC/{version}/rundirs/{rundir}'
    if variable=='OH':
        ds = xr.open_mfdataset( find_file_list(path, [f'ConcAfterChem.'])[:], combine='by_coords' )
        ds = ds[f'OHconcAfterChem'].isel(lev=lev) 
    else:
        try:
            ds = xr.open_mfdataset( find_file_list(path, [f'SpeciesConc.'])[:], combine='by_coords' )
        except:
            ds = xr.open_mfdataset( find_file_list(path, [f'SpeciesConc.'])[:-1], combine='by_coords' )

    lat = ds['lat']
    lon = ds['lon']
    return ds, lat, lon

def surface_data(rundirs, xy, version='14.0.1'):
    add_=[]
    for rundir in rundirs:
        ds, lat, lon = get_data_as_xr(rundir, version=version)
        ds = ds.isel( lev=0 )

        abslat = np.abs(ds.lat-xy[1])
        abslon = np.abs(ds.lon-xy[0])
        c = np.maximum(abslon, abslat)
        ([xloc], [yloc]) = np.where(c == np.min(c))
        ds = ds.isel(lon=xloc, lat=yloc)
        add_.append( ds )
    return add_

def time_plot(species, xy, ds, rundirs, labels=None, sname=''):
    fig=plt.figure()#figsize=(6,6))
    ax = fig.add_subplot(111)#, projection=ccrs.EqualEarth(), aspect='auto')

    for d, rundir, label in zip(ds, rundirs,labels):

        d = d[f'SpeciesConcVV_{species}'] * 1e12
        d.plot( ax=ax, label=label )
    plt.legend()
    plt.title(f'{xy[1]}$^\circ$N, {xy[0]}$^\circ$E')
    plt.tight_layout()
    plt.savefig( f'plots/{sname}.png' , dpi=300)
    plt.close()
    return

##----------------------------------------MAIN SCRIPT----------------------------------------##
def main():
    ## Specify the GC version, rundir and species to get surface plots of
    rundirs  = ['new_base_4x5','andersen_4x5']
    labels   = ['Base','NO3/RH langmuir']

    version = '14.1.0'
    species_=['HNO2']

    xy=[-24.9,16.9]
    
    ## Makes simple timeseries plot 
    ds_ = surface_data(rundirs, xy, version)
    for species in species_:
        time_plot( species, xy, ds_, rundirs, labels=labels, sname='TEST.14.1.0_NIThv' )
        
if __name__ == "__main__":
    main()
