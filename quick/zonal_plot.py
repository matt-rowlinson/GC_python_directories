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
import pandas as pd
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

def get_data_as_xr(rundir, year='2017', lev=0, version='13.1.2', variable=False):
    path=f'/mnt/lustre/users/mjr583/GC/{version}/rundirs/{rundir}'
    if variable=='OH':
        ds = xr.open_mfdataset( find_file_list(path, [f'ConcAfterChem.'])[:], combine='by_coords' )
        ds = ds[f'OHconcAfterChem'].isel(lev=lev) 
    else:
        try:
            ds = xr.open_mfdataset( find_file_list(path, [f'SpeciesConc.{year}'])[:], combine='by_coords' )
        except:
            ds = xr.open_mfdataset( find_file_list(path, [f'SpeciesConc.{year}'])[:-1], combine='by_coords' )
    return ds 

def surface_data(rundir, version='14.0.1'):
    add_=[]
    ds0 = get_data_as_xr(rundir, year='', version=version)
    ds0 = ds0.mean('time')
    ds0 = ds0.mean('lon')
    return ds0

def cart_plot(species, p, base, lim=200, sname=''):
    i = min(range(len(p)), key=lambda i: abs(p[i]-lim))

    fig=plt.figure()
    ax = fig.add_subplot(111)
    im = ax.pcolormesh( base['lat'], p[:i], base[:i,:] )
    plt.xlabel('Latitude') ; plt.ylabel('Pressure (hPa)')   
    cb = plt.colorbar(im, orientation='horizontal')
    cb.set_label( f'{species} mol mol-1' )
    plt.tight_layout()
    plt.gca().invert_yaxis()
    plt.savefig( f'plots/{sname}.png' , dpi=300)
    plt.close()
    return

##----------------------------------------MAIN SCRIPT----------------------------------------##
def main():
    ## Specify the GC version, rundir and species to get surface plots of
    rundir  = 'base_run'
    version = '14.0.1'
    species_=['OH','CH4']

    lim = 20 ## Upper pressure limit in hPa 
    
    ## Makes simple zonal latitude plot 
    base_ = surface_data(rundir, version)
    gc_p = pd.read_csv( '/users/mjr583/GC/info_files/GC_72_vertical_levels.csv')['Pressure (hPa)']
    gc_alts = pd.read_csv( '/users/mjr583/GC/info_files/GC_72_vertical_levels.csv')['Altitude (km)']
    print( gc_p )
    for species in species_:
        print( species )
        base = base_[f'SpeciesConc_{species}']
        cart_plot( species, gc_p, base, lim=lim, sname=f'v{version}.{rundir}.{species}_zonal')
        
if __name__ == "__main__":
    main()
