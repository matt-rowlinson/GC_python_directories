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

def get_data_as_xr(rundir, lev=0, dates='2017', version='13.1.2', variable=False):
    path=f'/mnt/lustre/users/mjr583/GC/{version}/rundirs/{rundir}'
    if variable=='OH':
        ds = xr.open_mfdataset( find_file_list(path, [f'ConcAfterChem.'])[:], combine='by_coords' )
        ds = ds[f'OHconcAfterChem'].isel(lev=lev) 
    else:
        try:
            ds = xr.open_mfdataset( find_file_list(path, [f'SpeciesConc.{dates}'])[:], combine='by_coords' )
        except:
            ds = xr.open_mfdataset( find_file_list(path, [f'SpeciesConc.{dates}'])[:-1], combine='by_coords' )

    lat = ds['lat']
    lon = ds['lon']
    return ds, lat, lon

def surface_data(rundirs, xy, version='14.0.1', dates='2017'):
    add_=[]
    for rundir in rundirs:
        ds, lat, lon = get_data_as_xr(rundir, version=version, dates=dates)
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
        print( rundir )
        if species=='nitrate':
            d = (d['SpeciesConcVV_NIT']   + d['SpeciesConcVV_NITs']   + d['SpeciesConcVV_NITD1'] +
                 d['SpeciesConcVV_NITD2'] + d['SpeciesConcVV_NITD3'] + d['SpeciesConcVV_NITD4']) * 1e9
        else:
            d = d[f'SpeciesConcVV_{species}'] * 1e12

        d = d.resample(time='M').mean()
        #sys.exit()
        d.plot( ax=ax, label=label )
    #if species == 'HNO2':
    #    ax.fill_between(d.time, 3.,5.,color='r', alpha=0.2, lw=0, label='Andersen et al. 2022')
    ax.set_ylabel( "HONO (ppt)" )
    plt.legend()
    plt.title(f'{xy[1]}$^\circ$N, {xy[0]}$^\circ$E')
    plt.tight_layout()
    plt.savefig( f'plots/{sname}.png' , dpi=300)
    plt.close()
    return

##----------------------------------------MAIN SCRIPT----------------------------------------##
def main():
    ## Specify the GC version, rundir and species to get surface plots of
    rundirs  = ['new_base_4x5','shah_4x5']#,'andersen_4x5','langmuir_4x5','langmuir_RH_4x5']
    labels=['Base','Shah et al.']#,'Andersen et al.','New langmuir', r'New langmuir $\times$ RH']#Andersen_RHbin','Shah et al.']

    version = '14.1.0'
    species_=['HNO2']#HNO2', 'O3', 'CO']
    dates   ='2015'

    xy=[-24.9,16.9]   ## Cape Verde
    xy=[-64.75,32.3]## Bermuda

    
    ## Makes simple timeseries plot 
    ds_ = surface_data(rundirs, xy, version, dates)
    for species in species_:
        print( species )
        time_plot( species, xy, ds_, rundirs, labels=labels, sname=f'bermuda.timeseries.14.1.0.{species}_NIThv' )
        
if __name__ == "__main__":
    main()
