#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#SBATCH --job-name=hemco_timeseries
#SBATCH --ntasks=1
#SBATCH --mem=1gb
#SBATCH --partition=nodes
#SBATCH --time=00:00:10
#SBATCH --output=Logs/hemco_timeseries_%A.log
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('agg')
import os
import numpy as np
import sys
sys.path.append('/users/mjr583/python_lib/')
import RowPy as rp
plt.style.use('seaborn-darkgrid')

def find_file_list(path, substrs):
    file_list =[]
    print( path, substrs )
    for root, directory, files in os.walk(path):
        for f in files:
            for s in substrs:
                if s in f:
                    file_list.append(os.path.join(root, f))
    file_list.sort()
    return file_list

regions={
    False : False,
    'Asia' : { 'minlon' : 70, 'maxlon' : 165, 'minlat' : 0., 'maxlat' : 50. },
    'US'   : { 'minlon' : -160, 'maxlon' : -60, 'minlat' : 15., 'maxlat' : 72. },
    'EU'   : { 'minlon' : -20, 'maxlon' : 40, 'minlat' : 36., 'maxlat' : 72. },
    'NH'   : { 'minlon' : -180, 'maxlon' : 180, 'minlat' : 0., 'maxlat' : 90. },
    'SH'   : { 'minlon' : -180, 'maxlon' : 180, 'minlat' : -90., 'maxlat' : 0. },
        }

def read_hemco(rundirs, versions):
    ds_=[]
    for n, rundir in enumerate(rundirs):
        path=f'/users/mjr583/scratch/GC/{versions[n]}/rundirs/{rundir}/OutputDir/'
        flist=find_file_list(path, ['HEMCO_diagnostics.2016'])
        ds = xr.open_mfdataset(flist)
        ds_.append( ds )
    return ds_

def var_and_convert(ds, species, AREA):
    days=np.array([31,28,31,30,31,30,31,31,30,31,30,31])
    #ds = ds.sum(dim='lev')
    ds = ds[f'Emis{species}_Anthro'].sum(dim='lev') * AREA * (3600 * 24 * days.mean()) * 1e-9
    #ds = ds[f'Emis{species}_Geo'] * AREA * (3600 * 24 * days.mean()) * 1e-9

    return ds

def regional_ems(ds, region, regions):
    xy = list(regions[region].values())
    ds = ds.where( ds.lon > xy[0] ).where( ds.lon < xy[1] ).sum(dim='lon').\
            where( ds.lat > xy[2] ).where( ds.lat < xy[3] ).sum(dim='lat')
    return ds

def main():
    rundirs  = ['ceds_only', 'geo_only', 'asia-meic_scale', 'all_scaled']
    versions = ['13.1.2', '13.4.0','13.1.2', '13.4.0']
    species=['C2H6']

    ds_ = read_hemco( rundirs, versions )
    AREA = ds_[0].AREA.mean(dim='time')

    for s in species:
        print( s )
        fig, ax = plt.subplots(figsize=(10,6))
        for n, d_s in enumerate( ds_ ):
            ds = var_and_convert(d_s, s, AREA)
            #ds = ds.resample(time='Y').sum()
             
            globe = ds.sum(dim='lat').sum(dim='lon')
            print( globe.values.sum(), 'Tg' )
            europ = regional_ems(ds, 'EU',  regions )
            namer = regional_ems(ds, 'US',  regions )
            easia = regional_ems(ds, 'Asia',regions )

            ax.plot(  ds.time, globe, label=versions[n]   )
            #ax.plot(        ds.time, europ, c='grey', label='Global'              )
            #ax.plot(        ds.time, europ, c='blue', label=f'v{versions[n]} Europe'              )
            #ax.plot(        ds.time, namer, c='red',  label=f'v{versions[n]} N. America'          )
            #ax.plot(        ds.time, easia, c='y',    label=f'v{versions[n]} E. Asia'             )
        
        ax.legend()
        ax.set_ylabel( f'{s} (Tg yr$^{-1}$)' )

        plt.savefig( f'plots/TEST.HEMCO_timeseries_{s}.png')
        plt.close()



if __name__ == "__main__":
    main()
