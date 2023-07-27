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
import pandas as pd
plt.style.use('seaborn-darkgrid')

def find_file_list(path, substrs):
    file_list =[]
    for root, directory, files in os.walk(path):
        for f in files:
            for s in substrs:
                if s in f:
                    file_list.append(os.path.join(root, f))
    file_list.sort()
    return file_list

def read_hemco_output( rundirs, versions, v, sector='Total', years='2017', lonlat=False):
    timeseries=[] 
    for n, rundir in enumerate(rundirs):
        print( rundir )
        if len( versions ) == 1:
            version=versions[0]
        else:
            version=versions[n]
        ds, time = get_hemco_as_xr(rundir, year=years, variable=v, sector=sector, version=version, lonlat=lonlat)
        ds = pd.DataFrame({rundir:ds}, index=time)#.resample('M').mean()
        timeseries.append( ds )
    return timeseries

def get_hemco_as_xr(rundir, year='', variable=False, sector='Total',  version='13.1.2', lonlat=False):
    days=[31,28,31,30,31,30,31,31,30,31,30,31]
    path=f'/mnt/lustre/users/mjr583/GC/{version}/rundirs/{rundir}'
    ds = xr.open_mfdataset( find_file_list(path, [f'HEMCO_diagnostics.{year}']), combine='by_coords' )
    AREA = ds.AREA
    time = ds.time
    
    ds = ds[f'Emis{variable}_{sector}'].sum(dim='lev')
    if lonlat:
        ds = by_region(ds, lonlat)

    ds_=[]
    for dt in ds.time:
        n   =  pd.Timestamp(dt.values).month
        em = ds.sel(time=dt) * AREA.sel(time=dt) * (3600 * 24 * days[n-1]) * 1e-9 
        ds_.append( np.nansum(em.values) )
    return ds_, time

def by_region(ds, lonlat):
    ds = ds.where( ds.lon >= lonlat[0] ).where( ds.lon <= lonlat[1] )\
           .where( ds.lat >= lonlat[2] ).where( ds.lat <= lonlat[3] )
    return ds

##----------------------------------------MAIN SCRIPT----------------------------------------##
def main():
    region=False
    rundirs  = ['ceds_only','asia-meic_scale']#base_run','geo','scale_only','scale_all']
    versions = ['13.1.2']
    
    rundirs   = ['base_run','geo','scale_only','scale_all']
    versions = ['14.0.1']

    year     = '2016'
    variable = 'C2H6'
    sector   = 'Anthro'
    lonlats=[ [-180, 180, 0, 30.],
              [-180, 180, 30., 60.],
              [-180, 180, 60., 90.],
              [-180, 180, 0. , 90.],
              [-180, 180, -90., 0.],
              [-139.95, -50.05, 20.05, 59.95],
              [70.0, 150., 10., 60],
              [-29.95, 89.95, 50.05, 89.95],
              [-29.95, 73.05, 30.05, 50.05] ]
    regions=['0-30','30-60','60-90','NH','SH', 'USA','China','EU1','EU2']

    lonslats=[[66., 97., 5., 33.]]
    regions = ['India']
    
    for lonlat, region in zip( lonlats, regions ):
        ds = read_hemco_output( rundirs, versions, variable, sector, years=year, lonlat=lonlat)
        
        f, ax = plt.subplots(figsize=(10,4))
        styles = ['-','--','-','--']
        labels = ['Base','Geo only','Scale only','Scale and Geo'] ; labs=[]
        for df, style, label in zip(ds, styles, labels):
            df.plot(drawstyle='steps-mid', linestyle=style, ax=ax)
            labs.append(f'{label} ({df.values.sum().round(2)} Tg yr-1)')
        
        ax.set_ylabel( f'{sector} emissions of {variable} (Tg)' )
        ax.legend(labs)
        if region:
            plt.savefig( f'plots/HEMCO.{region}.{sector}_{variable}.png')
        else:
            plt.savefig( f'plots/HEMCO.{sector}_{variable}.png')
        plt.close()    

if __name__ == "__main__":
    main()
