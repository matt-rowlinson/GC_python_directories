#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#SBATCH --job-name=o3bars
#SBATCH --ntasks=1
#SBATCH --mem=156Gb
#SBATCH --partition=himem
#SBATCH --time=00:30:00
#SBATCH --output=Logs/o3_bars_%A.log
import warnings
warnings.simplefilter(action='ignore', category=RuntimeWarning)
import sys
import os
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('agg')
plt.style.use('seaborn-darkgrid')
import numpy as np
import pandas as pd
import yaml

def find_file_list(path, substrs):
    file_list =[]
    for root, directory, files in os.walk(path):
        for f in files:
            for s in substrs:
                if s in f:
                    file_list.append(os.path.join(root, f))
    file_list.sort()
    return file_list

def get_data_as_xr(rundir, year='', lev=0, collection='SpeciesConc', version='13.1.2', variable=False):
    path=f'/mnt/lustre/users/mjr583/GC/{version}/rundirs/{rundir}'
    print( path, collection, year )
    if variable=='OH':
        ds = xr.open_mfdataset( find_file_list(path, [f'ConcAfterChem.{year}'])[:], combine='by_coords' )
        ds = ds[f'OHconcAfterChem'].isel(lev=lev) 
    else:
        try:
            ds = xr.open_mfdataset( find_file_list(path, [f'{collection}.{year}'])[:], combine='by_coords' )
        except:
            ds = xr.open_mfdataset( find_file_list(path, [f'{collection}.{year}'])[:-1], combine='by_coords' )

    lat = ds['lat']
    lon = ds['lon']
    return ds, lat, lon

def mean_surface_data(rundirs, versions=['13.1.2']):
    add_=[]
    for n, rundir in enumerate(rundirs):
        print( rundir )
        ds0, lat, lon = get_data_as_xr(rundir, year='', version=versions[0])
        ds0 = ds0.mean('lon').mean('lat')
        ds0 = ds0.isel( lev=0, drop=True )
        add_.append( ds0['SpeciesConc_O3'] * 1e9 )
    return add_

def burden_data(rundirs, species, versions=['13.1.2']):
    spec_db = yaml.load(open("./species_database.yml"), Loader=yaml.FullLoader)

    add_=[] ; y_=[]
    for n, rundir in enumerate(rundirs):
        print( rundir, rundir[-4:] )
        if "scale-25" in rundir:
            year='2015'
        else:
            year=rundir[-4:]

        ds, lat, lon = get_data_as_xr(rundir, year=year)
        CaC, lat, lon = get_data_as_xr('base_run_2000', collection='ConcAfterChem', year='')
        PrLo, lat, lon = get_data_as_xr('base_run_2000', collection='ProdLoss', year='')
        Met, lat, lon = get_data_as_xr('base_run_2000', collection='StateMet', year='')

        trop = Met['Met_TropP'].mean( dim='time' )
        AirMass = Met['Met_AD'].mean( dim='time' )
        AirDen = Met['Met_AIRDEN'].mean( dim='time' )
        pmid_press = Met['Met_PMID'].mean( dim='time' )
        MASK =  ( pmid_press > trop )
        
        if species=='NIT' or species =='SO4':
            var = ds[f'SpeciesConc_{species}'] + ds[f'SpeciesConc_{species}s'] + ds[f'SpeciesConc_{species}D1'] \
                  + ds[f'SpeciesConc_{species}D2'] + ds[f'SpeciesConc_{species}D3'] + ds[f'SpeciesConc_{species}D4']
        else:
            var = ds[f'SpeciesConc_{species}']
        time = ds['time']

        years = pd.DatetimeIndex( time.values )
        years = years.year.drop_duplicates().tolist()[:]
        for y in years:
            print( y )
            v = var.sel(time=str(y))
            v = v.mean(dim='time')
            v = v.where(MASK)
            RMM_air = (.78*(2.*14.)+.22*(2.*16.))
            conversion_factor = (AirMass*1E3 / RMM_air)
            v = v * conversion_factor
            v = v * float(spec_db[species]['MW_g'] * 1e-3 ) / 1E9
            
            add_.append( v.sum().values ) 
            y_.append( y )


    return add_, y_

##----------------------------------------MAIN SCRIPT----------------------------------------##
def main():
    global d, site, rundirs, v
    versions = ['13.1.2']
    interval = [10, 10, 5, 5, 5, 5, 5]
    v='$O_3$'

    rundirs  = ['base_run_1751','base_run_1980', 'base_run_1990', 'base_run_2000', 'base_run_2005',
                                'base_run_2010', 'base_run_2015', 'base_run_2020']   
    rundirs = ['nitrate_photol_2015_all_scale-25']
    species='O3'
    j100_, years = burden_data( rundirs, species, versions )

    df = pd.DataFrame({species : j100_}, index=years)
    print( df )
    
    savedf =  df#.resample( 'Y' ).mean()
    savedf.to_csv(f'csv/2015_j25_burden_{species}.csv') 
        
if __name__ == "__main__":
    main()
