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

def burden_data(rundirs, species, versions=['13.1.2']):
    spec_db = yaml.load(open("./species_database.yml"), Loader=yaml.FullLoader)
    Met, lat, lon = get_data_as_xr('base_run_2000', collection='StateMet', year='')

    add_=[] ; y_=[]
    for n, rundir in enumerate(rundirs):
        print( rundir, rundir[-4:] )
        year=rundir[-4:]

        ds, lat, lon = get_data_as_xr(rundir, year=year)

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
            if str(y) == year:
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


def bar_plot(df, save_name='TEST'):
    f,ax=plt.subplots(figsize=(8,4))
    c = ['grey','#e41a1c','#377eb8','#984ea3','#ff7f00']
    ax.bar( df.index, np.squeeze(df.values), width=.7, color=c )
    ax.set_ylabel('Tropospheric O$_3$ burden (Tg)')
    ax.set_ylim(top=378)
    plt.tight_layout()
    plt.savefig(f'plots/{save_name}.png')
    plt.close()
    return 

##----------------------------------------MAIN SCRIPT----------------------------------------##
def main():
    species='O3'
    
    # Burden plot for present-day (2015)
    rundirs     = ['base_run_2015', 'j25_run_2015','viral_run_2015', 'Ander22b_run_2015']
    run_nam     = ['Base', 'j25 (Kasibhatla et al. 2018)','Shah et al. 2022', 'Andersen et al. 2022']
    PD, years   = burden_data( rundirs, species, versions=['13.1.2'] )
    PD          = pd.DataFrame({species : PD }, index=run_nam).astype(float)
    bar_plot(PD, 'PD_tropospheric_o3_burdens.png')

    # Burden plot for preindustrial (1750)
    rundirs     = ['base_run_1751', 'base_run_1751','viral_run_1751', 'Ander22b_run_1751']
    run_nam     = ['Base', 'j25 (Kasibhatla et al. 2018)','Shah et al. 2022', 'Andersen et al. 2022']
    PI, years   = burden_data( rundirs, species, versions=['13.1.2'] )
    PI          = pd.DataFrame({species : PI }, index=run_nam).astype(float)
    bar_plot(PI, 'PI_tropospheric_o3_burdens.png')


        
if __name__ == "__main__":
    main()
