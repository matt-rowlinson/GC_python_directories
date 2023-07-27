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
import numpy as np

import warnings
warnings.filterwarnings("ignore")

def get_cosweights(lat):
    cosweights = np.cos(np.deg2rad(lat))
    cosweights = np.swapaxes((np.tile(cosweights,72)).reshape(72,46),0,1)
    return cosweights

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
    rundirs = ['shah_4x5', 'andersen_4x5', 'langmuir_4x5', 'langmuir_RH_4x5',
                'langmuir_RH_90_4x5','langmuir_RH_95_4x5','langmuir_RH_90_10_4x5',]
    
    for n, rundir in enumerate(rundirs):
        print( rundir )
        Jvalues = get_data_as_xr(rundir, collection='JValues', year='2015')
        Met     = get_data_as_xr(rundirs[0], collection='StateMet', year='2015')
        AirMass = Met['Met_AD'].mean( dim='time' )


        jHNO3 = Jvalues['Jval_HNO3'].mean(dim='time')
        jNIT  = Jvalues['Jval_NIT'].mean(dim='time')
        EF    = jNIT / jHNO3

        #cosw = get_cosweights( jHNO3.lat )
        #print( cosw.max(), cosw.min() )
        #print( cosw.shape )
        
        surf = EF.isel(lev=0) #* cosw
        trop = EF.mean(dim='lon')[:35]
        
        xy = [-24.9,16.9 ]
        abslat = np.abs(EF.lat-xy[1])
        abslon = np.abs(EF.lon-xy[0])
        c = np.maximum(abslon, abslat)
        ([xloc], [yloc]) = np.where(c == np.min(c))
        cvao = EF.isel(lon=xloc, lat=yloc)


        #print( surf.shape )
        print( surf.mean().values, 'surface EF'  )
        print( trop.mean().values, 'trop EF'  )
        print( cvao.mean().values, 'CVAO EF\n'  )



        
        
if __name__ == "__main__":
    main()
