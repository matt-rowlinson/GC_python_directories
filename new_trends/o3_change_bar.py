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
        ds = xr.open_mfdataset( find_file_list(path, [f'ConcAfterChem.{year}'])[:], combine='by_coords' )
        ds = ds[f'OHconcAfterChem'].isel(lev=lev) 
    else:
        try:
            ds = xr.open_mfdataset( find_file_list(path, [f'SpeciesConc.{year}'])[:], combine='by_coords' )
        except:
            ds = xr.open_mfdataset( find_file_list(path, [f'SpeciesConc.{year}'])[:-1], combine='by_coords' )

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
    add_=[]
    for n, rundir in enumerate(rundirs):
        print( rundir )
        #ds0, lat, lon = get_data_as_xr(rundir, year='', version=versions[0])
        #ds0 = ds0.mean('lon').mean('lat')
        #ds0 = ds0.isel( lev=0, drop=True )
        #add_.append( ds0['SpeciesConc_O3'] * 1e9 )

        ds  = get_data_as_xr(rundir, year='')
        CaC = get_data_as_xr(rundirs[0], collection='ConcAfterChem', year='')
        PrLo = get_data_as_xr(rundirs[0], collection='ProdLoss', year='')
        Met = get_data_as_xr(rundirs[0], collection='StateMet', year='')
        trop = Met['Met_TropP'].mean( dim='time' )
        AirMass = Met['Met_AD'].mean( dim='time' )
        AirDen = Met['Met_AIRDEN'].mean( dim='time' )
        pmid_press = Met['Met_PMID'].mean( dim='time' )
        MASK =  ( pmid_press > trop )#.mean( dim='time' )

        var = ds[f'SpeciesConc_{species}']
        var = var.mean(dim='time')

        var = var.where(MASK)
        RMM_air = (.78*(2.*14.)+.22*(2.*16.))
        conversion_factor = (AirMass*1E3 / RMM_air)
        var = var * conversion_factor
        var = var * float(spec_db[Ycore_spec[nn]]['MW_g'] * 1e-3 ) / 1E9
        print( var )
        sys.exit()
        a.append( var.sum().values ) 
        a.append( ds[f'SpeciesConc_{spec}'].mean(dim='time').isel(lev=0).mean().values * 1e9 )


    return add_

def ts_plot(rundirs, base_, sname=''):
    fig, ax = plt.subplots(figsize=(12,6))
    c = ['#edf8b1','#c7e9b4','#7fcdbb','#41b6c4','#1d91c0','#225ea8','#0c2c84']
    c = ['#ece2f0','#d0d1e6','#a6bddb','#67a9cf','#3690c0','#02818a','#016450']
    c = ['#d0d1e6','#a6bddb','#67a9cf','#3690c0','#02818a','#016c59','#014636']
    for n, run in enumerate(rundirs):
        # drop leap year day to plot all years alongside each other
        base_[n] = base_[n].sel(time=~((base_[n].time.dt.month == 2) & (base_[n].time.dt.day == 29)))
        
        z = j100_[n] #/ base_[n] * 100) - 100
        m = str(np.round( z.mean().values, 2))
        ax.plot( base_[0]['time'], z, c=c[n], label=f'{run[-4:]} ({m}%)' )

    import matplotlib.dates as mdates
    myFmt = mdates.DateFormatter('%b')
    ax.xaxis.set_major_formatter(myFmt)

    ax.set_ylabel( '$O_3$ (ppbv)' )
    plt.legend()
    plt.savefig(f'plots/TEST.{sname}.png')
    plt.close()
    return

##----------------------------------------MAIN SCRIPT----------------------------------------##
def main():
    global d, site, rundirs, v
    versions = ['13.1.2']
    interval = [10, 10, 5, 5, 5, 5, 5]
    v='$O_3$'

    yearly = pd.read_csv('csv/1980-2019_all_base_mean_surface_o3.csv', index_col=0)
    yearly.index = pd.to_datetime( yearly.index, format='%Y-%m-%d')

    base = pd.read_csv('csv/1750-2020_base_mean_surface_o3_SPUN.csv', index_col=0).dropna()
    base.index = pd.to_datetime( base.index, format='%Y-%m-%d')

    j100 = pd.read_csv('csv/1750-2020_j100_mean_surface_o3_SPUN.csv', index_col=0).dropna()
    j100.index = pd.to_datetime( j100.index, format='%Y-%m-%d')

    j25 = pd.read_csv('csv/2015_j25_mean_surface_o3.csv', index_col=0).dropna()

    viral = pd.read_csv('csv/2015_viral_mean_surface_o3.csv', index_col=0).dropna()
    viral.index = pd.to_datetime( viral.index, format='%Y-%m-%d') 
    
    a = base.iloc[0]   # 1750
    z = base.iloc[-1]  # 2020
    ja = j100.iloc[0]   # 1750 with j100

    x = yearly.index.year.drop_duplicates().tolist()
    jx = j100.iloc[1:].index.year.drop_duplicates().tolist()
    
    c= '#ff7f00'
    
    f, ax = plt.subplots(figsize=(12,6))
    ax.bar( 1977, a, color=c)
    ax.bar( 2020, z, color=c)
    zz = yearly.copy()

    yearly['O3'].loc['1980'] = np.nan
    yearly['O3'].loc['1990'] = np.nan
    yearly['O3'].loc['2000'] = np.nan
    yearly['O3'].loc['2005'] = np.nan
    yearly['O3'].loc['2010'] = np.nan
    yearly['O3'].loc['2015'] = np.nan
    yearly['O3'].loc['2020'] = np.nan
    ax.bar( x, yearly["O3"], color='slategrey',zorder=3, label='Base')
    
    ax.bar( x, zz['O3'], color=c, zorder=2, label='Base (year run with NIThv)')

    ax.bar( 1977, ja, fill=False, zorder=0, width=.7, ls='--', edgecolor=c,  linewidth=1.6)
    ax.bar( jx, j100['O3'].iloc[1:], width=0.7, ls='--', linewidth=1.6,
                edgecolor=c, fill=False, zorder=0, label='NIThv at J100' )
    ax.bar( 2015, viral['O3'][1], fill=False, zorder=1, width=.6, ls='--', linewidth=2, 
            edgecolor='#984ea3', label='Shah parameterisation')
    ax.bar( 1977, viral['O3'][0], fill=False, zorder=1, width=.6, ls='--', linewidth=2,
            edgecolor='#984ea3', label='_Shah parameterisation')
    ax.bar( 2015, j25.O3[0], fill=False, zorder=1, width=.6, ls='--', linewidth=2,
            edgecolor='r', label='NIThv at j25')

    x.insert(0,1977)
    x.insert(41,2020)
    plt.ylim( top=39.5 )
    labels=['1750','1980','','','','','1985','','','','',
            '1990','','','','','1995','','','','',
            '2000','','','','','2005','','','','',
            '2010','','','','','2015','','','','',
            '2020']
    ax.set_xticks(x,labels)

    plt.legend(ncol=5)
    plt.ylabel( "Global mean surface $O_3$ (ppbv)" )
    plt.savefig(f'plots/o3_mean_bar.png')
    plt.close()
    
        
if __name__ == "__main__":
    main()
