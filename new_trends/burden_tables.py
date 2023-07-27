#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#SBATCH --job-name=table_stats
#SBATCH --ntasks=1
#SBATCH --mem=16Gb
#SBATCH --partition=test
#SBATCH --time=00:30:00
#SBATCH --output=Logs/table_%A.log
#import warnings
#warnings.simplefilter(action='ignore', category=FutureWarning)
import sys
import os
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import pandas as pd
matplotlib.use('agg')
import yaml
from matplotlib.backends.backend_pdf import PdfPages
plt.style.use('seaborn-darkgrid')

def find_file_list(path, substrs):
    file_list =[]
    for root, directory, files in os.walk(path+'/OutputDir/'):
        for f in files:
            for s in substrs:
                if s in f:
                    file_list.append(os.path.join(root, f))
    file_list.sort()
    return file_list

def get_data_as_xr(rundir, year='', collection='SpeciesConc', archive=False, archive_dir='',variable=False, Chem=False, version='13.1.2'):
    if archive:
        path=f'/mnt/lustre/users/mjr583/GC/{version}/rundirs/{rundir}/{archive_dir}'
    else:
        path=f'/mnt/lustre/users/mjr583/GC/{version}/rundirs/{rundir}'

    if variable=='OH' or variable=='HO2':
        collection='ConcAfterChem'
        Chem=True
    try:
        if year=='201001':
            ds = xr.open_mfdataset( find_file_list(path, [f'{collection}.{year}']), combine='by_coords' )
        else:
            ds = xr.open_mfdataset( find_file_list(path, [f'{collection}.{year}']), combine='by_coords' )
    except:
        try:
            ds = xr.open_mfdataset( find_file_list(path, [f'{collection}.{year}']), combine='by_coords' )
        except:
            ds = xr.open_mfdataset( find_file_list(path, [f'{collection}.{year}']), combine='by_coords' )
    return ds


def map_plot(species, r, rundir, n, levels=None):
    plt.figure(figsize=(10,6))
    cbar_kwargs = {'orientation':'horizontal', 'label':f'[{species}]: dev / base','spacing': 'proportional'}#, 'ticks':levels[::10]}
    ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=0))
    im = r.plot.imshow(x='lon', y='lat', ax=ax, cmap='bwr', levels=levels, 
                              transform=ccrs.PlateCarree(), cbar_kwargs=cbar_kwargs
                              )
    ax.coastlines()
    plt.title(f'v13 {rundir} / Base') 
    plt.savefig( f'plots/{rundir}_{species}_pc.png' )
    plt.close()
    return

##----------------------------------------MAIN SCRIPT----------------------------------------##
def main():
    spec_db = yaml.load(open("./species_database.yml"), Loader=yaml.FullLoader)

    core_spec=['O3','CO','NO','NO2']
    Ycore_spec=['O3','CO','NO','NO2']
    rundirs  = ['Ander22b_run_2015']#,'viral_run_2015','Ander22b_run_2015']
    labels=['And22b et al. 2022']#,'Shah et al. 2022','Andersen et al. 2022']
    b=[]
    for n, rundir in enumerate(rundirs[:]):
        print( rundir )
        ds  = get_data_as_xr(rundir, year='2015')
        CaC = get_data_as_xr(rundirs[0], collection='ConcAfterChem', year='2015')
        PrLo = get_data_as_xr(rundirs[0], collection='ProdLoss', year='2015')
        Met = get_data_as_xr(rundirs[0], collection='StateMet', year='2015')
        print( 'a')
        trop = Met['Met_TropP'].mean( dim='time' )
        AirMass = Met['Met_AD'].mean( dim='time' )
        AirDen = Met['Met_AIRDEN'].mean( dim='time' )
        pmid_press = Met['Met_PMID'].mean( dim='time' )
        MASK =  ( pmid_press > trop )#.mean( dim='time' )
        print( 'b' )
        a=[]
         
        for nn, spec in enumerate(core_spec):
            print( spec )
            var = ds[f'SpeciesConc_{spec}']
            var = var.mean(dim='time')

            var = var.where(MASK)
            RMM_air = (.78*(2.*14.)+.22*(2.*16.))
            conversion_factor = (AirMass*1E3 / RMM_air)
            var = var * conversion_factor
            var = var * float(spec_db[Ycore_spec[nn]]['MW_g'] * 1e-3 ) / 1E9

            a.append( var.sum().values ) 
            print( var.sum().values )
            a.append( ds[f'SpeciesConc_{spec}'].mean(dim='time').isel(lev=0).mean().values * 1e9 )
        
        AirDen = Met['Met_AIRDEN'].mean( dim='time' )
        AirVol = Met['Met_AIRVOL'].mean( dim='time' )
        OHmass = CaC['OHconcAfterChem'].mean(dim='time') * AirMass
        print( OHmass )
        OHmass_sum = np.nansum( OHmass ) / np.nansum( AirMass )
        print( OHmass_sum )
        a.append( OHmass_sum * 1e-6 )

        ## ADD Prod-Loss Ox
        '''
        AVO= 6.0221413E23
        loss = PrLo['Loss_Ox'].mean(dim='time')
        prod = PrLo['Prod_Ox'].mean(dim='time')
        loss *= AirVol * 1E6 / AVO
        loss *= spec_db['O3']['MW_g'] / 1E12
        loss = loss*60*60*24*365
        units = 'Tg/year'
        loss = loss.sum()[var].values
        sys.exit()'''
        b.append( a )

    c = np.array( b )
    df = pd.DataFrame( c, index=labels, 
            columns= ['O3 burden (Tg)','O3 surface (ppbv)', 'CO burden (Tg)','CO surface (ppbv)',
                      'NO burden (Tg)','NO surface (ppbv)','NO2 burden (Tg)','NO2 surface (ppbv)',
                      'OH mean (1e6 molec cm-3)'] )
    df = df.transpose().round(2)
    df.to_csv(f"Output/{rundirs}_Limits_stats.csv")

    fig, ax =plt.subplots(figsize=(7.5,2.5))
    ax.axis('tight')
    ax.axis('off')
    rcolors = plt.cm.BuPu(np.full(len(df.index), 0.1))
    ccolors = plt.cm.BuPu(np.full(len(df.columns), 0.1))

    the_table = ax.table(cellText=df.values,
            rowLabels=df.index,
            rowColours=rcolors,
            rowLoc='right',
            colColours=ccolors,  
            colLabels=df.columns,
            loc='center')
    plt.tight_layout
    plt.savefig(f'Output/{rundirs[0]}_NIThv_tbl_stats.png', dpi=300, bbox_inches='tight')
    plt.close()

    #pp = PdfPages("Output/Limits_isotherm_tbl_stats.pdf")
    #pp.savefig(fig, bbox_inches='tight')
    #pp.close()

if __name__ == "__main__":
    main()
