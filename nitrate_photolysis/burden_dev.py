#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#SBATCH --job-name=jscale_time
#SBATCH --ntasks=1
#SBATCH --mem=18Gb
#SBATCH --partition=test
#SBATCH --time=00:30:00
#SBATCH --output=Logs/timeseries_%A.log
#import warnings
#warnings.simplefilter(action='ignore', category=FutureWarning)
import sys
import os
import xarray as xr
import matplotlib.pyplot as plt
#import cartopy.crs as ccrs
import matplotlib
#from matplotlib.colors import LogNorm
#from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
import pandas as pd
matplotlib.use('agg')
#sys.path.append('/users/mjr583/python_lib')
#from sites_dicts import GAW_dict as sites
#from CVAO_dict import CVAO_dict as d
#import RowPy as rp
#import argparse
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

def get_data_as_xr(rundir, year='', collection='SpeciesConc', archive=False, archive_dir='',variable=False, Chem=False, version='14.1.0'):
    path=f'/mnt/lustre/users/mjr583/GC/{version}/rundirs/{rundir}'

    if variable=='OH' or variable=='HO2':
        collection='ConcAfterChem'
        Chem=True
    try:
        ds = xr.open_mfdataset( find_file_list(path, [f'{collection}.{year}']), combine='by_coords' )
    except:
        ds = xr.open_mfdataset( find_file_list(path, [f'{collection}.{year}'])[:-1], combine='by_coords' )

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

    core_spec=['O3','CO','NO','NO2', 'HNO2']
    Ycore_spec=['O3','CO','NO','NO2', 'HNO2']
    scaler = [1e9,1e9,1e12,1e12,1e12]
    rundirs  = ['new_base_4x5',
                'DD_new_base_4x5',
                'langmuir_4x5',
                'DD_langmuir_4x5',
                'langmuir_RH_90_4x5',
                'DD_langmuir_RH_90_4x5'                ]
    labels=['Base',
            'Base with HNO2 dep',
            'Langmuir',
            'Langmuir with HNO2 dep',
            'Langmuir max90',
            'Langmuir max90, with HNO2 dep']
    b=[]
    for n, rundir in enumerate(rundirs[:]):
        print( rundir )
        ds  = get_data_as_xr(rundir, year='2015')
        #print( ds )
        #sys.exit()
        CaC = get_data_as_xr(rundir, collection='ConcAfterChem', year='2015')
        PrLo = get_data_as_xr(rundir, collection='ProdLoss', year='2015')
        Met = get_data_as_xr(rundir, collection='StateMet', year='2015')
        trop = Met['Met_TropP'].mean( dim='time' )
        AirMass = Met['Met_AD'].mean( dim='time' )
        AirDen = Met['Met_AIRDEN'].mean( dim='time' )
        pmid_press = Met['Met_PMID'].mean( dim='time' )
        MASK =  ( pmid_press > trop )#.mean( dim='time' )
        
        a=[]
        for nn, spec in enumerate(core_spec):
            var = ds[f'SpeciesConcVV_{spec}']
            var = var.mean(dim='time')

            var = var.where(MASK)
            RMM_air = (.78*(2.*14.)+.22*(2.*16.))
            conversion_factor = (AirMass*1E3 / RMM_air)
            var = var * conversion_factor
            if spec=='HNO2':
                var = var * float(spec_db[Ycore_spec[nn]]['MW_g'] * 1e-3 ) / 1e6
            else:
                var = var * float(spec_db[Ycore_spec[nn]]['MW_g'] * 1e-3 ) / 1e9

            a.append( var.sum().values ) 
            a.append( ds[f'SpeciesConcVV_{spec}'].mean(dim='time').isel(lev=0).mean().values * scaler[nn] )
        
        AirDen = Met['Met_AIRDEN'].mean( dim='time' )
        AirVol = Met['Met_AIRVOL'].mean( dim='time' )
        OHmass = CaC['OHconcAfterChem'].mean(dim='time') * AirMass
        OHmass_sum = np.nansum( OHmass ) / np.nansum( AirMass )
        a.append( OHmass_sum * 1e-6 )
        print( OHmass_sum * 1e-6 ) 

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
    #print( c )
    #sys.exit()
    df = pd.DataFrame( c, index=labels, 
            columns= ['O3 burden (Tg)','O3 surface (ppbv)', 'CO burden (Tg)','CO surface (ppbv)',
                      'NO burden (Tg)','NO surface (pptv)','NO2 burden (Tg)','NO2 surface (pptv)',
                      'HNO2 burden (Gg)','HNO2 surface (pptv)',
                      'mw mean OH (1e6 molec cm-3)'] )
    df = df.transpose().round(2)
    df.to_csv("Output/Limits_stats.csv")

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
    plt.savefig('Output/wDep_tbl_stats.png', dpi=300, bbox_inches='tight')
    plt.close()

    #pp = PdfPages("Output/Limits_isotherm_tbl_stats.pdf")
    #pp.savefig(fig, bbox_inches='tight')
    #pp.close()

if __name__ == "__main__":
    main()
