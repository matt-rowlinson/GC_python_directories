#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#SBATCH --job-name=table_stats
#SBATCH --ntasks=1
#SBATCH --mem=52Gb
#SBATCH --partition=nodes
#SBATCH --time=00:30:00
#SBATCH --output=Logs/table_%A.log
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
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
    cmap = ['#377eb8', '#ff7f00', '#4daf4a','#f781bf', '#a65628', '#984ea3','#999999', '#e41a1c', '#dede00']

    df = pd.read_csv('csvs/scalers.csv', index_col=0)  
    af = df[df>=1.].fillna(0.)
    bf = -(1/df[df<=1.]).fillna(0.)
    df = af + bf 
    df = df.replace( [np.inf, -np.inf], 0.)#, inplace=True )
    df['Mean SF'] = df.mean(axis=1)
    
    df.plot(kind='bar', figsize=(6,5.5), rot=45, color=cmap, hatch=(df['Mean SF'] >= 0.).map({True:None, False:'x'}))
    plt.ylabel( 'Emission scaling' )
    plt.title( 'VOC Emission Scale Factors')
    plt.tight_layout()
    plt.savefig('Output/scale_factor_bars.png')
    plt.close()
    sys.exit()

    core_spec=['O3','C3H8','C2H6','ALK4','BENZ','EOH','CH2O','TOLU','MEK','PRPE','XYLE']
    unit_scal=[1e9 , 1e12 , 1e12 , 1e12 , 1e12 , 1e12, 1e12 , 1e12 , 1e12, 1e12 , 1e12 ]
    rundirs  = ['geo','scale_all']
    labels=['Base', 'Scaled CEDS']
    version = '14.0.1'
    b=[] 
    for n, rundir in enumerate(rundirs[:]):
        print( rundir )
        ds  = get_data_as_xr(rundir, version=version, year='2016')
        Met = get_data_as_xr('base_run_2010', version='13.1.2', collection='StateMet', year='2010')
        trop = Met['Met_TropP'].mean( dim='time' )
        AirMass = Met['Met_AD'].mean( dim='time' )
        pmid_press = Met['Met_PMID'].mean( dim='time' )
        MASK =  ( pmid_press > trop )
        a=[] ; col=[]
         
        for nn, spec in enumerate(core_spec):
            if unit_scal[nn] == 1e9:
                unit='ppbv'
            else:
                unit='pptv'
            var = ds[f'SpeciesConc_{spec}']
            var = var.mean(dim='time')

            var = var.where(MASK)
            RMM_air = (.78*(2.*14.)+.22*(2.*16.))
            conversion_factor = (AirMass*1E3 / RMM_air)
            var = var * conversion_factor
            var = var * float(spec_db[core_spec[nn]]['MW_g'] * 1e-3 ) / 1E9

            a.append( (var.sum().values) ) 
            col.append( spec+' burden (Tg yr$^{-1}$)')

        b.append( a )

    c = np.array( b )
    df = pd.DataFrame( c, index=labels, 
            columns= col )
    df = df.transpose()
    df['% difference'] = ((df['Scaled CEDS'].values - df['Base'].values )/ df['Base'].values ) * 100
    df['Abs difference'] = df['Scaled CEDS'].values - df['Base'].values  
    df = df.round(2)
    #df.to_csv(f"Output/VOC_burdens_4x5_stats.csv")

    fig, ax =plt.subplots(figsize=(5.5,2.5))
    ax.axis('tight')
    ax.axis('off')
    rcolors = plt.cm.BuPu(np.full(len(df.index), 0.1))
    ccolors = plt.cm.BuPu(np.full(len(df.columns), 0.1))

    the_table = ax.table(cellText=df.values,
            rowLabels=df.index,
            rowColours=rcolors,
            rowLoc='center',
            colLoc='center',
            cellLoc='center',
            colColours=ccolors,  
            colLabels=df.columns,
            loc='center')
    plt.tight_layout
    plt.savefig(f'Output/VOC_burdens_4x5_tbl.png', dpi=300, bbox_inches='tight')
    plt.close()

    df = df.transpose()
    df = df[[col for col in df.columns if 'burden' in col or 'OH' in col]]
    df.columns = core_spec
    df = df.transpose() 
    df['positive'] = df['% difference'] > 0
    df['% difference'].plot(kind='bar', rot=45, color=df.positive.map({True:'r', False:'b'}))
    plt.ylabel( '% difference after CEDS scaling' )
    plt.title( 'VOC burden difference')
    plt.tight_layout()
    plt.savefig('Output/voc_burdens_bars.png')
    plt.close()



if __name__ == "__main__":
    main()
