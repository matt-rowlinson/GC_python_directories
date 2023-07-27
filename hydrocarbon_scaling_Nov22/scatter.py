#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#SBATCH --job-name=jscale_time
#SBATCH --ntasks=1
#SBATCH --mem=18Gb
#SBATCH --partition=interactive
#SBATCH --time=00:30:00
#SBATCH --output=Logs/timeseries_%A.log
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import sys
import os
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.colors import LogNorm
import numpy as np
import pandas as pd
matplotlib.use('agg')
sys.path.append('/users/mjr583/python_lib')
from sites_dicts import GAW_dict as sites
import RowPy as rp
plt.style.use('seaborn-darkgrid')

import yaml
from sklearn.metrics import mean_squared_error
from math import sqrt

def find_file_list(path, substrs):
    file_list =[]
    for root, directory, files in os.walk(path+'/OutputDir/'):
        for f in files:
            for s in substrs:
                if s in f:
                    file_list.append(os.path.join(root, f))
    file_list.sort()
    return file_list

def open_family_yaml():
    with open("family_variables.yml", "r") as stream:
        try:
            d = yaml.safe_load(stream)
        except yaml.YAMLError as exc:
            print(exc)
    return d

def open_species_dict():
    sys.path.append('/users/mjr583/python_lib')
    from CVAO_dict import CVAO_dict as d
    return d

def get_data_as_xr(rundir, year='', collection='SpeciesConc', variable=False, Chem=False, version='13.1.2'):
    family_dict  = open_family_yaml()
    species_dict = open_species_dict()

    path=f'/mnt/lustre/users/mjr583/GC/{version}/rundirs/{rundir}'
    try:
        ds = xr.open_mfdataset( find_file_list(path, [f'{collection}.{year}']), combine='by_coords' )
    except:
        print( 'Incomplete run' )
        ds = xr.open_mfdataset( find_file_list(path, [f'{collection}.{year}'])[:-1], combine='by_coords' )

    if f"{collection}_{variable}" in ds.variables.keys():
        ### Read single variable from SpeciesConc           ###
        ds = ds[f'SpeciesConc_{variable}'] * float(species_dict[variable]['scale'])
    elif f"{variable}concAfterChem" in ds.variables.keys():
        ### Read single variable from ConcAfterChem         ###
        ds = ds[f'SpeciesConc_{variable}'] * float(species_dict[variable]['scale'])
    elif variable in family_dict.keys():
        ### Read family of variables to read i.e. NOx, Bry  ###
        ds_=0
        for v in family_dict[variable]:
            ds_ += ds[f'SpeciesConc_{v}'] 
        ds = ds_ * float(species_dict[variable]['scale'])
    return ds

def site_data(ds, lat=16.9, lon=-24.9, lev=0):
    x = rp.find_nearest(ds.lon, lon)-1
    y = rp.find_nearest(ds.lat, lat)-1
    if type(lev)==int:
        data = ds.isel( lon=x, lat=y, lev=lev )
    else:
        data = ds.isel( lon=x, lat=y ).isel( lev=slice(lev)).mean(dim='lev', keep_attrs=True)
    return data

def load_observations( site, variable, year ):
    species_dict = open_species_dict()

    df = pd.read_csv( site['filepath'], index_col=0, dtype={"Airmass":str, "New_Airmass":str}, low_memory=False)
    if variable=='NOx':
        df = df['NO_pptV'] + df['NO2_pptV']
    elif variable=='ALK4':
        iB = df['iso_butane'].fillna(0)
        nB = df['n_butane'].fillna(0)
        iP = df['iso_pentane'].fillna(0)
        nP = df['n_pentane'].fillna(0)
        df = iB + nB + iP + nP
    else:
        df = df[species_dict[variable]['merge_name']]

    if site['unit_conv'] == True:
        df = df / 1.96
    df.index = pd.to_datetime( df.index, format="%Y-%m-%d %H:%M:%S")
    return df[year]

def plot(ax, times, obs=False,
                    labels='',
                    sname="Jscale_25.png",
                    cs=['#1b9e77','#e7298a','#d95f02','#7570b3','#e7298a','#66a61e','#e6ab02'],
                    style=['-','--','--',':',':']):
    for n in range(len(times)):
        ax.plot( times[n].index, times[n][site["save_name"]], label=labels[n],
                 c=cs[n], zorder=2, ls=style[n])
    if type(obs)!=bool:
        obs.plot(c='k',zorder=1, alpha=.5, label=site['site_name'], ax=ax)

    plt.ylabel(f'{v} pptv')
    plt.legend(loc=0, ncol=2)
    plt.tight_layout()
    plt.savefig( f"plots/TEST.Ander22b.{site['save_name']}_{v}_Limits.png" )
    plt.close()
    return

def read_model_output( rundirs, versions, v, years='2017'):
    timeseries=[] 
    for n, rundir in enumerate(rundirs):
        print( rundir )
        if len( versions ) == 1:
            version=versions[0]
        else:
            version=versions[n]
        ds = get_data_as_xr(rundir, year=years, variable=v, version=version)
        ds = site_data( ds, lon=site['longitude'], lat=site['latitude'], lev=0 )
        ds = pd.DataFrame({rundir:ds.values}, index=ds['time'].values)#.resample('M').mean()
        timeseries.append( ds )
    return timeseries

def calc_stats(obs, model):
    obs = obs.dropna()
    model - model.reindex( obs.index )
    absError = obs.values - model.values
    SE_error = np.square( absError )
    MSE = np.mean( SE_error )
    RMSE=np.round( np.sqrt(MSE), 2)
    R2 = np.round( 1. - (np.var(absError) / np.var(obs)), 2)
    return RMSE, R2

##----------------------------------------MAIN SCRIPT----------------------------------------##
def main():
    global site
    site     = sites['CVO']
    rundirs  = ['base_run','geo','scale_only','scale_all']
    versions = ['14.0.1']
    year     = '2017'
    variable = 'PRPE'

    ds = read_model_output( rundirs, versions, variable, years=year)
    df = load_observations(site, variable, year).resample("D").mean()
    
    f, ax = plt.subplots(1,1,figsize=(8,8))

    line = matplotlib.lines.Line2D([0, 1], [0, 1], color='darkgrey', alpha=1.)
    transform = ax.transAxes
    line.set_transform(transform)
    ax.add_line(line)

    labels=['Base','Geo only','Scale only','Scale + Geo']
    cs = ['#377eb8','#4daf4a','#984ea3','#ff7f00']
    for n,nn in enumerate(ds):
        ds[n].index = df.index
        rmse, r2 = calc_stats( df, ds[n] )
        dff = df.dropna()
        d  = ds[n].reindex(dff.index)
        a, _, _, _ = np.linalg.lstsq(dff.values[:,np.newaxis], d, rcond=None)
        rmse2, r2 = calc_stats( dff, a[0] * dff )
        ax.scatter( df, ds[n], c=cs[n], label=f'{labels[n]} (RMSE={rmse} ppt, R$^2$={r2})')
        ax.plot(dff, a[0] * dff, '-', c=cs[n])

    MAX = np.nanmax( [np.nanmax(df.values), np.max(ds)] )
    ax.set_ylim( 0, MAX*1.01 ) ; ax.set_xlim( 0, MAX*1.01 )
    ax.set_ylabel( f'Modelled {variable} (pptv)' )
    ax.set_xlabel( f'Observed {variable} (pptv)' )

    plt.legend()
    plt.tight_layout()
    plt.savefig(f'plots/scatter.{site["save_name"]}.{variable}.png')
    plt.close()
    
if __name__ == "__main__":
    main()
