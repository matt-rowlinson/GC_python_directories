#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#SBATCH --job-name=jscale_time
#SBATCH --ntasks=1
#SBATCH --mem=1000Gb
#SBATCH --partition=interactive
#SBATCH --time=01:30:00
#SBATCH --output=Logs/jscale_time_plots_%A.log
import sys
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import pandas as pd
matplotlib.use('agg')
sys.path.append('/users/mjr583/python_lib')
from sites_dicts import GAW_dict as sites
from CVAO_dict import CVAO_dict as d
import RowPy as rp
plt.style.use('seaborn-darkgrid')
import os
import xarray as xr

def find_file_list(path, substrs):
    file_list =[]
    for root, directory, files in os.walk(path+'/OutputDir/'):
        for f in files:
            for s in substrs:
                if s in f:
                    file_list.append(os.path.join(root, f))
    file_list.sort()
    return file_list

def get_data_as_xr(rundir, year='', collection='SpeciesConc', variable=False, Chem=False, version='13.1.2'):
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
            ds = xr.open_mfdataset( find_file_list(path, [f'{collection}.{year}'])[:-1], combine='by_coords' )

    if variable:
        try:
            ds = ds[f'SpeciesConc_{variable}'].isel(lev=0) * float(d[variable]['scale'])
        except:
            ds = ds[f'SpeciesConc_{variable}'].isel(lev=0) * 1e12
    return ds

def site_data(ds1, ds, lat=16.9, lon=-24.9, lev=0):
    #if site['save_name']=='cvao':
    #    x = rp.find_nearest(ds.lon, lon)-1
    #    y = rp.find_nearest(ds.lat, lat)+1
    #else:
    x = rp.find_nearest(ds.lon, lon)-1
    y = rp.find_nearest(ds.lat, lat)+1

    if type(lev)==int:
        data = ds1.isel( lon=x, lat=y )
    else:
        data = ds1.isel( lon=x, lat=y ).isel( lev=slice(lev)).mean(dim='lev', keep_attrs=True)
    return data

def find_model_output_for_site(ds, rundirs, versions, v, years='2017'):
    timeseries=[] ; add_lowes=[]
    for n, rundir in enumerate(rundirs):
        #print( rundir )
        if len( versions ) == 1:
            version=versions[0]
        else:
            version=versions[n]
        ds0 = get_data_as_xr(rundir, year=years, variable=v, version=versions[0])

        ds0 = site_data( ds0, ds )
        ds0 = pd.DataFrame({'CVAO':ds0.values}, index=ds0['time'].values)#.resample('M').mean()
        #if site['save_name']=='cvao':
        #    ds0=ds0['2006-10-01':]
        timeseries.append( ds0 )
    return timeseries

def main():
    #global d, site, rundirs, v
    sys.path.append('/users/mjr583/cvao/')
    from lowess_function import loess
    var='O3'
    alpha=.35 ; poly=1
    
    df0 = pd.read_csv('/mnt/lustre/users/mjr583/NCAS_CVAO/CVAO_datasets/Ebas_220512_1444/Ebas_ozone_20220512.csv',
                    index_col=0) / 1.9957
    df0.index = pd.to_datetime( df0.index, format="%Y-%m-%d %H:%M:%S").strftime('%Y-%m-%d %H:%M:%S')
    df0.index = pd.to_datetime( df0.index, format="%Y-%m-%d %H:%M:%S")
    df0['2009-07-01':'2009-09-01']=np.nan
    
    df2 = pd.DataFrame( {"XValue" : np.arange(len(df0.Value.values)), "YValue" : df0.Value.values   }, index=df0.index  )
    df2 = df2.resample('W').mean()
    df2=df2.dropna()
    eval_DF = loess("XValue", 'YValue', data = df2, alpha=alpha, poly_degree=poly )

    ds = get_data_as_xr('nitrate_photol_control', year='201001')
    ds = find_model_output_for_site(ds, ['nitrate_photol_control'], ['13.1.2'], var, years='20')[0]['2006-10':]

    ds1 = pd.DataFrame( {"XValue" : np.arange(len(ds.values)), "YValue" : ds.CVAO.values   }, index=ds.index  )
    ds1 = ds1.resample('W').mean()
    eval_DF2 = loess("XValue", 'YValue', data = ds1, alpha=alpha, poly_degree=poly )

    fig = plt.figure(figsize=(10,5))
    ax1 = fig.add_subplot(111)

    ds25 = df0.resample('M').quantile(.25)
    ds75 = df0.resample('M').quantile(.75)
    ax1.fill_between(ds25.index, ds25.Value.values, ds75.Value.values, color='k', alpha=.2)

    #ax1.scatter(df2.index, df2["YValue"], color="k", marker="o", alpha=0.3,s=5, label="_nolegend_")
    #ax1.plot(df0.Value.resample('M').mean().index, df0.Value.resample('M').mean(), color="k", alpha=0.7, label="_nolegend_")
    #ax1.plot(df2.index, eval_DF['g'][1:], color='k', linewidth= 3, label="CVAO")

    #ax1.scatter( ds.index, ds.CVAO, color="green", marker="o", alpha=0.3, s=5, label="_nolegend_")
    #ax1.plot(ds.resample('M').mean().index, ds.resample('M').mean(), color="green", alpha=0.7, label="_nolegend_")
    ds25 = ds.resample('M').quantile(.25)
    ds75 = ds.resample('M').quantile(.75)
    ax1.fill_between(ds25.index, ds25.CVAO.values, ds75.CVAO.values, color='green', alpha=.2)
    
    
    ax1.plot(df0.Value.resample('M').mean().index, df0.Value.resample('M').mean(), color="k", alpha=0.7, label="_nolegend_")
    ax1.plot(ds.resample('M').mean().index, ds.resample('M').mean(), color="green", alpha=0.7, label="_nolegend_")

    ax1.plot(df2.index, eval_DF['g'][1:], color='k', linewidth= 2.5, label="CVAO")
    ax1.plot(ds1.index, eval_DF2['g'][1:], color='green', linewidth= 2.5, label="GEOS-Chem")


    plt.xlabel(None), plt.ylabel(f"{var} / ppb")
    plt.legend()
    plt.tight_layout()
    plt.savefig(f'plots/cvao_{var}_lowess.png')
    plt.close()    

        
if __name__ == "__main__":
    main()
