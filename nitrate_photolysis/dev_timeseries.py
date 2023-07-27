#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#SBATCH --job-name=devtimeseries
#SBATCH --ntasks=1
#SBATCH --mem=18Gb
#SBATCH --partition=interactive
#SBATCH --time=03:30:00
#SBATCH --output=Logs/devtimeseries_%A.log
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import sys
import os
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
import pandas as pd
matplotlib.use('agg')
sys.path.append('/users/mjr583/python_lib')
from sites_dicts import GAW_dict as sites
from CVAO_dict import CVAO_dict as d
import RowPy as rp
import argparse
from scipy.optimize import curve_fit
from sklearn.metrics import mean_squared_error, r2_score
plt.style.use('seaborn-darkgrid')

def get_arguments():
    parser = argparse.ArgumentParser(description="Parse arguments to pass to GC processing scripts")
    parser.add_argument("-r", "--rundir", type=str, 
                        help='Name of desired GC rundir')
    parser.add_argument("-v", "--var", type=list,
                        default=["O3"],
                        help="Name of GC variable")
    parser.add_argument("-V", "--version", type=str,
                        default='13.1.2',
                        help="Version of GEOS-Chem")
    parser.add_argument("-s", "--site", type=str,
                        default="CVO",
                        help="GAW site of interest")
    args=parser.parse_args()
    return args.rundir, args.var, args.version, args.site

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
        if variable == 'NOx':
            ds = ( ds[f'SpeciesConc_NO'] + ds[f'SpeciesConc_NO2'] ).isel(lev=0) * 1e12
        elif variable == 'all_nitrate':
            try:
                ds = ( ds[f'SpeciesConc_NIT']   + ds[f'SpeciesConc_NITs']  + \
                       ds[f'SpeciesConc_NITD1'] + ds[f'SpeciesConc_NITD2'] + \
                       ds[f'SpeciesConc_NITD3'] + ds[f'SpeciesConc_NITD4']).isel(lev=0) * 1e12
            except:
                ds = ( ds[f'SpeciesConc_NIT'] + ds[f'SpeciesConc_NITs'] ).isel(lev=0) * 1e12
        elif variable == 'all_sulphate':
            try:
                ds = ( ds[f'SpeciesConc_SO4']   + ds[f'SpeciesConc_SO4s']  + \
                       ds[f'SpeciesConc_SO4D1'] + ds[f'SpeciesConc_SO4D2'] + \
                       ds[f'SpeciesConc_SO4D3'] + ds[f'SpeciesConc_SO4D4']).isel(lev=0) * 1e12
            except:
                ds = ( ds[f'SpeciesConc_SO4'] + ds[f'SpeciesConc_SO4s'] ).isel(lev=0) * 1e12
        elif variable=='NOy':
            ds = ( ds[f'SpeciesConc_NO'] + ds[f'SpeciesConc_NO2'] + ds[f'SpeciesConc_N2O5'] \
                 + ds[f'SpeciesConc_ClNO2'] + ds[f'SpeciesConc_HNO3'] + ds[f'SpeciesConc_NIT']  \
                 + ds[f'SpeciesConc_NITs'] + ds[f'SpeciesConc_HNO2'] + ds[f'SpeciesConc_NH3']   \
                 + ds[f'SpeciesConc_PAN'] + ds[f'SpeciesConc_HNO4'] \
                         ).isel(lev=0) * 1e12
        elif variable=='DST':
            ds = (  ds[f'SpeciesConc_DST1'] + ds['SpeciesConc_DST2'] + ds['SpeciesConc_DST3'] + \
                    ds['SpeciesConc_DST4'] ).isel(lev=0) * 1e12
        elif variable == 'all_nitrate':
            ds = ( ds[f'SpeciesConc_NIT']   + ds[f'SpeciesConc_NITs']  + \
                   ds[f'SpeciesConc_NITD1'] + ds[f'SpeciesConc_NITD2'] + \
                   ds[f'SpeciesConc_NITD3'] + ds[f'SpeciesConc_NITD4']).isel(lev=0) * 1e9

        elif variable == 'all-nitrate-plus-HNO3':
            ds = ( ds[f'SpeciesConc_NIT']   + ds[f'SpeciesConc_NITs']  + \
                   ds[f'SpeciesConc_NITD1'] + ds[f'SpeciesConc_NITD2'] + \
                   ds[f'SpeciesConc_NITD3'] + ds[f'SpeciesConc_NITD4'] + ds[f'SpeciesConc_HNO3']).isel(lev=0) * 1e9

        elif variable == 'all_sulphate':
            try:
                ds = ( ds[f'SpeciesConc_SO4']   + ds[f'SpeciesConc_SO4s']  + \
                       ds[f'SpeciesConc_SO4D1'] + ds[f'SpeciesConc_SO4D2'] + \
                       ds[f'SpeciesConc_SO4D3'] + ds[f'SpeciesConc_SO4D4']).isel(lev=0) * 1e9
            except:
                ds = ( ds[f'SpeciesConc_SO4'] + ds[f'SpeciesConc_SO4s'] ).isel(lev=0) * 1e9
        elif variable=='Bry':
            ds = ( ds[f'SpeciesConc_Br'] + (ds[f'SpeciesConc_Br2']*2) + ds[f'SpeciesConc_HOBr'] \
                 + ds[f'SpeciesConc_BrO'] + ds[f'SpeciesConc_HBr'] + ds[f'SpeciesConc_BrNO2']  \
                 + ds[f'SpeciesConc_BrNO3'] + ds[f'SpeciesConc_IBr'] + ds[f'SpeciesConc_BrCl']   \
                         ).isel(lev=0) * 1e12
        elif variable=='Cly':
            ds = ( ds[f'SpeciesConc_Cl'] + (ds[f'SpeciesConc_Cl2']*2) + ds[f'SpeciesConc_HOCl'] \
                 + ds[f'SpeciesConc_ClO'] + ds[f'SpeciesConc_HCl'] + ds[f'SpeciesConc_ClNO2']  \
                 + ds[f'SpeciesConc_ClNO3'] + ds[f'SpeciesConc_ICl'] + ds[f'SpeciesConc_BrCl']   \
                 + ds[f'SpeciesConc_ClOO'] + ds[f'SpeciesConc_OClO'] + (ds[f'SpeciesConc_Cl2O2']*2)   \
                         ).isel(lev=0) * 1e12
        elif variable=='Iy':
            ds = ( ds[f'SpeciesConc_I'] + (ds[f'SpeciesConc_I2']*2) + ds[f'SpeciesConc_HOI'] \
                 + ds[f'SpeciesConc_IO'] + ds[f'SpeciesConc_OIO'] + ds[f'SpeciesConc_HI']  \
                 + ds[f'SpeciesConc_INO'] \
                 + (ds[f'SpeciesConc_I2O2']*2) + (ds[f'SpeciesConc_I2O3']*2) + (ds[f'SpeciesConc_I2O4']*2)   \
                         ).isel(lev=0) * 1e12
        elif Chem:
            ds = ds[f'{variable}concAfterChem'].isel(lev=0) 
        else:
            try:
                ds = ds[f'SpeciesConc_{variable}'].isel(lev=0) * float(d[variable]['scale'])
            except:
                ds = ds[f'SpeciesConc_{variable}'].isel(lev=0) * 1e12
    return ds

def site_data(ds1, ds, lat=16.9, lon=-24.9, lev=0):
    if site['save_name']=='cvao':
        x = rp.find_nearest(ds.lon, lon)-1
        y = rp.find_nearest(ds.lat, lat)+1
    else:
        x = rp.find_nearest(ds.lon, lon)
        y = rp.find_nearest(ds.lat, lat)

    if type(lev)==int:
        data = ds1.isel( lon=x, lat=y )
    else:
        data = ds1.isel( lon=x, lat=y ).isel( lev=slice(lev)).mean(dim='lev', keep_attrs=True)
    return data

def load_observations( site, variable ):
    if site['site_name'] == 'CVAO':
        df = pd.read_csv( site['filepath'], index_col=0, dtype={"Airmass":str, "New_Airmass":str}, low_memory=False)
        if variable=='NOx':
            df = df['NO_pptV'] + df['NO2_pptV'] 
        else:
            df = df[d[variable]['merge_name']]
    elif site['site_name'] == 'Cape Grim':
        df = pd.read_csv( site['filepath']+site[f'{variable}_file'],index_col=0)
        df = df['Value']
        df = df / 1.9957
        df.index = pd.to_datetime( df.index, format="%d/%m/%Y %H:%M")
    else:
        df = pd.read_csv( site['filepath']+site[f'{variable}_file'],index_col=0)
        df = df['Value']
    if site['unit_conv'] == True:
        df = df / 1.96
    df.index = pd.to_datetime( df.index, format="%Y-%m-%d %H:%M:%S")
    return df

def plot(ax, times, obs=False,
                    labels='',
                    sname="",
                    cs=['#1b9e77','#e7298a','#d95f02','#7570b3','#e7298a','#66a61e','#e6ab02'],
                    style=['-','--',':','--',':','--',':']):
    for n in range(len(times)):
        ax.plot( times[n].index, times[n][site["save_name"]], label=labels[n],
                 c=cs[n], zorder=2, ls=style[n])
    if type(obs)!=bool:
        obs.plot(c='k',zorder=1, alpha=.75, label=site['site_name'], ax=ax)

    plt.ylabel(f'{v} pptv')
    plt.legend(loc=0, ncol=2)
    plt.tight_layout()
    plt.savefig( f"plots/{site['save_name']}_{v}_{sname}.png" )
    plt.close()
    return

def find_model_output_for_site(ds, rundirs, versions, v, years='2017'):
    timeseries=[] ; add_lowes=[]
    for n, rundir in enumerate(rundirs):
        #print( rundir )
        if len( versions ) == 1:
            version=versions[0]
        else:
            version=versions[n]
        ds0 = get_data_as_xr(rundir, year=years, variable=v, version=versions[0])

        ds0 = site_data( ds0, ds, lon=site['longitude'], lat=site['latitude'] )
        ds0 = pd.DataFrame({site['save_name']:ds0.values}, index=ds0['time'].values)#.resample('M').mean()
        if site['save_name']=='cvao':
            ds0=ds0['2006-10-01':]
        timeseries.append( ds0 )
    return timeseries

##----------------------------------------MAIN SCRIPT----------------------------------------##
def main():
    global d, site, versions, rundirs, v
    ds = get_data_as_xr('nitrate_photol_control', year='201001')
    rundir, variables, version, site = get_arguments()
    
    all_sites=['CVO','MAC','NEU','RAG','FIREX','FIREX2']
    for s in all_sites:
        print( s )
        snames=['isotherm','isotherm_Limits']
        for sname in snames:
            print( sname)
            site = sites[s]
            if "Limits" in sname:
                rundirs  = ['dev_new_base','dev_onlyNO-HO2','dev_onlyNIThv','dev_both','isotherm_ON','dev_lim_HO2nNO']
                labels=['Base','._NO+HO2','._NIThv','._NIThv_NO+HO2','._NIThvLimit','._NIThvLimit_NO+HO2']
            else:
                rundirs  = ['dev_new_base','dev_onlyNO-HO2','dev_onlyNIThv','dev_both']
                labels=['Base','._NO+HO2','._NIThv','._NIThv_NO+HO2']
            
            variables=['O3','NOx','HNO2','HNO3','all_nitrate', 'all_sulphate','SO4','SO2','CO','OH']
            years='2018'
            cs=['#1b9e77','#d95f02','#7570b3','#e7298a','#66a61e','#e6ab02']
            for v in variables:
                add_times = find_model_output_for_site(ds, rundirs, ['13.1.2'], v, years=years)
                f, ax = plt.subplots(1,1,figsize=(10,4))
                try:
                    df = load_observations(site, v)[:]#.resample("M").mean()
                    df=df['2018'].resample('D').mean()
                    plot( ax, add_times, labels=labels, obs=df, cs=cs, sname=sname)
                except:
                    plot(ax, add_times, labels=labels, cs=cs, sname=sname)

if __name__ == "__main__":
    main()
