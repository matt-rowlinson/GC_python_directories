#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#SBATCH --job-name=jscale_time
#SBATCH --ntasks=1
#SBATCH --mem=18Gb
#SBATCH --partition=interactive
#SBATCH --time=00:30:00
#SBATCH --output=Logs/jscale_time_plots_%A.log
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
    for root, directory, files in os.walk(path):
        for f in files:
            for s in substrs:
                if s in f:
                    file_list.append(os.path.join(root, f))
    file_list.sort()
    return file_list

def get_data_as_xr(rundir, year='', variable=False, version='13.1.2'):
    path=f'/mnt/lustre/users/mjr583/GC/{version}/rundirs/{rundir}'
    ds = xr.open_mfdataset( find_file_list(path, [f'SpeciesConc.{year}']), combine='by_coords' )
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
        else:
            ds = ds[f'SpeciesConc_{variable}'].isel(lev=0) * 1e12# float(d[variable]['scale'])
    return ds

def site_data(ds1, ds, lat=16.9, lon=-24.9, lev=0):
    x = rp.find_nearest(ds.lon, lon)-1
    y = rp.find_nearest(ds.lat, lat)+1
    if type(lev)==int:
        data = ds1.isel( lon=x, lat=y )
    else:
        data = ds1.isel( lon=x, lat=y ).isel( lev=slice(lev)).mean(dim='lev', keep_attrs=True)
    return data

def plot(ax, times, obs=False,
                    sname="Jscale_25.png",
                    cs=['#1b9e77','#e7298a','#d95f02','#7570b3','#e7298a','#66a61e','#e6ab02'],
                    style=['-','-','--',':',':']):
    for n in range(len(times)):
        ax.plot( times[n].index, times[n][site["save_name"]], label=versions[n].replace('nitrate_photol_','').title(), 
                 c=cs[n], zorder=2, ls=style[n])
    if type(obs)!=bool:
        obs.plot(c='k',zorder=1, alpha=.5, label=site['site_name'], ax=ax)

    plt.ylabel(f'{v} pptv')
    plt.legend(loc=0, ncol=2)
    plt.tight_layout()
    plt.savefig( f"plots/{site['save_name']}_{v}_{sname}" )
    plt.close()
    return

def find_model_output_for_site(ds, rundirs, versions, v, years='2017'):
    timeseries=[] ; add_lowes=[]
    for n, rundir in enumerate(rundirs):
        ds0 = get_data_as_xr(rundir, year=years, variable=v, version=versions[n])
      
        ds0 = site_data( ds0, ds, lon=site['longitude'], lat=site['latitude'] )
        ds0 = pd.DataFrame({site['save_name']:ds0.values}, index=ds0['time'].values)#.resample('M').mean()
        
        if site['save_name']=='cvao':
            ds0=ds0['2006-10-01':]
        ds0=ds0['2017':'2018']
        timeseries.append( ds0 )
    return timeseries

##----------------------------------------MAIN SCRIPT----------------------------------------##
def main():
    global d, site, versions, rundirs, v
    ds = get_data_as_xr('nitrate_photol_control', year='201001')
    rundir, variables, version, site = get_arguments()
    sys.path.append( '/users/mjr583/cvao' )
    from lowess_smoother import loess

    site = sites[site]
    rundirs  = ['nitrate_photol_2017_control', 'nitrate_photol_control']
    versions = ['13.1.2','13.1.2_DUST']
    sname='DUST'
    years='2018'
    cs=['#1b9e77','#e7298a','#d95f02','#7570b3','#e7298a','#66a61e','#e6ab02']
    
    cNH3 = find_model_output_for_site(ds, ['nitrate_photol_2017_control'], ['13.1.2'], 'SO4', years=years)[0] 
    cNH4 = find_model_output_for_site(ds, ['nitrate_photol_2017_control'], ['13.1.2'], 'SO4s', years=years)[0]
    cTOT = cNH3 + cNH4

    dNH3 = find_model_output_for_site(ds, ['nitrate_photol_control'], ['13.1.2_DUST'], 'SO4', years=years)[0]
    dNH4 = find_model_output_for_site(ds, ['nitrate_photol_control'], ['13.1.2_DUST'], 'SO4s', years=years)[0]
    dSO4D1 = find_model_output_for_site(ds, ['nitrate_photol_control'], ['13.1.2_DUST'], 'SO4D1', years=years)[0]
    dSO4D2 = find_model_output_for_site(ds, ['nitrate_photol_control'], ['13.1.2_DUST'], 'SO4D2', years=years)[0]
    dSO4D3 = find_model_output_for_site(ds, ['nitrate_photol_control'], ['13.1.2_DUST'], 'SO4D3', years=years)[0]
    dSO4D4 = find_model_output_for_site(ds, ['nitrate_photol_control'], ['13.1.2_DUST'], 'SO4D4', years=years)[0]
    dTOT = dNH3 + dNH4 + dSO4D1 + dSO4D2 + dSO4D3 + dSO4D4


    cs=['#d95f02','#7570b3','#66a61e','#e6ab02']
    f,ax=plt.subplots(figsize=(12,6))
    #ax.plot( cNH3.index, cNH3.values, label='SO4 v13.1.2', c=cs[0], zorder=2, ls='-')
    #ax.plot( cNH4.index, cNH4.values, label='SO4s v13.1.2', c=cs[1], zorder=2, ls='-')
    ax.plot( cNH4.index, cTOT.values, label='Total SO4 (no dust) v13.1.2', c=cs[2], zorder=2, ls='-')

    #ax.plot( dNH3.index, dNH3.values, label='SO4 v13.1.2_DUST', c=cs[0], zorder=2, ls='--')
    #ax.plot( dNH4.index, dNH4.values, label='SO4s v13.1.2_DUST', c=cs[1], zorder=2, ls='--')
    ax.plot( cNH4.index, dTOT.values, label='Total SO4 (inc. dust) v13.1.2', c=cs[2], zorder=2, ls='--')

    plt.ylabel(f' pptv')
    plt.legend(loc=0, ncol=2)
    plt.tight_layout()
    plt.savefig( f"plots/totalSO4_DUST.png" )
    plt.close()


if __name__ == "__main__":
    main()
