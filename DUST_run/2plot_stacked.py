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
        ds = ds[f'SpeciesConc_{variable}'].isel(lev=0) * 1e12
    return ds

def site_data(ds1, ds, lat=16.9, lon=-24.9, lev=0):
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
        df = df[variable]
    if site['unit_conv'] == True:
        df = df / 1.96
    df.index = pd.to_datetime( df.index, format="%Y-%m-%d %H:%M:%S")
    return df

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

def find_model_output_for_site(ds, rundirs, versions, v, years='2018'):
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
    rundirs  = ['nitrate_photol_control']
    versions = ['13.1.2']
    variables=['NIT','NITs']
    years='2017'
    cs=['#1b9e77','#d95f02','#7570b3','#e7298a','#66a61e','#e6ab02']
    colordict = {'NIT' : cs[0], 'NITs' : cs[1], 'NITD1' : cs[2], 'NITD2' : cs[3], 'NITD3' : cs[4], 'NITD4' : cs[5],}
    vs=[]
    for v in variables:
        add_times = find_model_output_for_site(ds, rundirs, versions, v, years=years)[0] 
        vs.append( pd.DataFrame({ v : np.squeeze(add_times.values) }, index=add_times.index ) )
    vs1 = pd.concat( vs, axis=1 ).resample( 'M' ).mean()
    
    variables=['NIT','NITs','NITD1','NITD2','NITD3','NITD4']
    vs=[]
    for v in variables:
        add_times = find_model_output_for_site(ds, ['nitrate_photol_control_DUST'], ['13.1.2'], v, years=years)[0] 
        vs.append( pd.DataFrame({ v : np.squeeze(add_times.values) }, index=add_times.index ) )
    vs2 = pd.concat( vs, axis=1 ).resample( 'M' ).mean()
     
    f, (ax1,ax2) = plt.subplots(1,2,figsize=(16,6), sharey=True)
    vs1.plot.area(ax=ax1,color=[colordict.get(x, '#333333') for x in vs1.columns])
    ax1.legend(ncol=6)
    vs2.plot.area(ax=ax2,color=[colordict.get(x, '#333333') for x in vs2.columns])
    ax2.legend(ncol=6)

    plt.ylabel( 'NIT aerosol (pptv)')
    plt.tight_layout()
    plt.savefig(f'plots/{site["site_name"].replace(" ","")}_stacked_NIT.png')
    plt.close()
    
    f, ax = plt.subplots(1,1,figsize=(12,6))

    plt.plot( vs1.index, vs1.NIT, 'k--', zorder=101)
    plt.plot( vs1.index, vs1.NIT.values + vs1.NITs.values, 'k--', label='No dust', zorder=101)

    plt.stackplot( vs2.index, vs2.NIT, vs2.NITs, vs2.NITD1, vs2.NITD2, vs2.NITD3, vs2.NITD4, labels=vs2.columns, colors=cs )  
    plt.stackplot( vs1.index, vs1.NIT, vs1.NITs, alpha=.5, colors=cs)
    plt.legend( ncol=7 )

    plt.ylabel( 'NIT aerosol (pptv)')
    plt.tight_layout()
    plt.savefig(f'plots/overlayed_{site["site_name"].replace(" ","")}_stacked_NIT.png')
    plt.close() 

if __name__ == "__main__":
    main()
