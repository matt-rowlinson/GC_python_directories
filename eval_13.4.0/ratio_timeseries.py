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
    ds = ds[f'SpeciesConc_{variable}'].isel(lev=0)
    return ds

def site_data(ds1, ds, lat=16.9, lon=-24.9, lev=0):
    x = rp.find_nearest(ds.lon, lon)-1
    y = rp.find_nearest(ds.lat, lat)+1
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
    plt.savefig( f"plots/{site['save_name']}_{v}.png" )
    plt.close()
    return

def find_model_output_for_site(ds, rundirs, versions, v, years='2017'):
    timeseries=[] ; add_lowes=[]
    for n, rundir in enumerate(rundirs):
        ds0 = get_data_as_xr(rundir, year=years[n], variable=v, version=versions[n])

        ds0 = site_data( ds0, ds, lon=site['longitude'], lat=site['latitude'] )
        ds0 = pd.DataFrame({v:ds0.values}, index=ds0['time'].values)#.resample('M').mean()
        
        timeseries.append( ds0 )
    return timeseries

##----------------------------------------MAIN SCRIPT----------------------------------------##
def main():
    global d, site, versions, rundirs, v
    ds = get_data_as_xr('nitrate_photol_control', variable='O3', year='201001')
    rundir, variables, version, site = get_arguments()
    sys.path.append( '/users/mjr583/cvao' )
    from lowess_smoother import loess

    site = sites[site]
    rundirs  = ['old_ceds_tropchem_merra_4x5','nitrate_photol_2017_control','TEST']
    versions = ['12.9.3','13.1.2_DUST','13.4.0']
    variables=['SALCAL','SALC']
    years=['2016','2018','2018']
    cs=['#1b9e77','#e7298a','#d95f02','#7570b3','#e7298a','#66a61e','#e6ab02'][::-1]
    vs=[]
    for v in variables:
        print( v )
        add_times = find_model_output_for_site(ds, rundirs, versions, v, years=years) 
        vs.append( add_times )
    
    v129 = vs[0][0].SALCAL.values / vs[1][0].SALC.values
    v131 = vs[0][1].SALCAL.values / vs[1][1].SALC.values
    v134 = vs[0][2].SALCAL.values / vs[1][2].SALC.values
    
    i2018=[]
    for dt in vs[0][0].index:
        dt = dt.replace( year=2018 )
        i2018.append( dt )
    vs[0][0].index = i2018
    
    f, ax = plt.subplots(1,1,figsize=(10,4))
    ax.plot( vs[0][0].index, v129, label='v12.9.3', c=cs[2] )
    ax.plot( vs[0][1].index, v131, label='v13.1.2', c=cs[3] )
    ax.plot( vs[0][1].index, v134, label='v13.4.0', c=cs[4] )
    plt.legend()
    plt.ylabel( "SALCAL / SALC" )
    plt.suptitle( site["site_name"]+": SALCAL / SALC" )
    plt.savefig( f'plots/{site["save_name"]}_SALC-SALCAL.png' )
    plt.close()

if __name__ == "__main__":
    main()
