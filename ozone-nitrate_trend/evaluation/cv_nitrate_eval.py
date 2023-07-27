#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#SBATCH --job-name=jscale_time
#SBATCH --ntasks=1
#SBATCH --mem=30Gb
#SBATCH --partition=interactive
#SBATCH --time=01:30:00
#SBATCH --output=Logs/jscale_time_plots_%A.log
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

def get_data_as_xr(rundir, year='', variable=False):
    path=f'/mnt/lustre/users/mjr583/GC/13.1.2/rundirs/{rundir}'
    ds = xr.open_mfdataset( find_file_list(path, [f'SpeciesConc.{year}']), combine='by_coords' )
    if variable:
        if variable == 'NOx':
            ds = ( ds[f'SpeciesConc_NO'] + ds[f'SpeciesConc_NO2'] ).isel(lev=0) * float(d[variable]['scale'])
        else:
            ds = ds[f'SpeciesConc_{variable}'].isel(lev=0) * float(d[variable]['scale'])
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
        df = pd.read_csv( site['filepath'], index_col=0, dtype={"Airmass":str, "New_Airmass":str})
        #df.index = pd.to_datetime( df.index, format='%Y-%m-%d %H:%M:%S' )
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

def plot(ax, df, times,
                    sname="Jscale_25.png",
                    cs=['#1b9e77','#e7298a','#d95f02','#7570b3','#e7298a','#66a61e','#e6ab02'],
                    style=['--','-',':']):
    for n in range(len(times)):
        ax.plot( times[n].index, times[n][site["save_name"]], label=rundirs[n].replace('nitrate_photol_','').title(), 
                 c=cs[n], zorder=2)
        #ax.plot( times[n].index, lowes[n]['g'][1:], c=cs[n])#, ls=style[n])
    
    df.plot(c='k',zorder=1, alpha=.5, label=site['site_name'], ax=ax)
    #ax.plot( df.index[:len(eval_DF['g'][1:])], eval_DF['g'][1:], c="k" )

    plt.ylabel(f'Nitrate (ug m-3)')
    plt.legend(loc=0, ncol=2)
    plt.tight_layout()
    plt.savefig( f"plots/{sname}" )
    plt.close()
    return


##----------------------------------------MAIN SCRIPT----------------------------------------##
def main():
    global site, rundirs
    dds = get_data_as_xr('nitrate_photol_control', year='201001')
    from CVAO_dict import CVAO_dict as d
    sys.path.append( '/users/mjr583/cvao' )
    from lowess_smoother import loess

    site = sites['CVO']
    rundirs = ['nitrate_photol_control']

    variables=['NIT', 'NITs']
    cs=['#1b9e77','#e7298a']
    f, ax = plt.subplots(1,1,figsize=(10,4))
    add_times=[] ; add_lowes=[]
    for n, rundir in enumerate(rundirs):
        ds0 = get_data_as_xr(rundir, year='20')
        ds = (ds0['SpeciesConc_NIT'] + ds0['SpeciesConc_NITs']).isel(lev=0)
        ds = ds * 1e9 * 2.54

        ds = site_data( ds, dds, lon=site['longitude'], lat=site['latitude'] )
        ds = pd.DataFrame({site['save_name']:ds.values}, index=ds['time'].values)#.resample('M').mean()
        ds = ds['2007':'2011'] 
        add_times.append( ds )

    # read_nitrate
    import read_nitrate
    df = read_nitrate.read_nitrate_dataset(convert=False)
    plot( ax, df, add_times,
          sname='cvao_nitrate.png')
       

if __name__ == "__main__":
    main()
