#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#SBATCH --job-name=jscale_time
#SBATCH --ntasks=1
#SBATCH --mem=18Gb
#SBATCH --partition=interactive
#SBATCH --time=00:30:00
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
        if variable=='NOx':
            ds = (ds[f'SpeciesConc_NO'] + ds[f'SpeciesConc_NO2'] ) * float(d['NO']['scale'])
        else:
            ds = ds[f'SpeciesConc_{variable}'] * 1e12#float(d[variable]['scale'])
    return ds

def site_data(ds1, ds, lat=16.9, lon=-24.9):
    x = rp.find_nearest(ds.lon, lon)
    y = rp.find_nearest(ds.lat, lat)
    
    data = ds1.isel( lon=x, lat=y )#.isel( lev=slice(lev)).mean(dim='lev', keep_attrs=True)
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

def plot(ax, rundirs, times, sname="a.png",
                    cs=['#1b9e77','#d95f02','#7570b3','#e7298a','#66a61e','#e6ab02'],
                    style=['--','-',':']):
    for n in range(len(times)):
        print( n, cs[n] )
        ax.plot( times[n][site["save_name"]], times[n].index, label=rundirs[n].replace('nitrate_photol_2015_','').title(), 
                 c=cs[n], alpha=.5,zorder=2)

    plt.ylabel('hPa')
    plt.xlabel(f'{v} {d[v]["unit"]}')
    plt.xlim( left=0 )
    plt.legend(loc=0, ncol=2)
    plt.tight_layout()
    plt.gca().invert_yaxis()
    plt.savefig( f"plots/{sname}_{site['save_name']}_{v}" )
    plt.close()
    return


##----------------------------------------MAIN SCRIPT----------------------------------------##
def main():
    global d, site, rundirs, v
    ds = get_data_as_xr('nitrate_photol_control', year='201001')
    rundir, variables, version, site = get_arguments()
    from CVAO_dict import CVAO_dict as d
    sys.path.append( '/users/mjr583/cvao' )
    from lowess_smoother import loess
    
    verticals = pd.read_csv('/users/mjr583/GC/info_files/GC_72_vertical_levels.csv', delimiter=',')
    top_lev=35
    pressure  = verticals['Pressure (hPa)'][:top_lev] 
    #site = sites[site]
    all_sites=['FIREX','FIREX2','CVO','MAC','Grim','NYC','LA','Revin']
    for site in all_sites:
        site=sites[site]
        rundirs = ['nitrate_photol_2015_control',     'nitrate_photol_2015_all_scale-25',
                   'nitrate_photol_2015_all_scale-50','nitrate_photol_2015_all_scale-100',
                   'nitrate_photol_2015_all_scale-1000']

        variables=['HNO2']#,'CO', 'C3H8','O3','SO4','HNO2','CO','NOx','NIT','NITs']
        cs=['#e41a1c','#377eb8','#4daf4a', '#984ea3','#ff7f00']
        for v in variables:
            f, ax = plt.subplots(1,1,figsize=(8,8))
            add_times=[] ; add_lowes=[]
            for n, rundir in enumerate(rundirs):
                print( n, rundir )
                ds0 = get_data_as_xr(rundir, year='2018', variable=v)
                ds0 = site_data( ds0, ds, lon=site['longitude'], lat=site['latitude'] )
                ds0 = ds0.mean( dim='time' )[:top_lev]
                ds0 = pd.DataFrame({site['save_name']:ds0.values}, index=pressure)
                add_times.append( ds0 )
            plot( ax, rundirs, add_times, cs=cs, sname=f'vert_{str(int(pressure.values[-1]))}')
        
if __name__ == "__main__":
    main()
