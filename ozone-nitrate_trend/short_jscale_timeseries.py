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
        if variable == 'NOx':
            ds = ( ds[f'SpeciesConc_NO'] + ds[f'SpeciesConc_NO2'] ).isel(lev=0) * 1e12#float(d[variable]['scale'])
        else:
            ds = ds[f'SpeciesConc_{variable}'].isel(lev=0) * float(d[variable]['scale'])
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
        df = pd.read_csv( site['filepath'], index_col=0, dtype={"Airmass":str, "New_Airmass":str})
        #df.index = pd.to_datetime( df.index, format='%Y-%m-%d %H:%M:%S' )
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
        #times[n] = times[n].resample('M').mean()
        ax.plot( times[n].index, times[n][site["save_name"]], label=rundirs[n].replace('nitrate_photol_','').title(), 
                 c=cs[n], zorder=2, ls=style[n])
        #ax.plot( times[n].index, lowes[n]['g'][1:], c=cs[n])#, ls=style[n])
    if type(obs)!=bool:
        obs.plot(c='k',zorder=1, alpha=.5, label=site['site_name'], ax=ax)
    #ax.plot( df.index[:len(eval_DF['g'][1:])], eval_DF['g'][1:], c="k" )

    plt.ylabel(f'{d[v]["abbr"]}  ({d[v]["unit"]})')
    plt.legend(loc=0, ncol=2)
    plt.tight_layout()
    plt.savefig( f"plots/{site['save_name']}_{v}_{sname}" )
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

    site = sites[site]

    rundirs = ['nitrate_photol_control','nitrate_photol_all_Jscale-25',
            'nitrate_photol_2015_all_scale-50','nitrate_photol_2015_all_scale-100', 'nitrate_photol_2015_all_scale-1000' ]
    sname='j1000'

    variables=['NOx','O3']
    cs=['#1b9e77','#e7298a','#d95f02','#7570b3','#e7298a','#66a61e','#e6ab02']
    for v in variables:
        print( v )
        f, ax = plt.subplots(1,1,figsize=(10,4))
        add_times=[] ; add_lowes=[]
        for n, rundir in enumerate(rundirs):
            print(rundir)
            ds0 = get_data_as_xr(rundir, year='2019', variable=v)
          
            ds0 = site_data( ds0, ds, lon=site['longitude'], lat=site['latitude'] )
            ds0 = pd.DataFrame({site['save_name']:ds0.values}, index=ds0['time'].values)#.resample('M').mean()
           
            if site['save_name']=='cvao':
                ds0=ds0['2006-10-01':]
            ds0=ds0[:'2019']
            
            data = pd.DataFrame( {"XValue" : np.arange(len(ds0[site['save_name']].values)), "YValue" : ds0[site['save_name']]},
                        index=ds0.index  )
            add_times.append( ds0 )
        
        try:
            df = load_observations(site, v)[:]#.resample("M").mean()
            df=df['2019']
            plot( ax, add_times, sname=sname, obs=df)
        except:
            plot(ax, add_times, sname=sname)

        
if __name__ == "__main__":
    main()
