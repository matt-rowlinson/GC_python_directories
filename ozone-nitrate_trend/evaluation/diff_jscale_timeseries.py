#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#SBATCH --job-name=global_ratios
#SBATCH --ntasks=1
#SBATCH --mem=8Gb
#SBATCH --partition=interactive
#SBATCH --time=00:45:00
#SBATCH --output=Logs/diff_jTims_%A.log
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

def get_data_as_xr(rundir, year=''):
    path=f'/mnt/lustre/users/mjr583/GC/13.1.2/rundirs/{rundir}'
    ds = xr.open_mfdataset( find_file_list(path, [f'SpeciesConc.{year}']), combine='by_coords' )
    return ds

def site_data(data, ds, lat=16.9, lon=-24.9, lev=0):
    x = rp.find_nearest(ds.lon, lon)
    y = rp.find_nearest(ds.lat, lat)
    if type(lev)==int:
        data = data.isel( lev=0, lon=x, lat=y )
    else:
        data = data.isel( lon=x, lat=y ).isel( lev=slice(lev)).mean(dim='lev', keep_attrs=True)
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
    df.index = pd.to_datetime( df.index, format="%Y-%m-%d %H:%M:%S")
    return df


##----------------------------------------MAIN SCRIPT----------------------------------------##
def main():
    rundir, variables, version, site = get_arguments()
    from CVAO_dict import CVAO_dict as d
    sys.path.append( '/users/mjr583/cvao' )
    from lowess_smoother import loess
    
    thesites=['CVO','MAC','Grim']
    cs = ['#1b9e77','#d95f02','#7570b3']
    f, ax = plt.subplots(1,1,figsize=(10,4))
    for n, site in enumerate(thesites):
        site = sites[site]
        variables=['O3']#'O3','CO','ethane'] 
        for v in variables:
            data = pd.read_csv(f'csv_outputs/{site["save_name"]}_{v}.csv', index_col=0 )
            data.index= pd.to_datetime( data.index, format="%Y-%m-%d")

            J25data = pd.read_csv(f'csv_outputs/{site["save_name"]}_{v}_J25.csv', index_col=0 )
            J25data.index= pd.to_datetime( J25data.index, format="%Y-%m-%d")

            J100data = pd.read_csv(f'csv_outputs/{site["save_name"]}_{v}_J100.csv', index_col=0 )
            J100data.index= pd.to_datetime( J100data.index, format="%Y-%m-%d")

            j25  = (J25data[site["save_name"]]  - data[site["save_name"]]).dropna()
            
            #j100 = (J100data[site["save_name"]] - data[site["save_name"]]).dropna()
            j100 = - ( 100 -  (J100data[site["save_name"]] / data[site["save_name"]] * 100 ) )
            j100 = j100.dropna()
            
            ax.scatter( j100.index, j100, label='__nolegend__', c=cs[n])
            ax.plot( j100[:12].index, j100[:12], label=site["save_name"], c=cs[n])
            ax.plot( j100[12:24].index, j100[12:24], label='__nolegend__', c=cs[n])
            ax.plot( j100[24:36].index, j100[24:36], label='__nolegend__', c=cs[n])
            ax.plot( j100[36:48].index, j100[36:48], label='__nolegend__', c=cs[n])
            ax.plot( j100[48:60].index, j100[48:60], label='__nolegend__', c=cs[n])



    plt.ylabel('delta O3 (%)')
    #plt.ylabel('delta O3 (ppbv)')
    plt.legend(loc=0)
    plt.tight_layout()
    plt.savefig( f"plots/diff_{v}_Jscale_lowess_%_.png" )
    plt.close()
        
if __name__ == "__main__":
    main()
