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

def plot(v, month, ax, rundirs, pressure, times, sname="a.png",
                    cs=['#1b9e77','#d95f02','#7570b3','#e7298a','#66a61e','#e6ab02'],
                    style=['--','-',':']):
    for n in range(len(times)):
        for x in range(len(times[n])):
            if x==0:
                ax.scatter( times[n][x], pressure, label=rundirs[n].replace('nitrate_photol_2015_','').title(), 
                 c=cs[n], alpha=.1,zorder=2)
            else:
                ax.scatter( times[n][x], pressure, c=cs[n], alpha=.1,zorder=2)
        ax.plot( np.mean(times[n],0), pressure, 
                  label=f'{rundirs[n].replace("nitrate_photol_2015_","").title()} {month} mean',
                  c=cs[n], alpha=.76,zorder=2)
        ax.plot( np.median(times[n],0), pressure, ls='--', 
                  label=f'{rundirs[n].replace("nitrate_photol_2015_","").title()} {month} median', 
                  c=cs[n], alpha=.76,zorder=2)

    plt.ylabel('hPa')
    plt.xlabel(f'{v} {d[v]["unit"]}')
    plt.ylim( bottom=100, top=1000)
    plt.title(month)

    box = ax.get_position()
    ax.set_position([box.x0, box.y0 - box.height * 0.065,
                     box.width, box.height * 0.85])
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.075),
                  fancybox=True, shadow=True, ncol=3)

    if "NIT" in v:
        ax.set_xscale('log')
        plt.xlim( left=1e-6, right=1e3)
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
    all_sites=['CVO']#,'MAC','Grim','NYC','LA','Revin']
    for site in all_sites:
        site=sites[site]
        rundirs = ['nitrate_photol_2015_control',
                   'nitrate_photol_2015_all_scale-25']

        variables=['NIT','HNO2']#, 'C3H8','O3','SO4','HNO2','CO','NOx','NIT','NITs']
        cs=['#1b9e77','#d95f02','#7570b3','#e7298a','#66a61e','#e6ab02']
        mns=['01','02','03','04','05','06','07','08','09','10','11','12']
        months=['jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nov','dec']
        for v in variables:
            for nn, mn in enumerate(mns):
                f, ax = plt.subplots(1,1,figsize=(6,6))
                add_times=[] ; add_lowes=[]
                for n, rundir in enumerate(rundirs):
                    print( mn, n, rundir )
                    ds0 = get_data_as_xr(rundir, year=f'2017{mn}', variable=v)
                    ds0 = site_data( ds0, ds, lon=site['longitude'], lat=site['latitude'] ).values[:,:35]
                    add_times.append( ds0 )
                cs[1]=cs[3]
                plot( v, months[nn], ax, rundirs, pressure, add_times, cs=cs, 
                        sname=f'scatvert_{months[nn]}_{str(int(pressure.values[-1]))}')
        
if __name__ == "__main__":
    main()
