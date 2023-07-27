#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#SBATCH --job-name=jscale_time
#SBATCH --ntasks=1
#SBATCH --mem=18Gb
#SBATCH --partition=interactive
#SBATCH --time=00:30:00
#SBATCH --output=Logs/timeseries_%A.log
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
    ds = xr.open_mfdataset( find_file_list(path, [f'{collection}.{year}']), combine='by_coords' )
    try:
        ds = ds[f'SpeciesConc_{variable}'].isel(lev=0) * float(d[variable]['scale'])
    except:
        ds = ds[f'SpeciesConc_{variable}'].isel(lev=0) * 1e12
    return ds

def site_data( ds, lat=16.9, lon=-24.9, lev=0):
    x = rp.find_nearest(ds.lon, lon)#-1
    y = rp.find_nearest(ds.lat, lat)#+1
    if type(lev)==int:
        data = ds.isel( lon=x, lat=y )
    else:
        data = ds.isel( lon=x, lat=y ).isel( lev=slice(lev)).mean(dim='lev', keep_attrs=True)
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
        if variable=='O3':
            df = df / 1.9957
            df.index = pd.to_datetime( df.index, format="%d/%m/%Y %H:%M")
        else:
            df.index = pd.to_datetime( df.index, format="%Y-%m-%d %H:%M:%S")

    else:
        df = pd.read_csv( site['filepath']+site[f'{variable}_file'],index_col=0)
        df = df['Value']
    if variable=='O3':
        if site['unit_conv'] == True:
            df = df / 1.96
    df.index = pd.to_datetime( df.index, format="%Y-%m-%d %H:%M:%S")
    return df

def plot(ax, times, obs=False,
                    labels='',
                    sname="Jscale_25.png",
                    cs=['#1b9e77','#e7298a','#d95f02','#7570b3','#e7298a','#66a61e','#e6ab02'],
                    style=['-','--','--',':',':']):
    for n in range(len(times)):
        ax.plot( times[n].index, times[n][site["save_name"]], label=labels[n],
                 c=cs[n], zorder=2, ls=style[n])
    if type(obs)!=bool:
        obs.plot(c='k',zorder=1, alpha=.5, label=site['site_name'], ax=ax)

    plt.ylabel(f'{v} pptv')
    plt.legend(loc=0, ncol=2)
    plt.tight_layout()
    plt.savefig( f"plots/TEST.Ander22b.{site['save_name']}_{v}_Limits.png" )
    plt.close()
    return

def find_model_output_for_site(rundirs, versions, v, years='2017', site=sites['CVO']):
    timeseries=[] ; add_lowes=[]
    for n, rundir in enumerate(rundirs):
        #print( rundir )
        if len( versions ) == 1:
            version=versions[0]
        else:
            version=versions[n]
        ds0 = get_data_as_xr(rundir, year=years, variable=v, version=versions[0])

        ds0 = site_data( ds0, lon=site['longitude'], lat=site['latitude'] )
        ds0 = pd.DataFrame({site['save_name']:ds0.values}, index=ds0['time'].values).resample('M').mean()
        if site['save_name']=='cvao':
            ds0=ds0['2006-10-01':]
        timeseries.append( ds0 )
    return timeseries

##----------------------------------------MAIN SCRIPT----------------------------------------##
def main():
    ## Set script constants
    rundirs  = ['base_run_2015','j25-SS_run_2015','viral_run_2015', 'Ander22b_run_2015']
    lab      = ['Base','j25 (Kasibhatla et al. 2018)','Shah et al. 2022', 'Andersen et al. 2022']
    versions = ['13.1.2']
    labels   = rundirs
    variable = 'CO'
    years    = '2015'
    
    ## Get model and observational data for each GAW site
    site_list = ['Cape Verde','Mace Head', 'Tudor Hill','American Samoa','Cape Grim', 'Cape Point']
    site_     = ['CVO'       ,'MAC'      , 'TUD'       ,'SMO'           ,'Grim'     , 'Point']
    obs=[] ; sim=[] ; xy=[]
    for ss in site_:
        print( ss )
        site = sites[ss]
        mf = find_model_output_for_site(rundirs, versions, variable, years=years, site=site)
        df = load_observations(site, variable)[years].resample('M').mean()
        obs.append( df )
        sim.append( mf )
        xy.append(  f"{int(site['longitude'])}$^\circ$E, {int(site['latitude'])}$^\circ$N")
    
    ## Make the plot
    c = ['#377eb8','#984ea3','#ff7f00','#e41a1c']
    fig=plt.figure(figsize=(15,9))
    handles=[]
    for n in range(len(obs))[::-1]:
        print(n)
        ax = fig.add_subplot(3,2,n+1, aspect='auto')
        ax.plot( obs[n].index, obs[n], lw=3, label="Observations", color='k')

        for s in range(len(rundirs)):
            ln = ax.plot( sim[n][s].index, np.squeeze(sim[n][s]), lw=2, c=c[s], label=lab[s])
            if n==0:
                handles.append( ln )
        if variable == 'O3':
            ax.set_ylabel( 'O$_3$ (ppb)')
        else:
            ax.set_ylabel( f'{variable} (ppb)')
        ax.set_title(f"{site_list[n]} ({xy[n]})", weight='bold')
    
    plt.tight_layout()
    plt.subplots_adjust(top=.9)
    plt.legend( bbox_to_anchor=[1.82,1.35], fontsize=12, ncol=5, frameon=True, facecolor='w', fancybox=True )
    plt.savefig(f'plots/6xtimeseries_{variable}.png')
    plt.close()

if __name__ == "__main__":
    main()
