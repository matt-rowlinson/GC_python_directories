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
from scipy.optimize import curve_fit
from sklearn.metrics import mean_squared_error, r2_score
plt.style.use('seaborn-darkgrid')

def get_cosweights(lat):
    cosweights = np.cos(np.deg2rad(lat))
    cosweights = np.swapaxes((np.tile(cosweights,72)).reshape(72,46),0,1)
    return cosweights

def find_file_list(path, substrs):
    file_list =[]
    for root, directory, files in os.walk(path):
        for f in files:
            for s in substrs:
                if s in f:
                    file_list.append(os.path.join(root, f))
    file_list.sort()
    return file_list

def site_data(ds, lat=16.9, lon=-24.9, lev=0):
    x = rp.find_nearest(ds.lon, lon)-1
    y = rp.find_nearest(ds.lat, lat)+1
    if type(lev)==int:
        data = ds.isel( lon=x, lat=y, lev=lev)
    elif type(lev)==bool:
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
        df = df / 1.9957
        df.index = pd.to_datetime( df.index, format="%d/%m/%Y %H:%M")
    else:
        df = pd.read_csv( site['filepath']+site[f'{variable}_file'],index_col=0)
        df = df[variable]
    if site['unit_conv'] == True:
        df = df / 1.96
    df.index = pd.to_datetime( df.index, format="%Y-%m-%d %H:%M:%S")
    return df


##----------------------------------------MAIN SCRIPT----------------------------------------##
def main(): 
    #site='MAC'
    site = sites['CVO']
    rundirs  = ['base_run_2015','j25_run_2015','viral_run_2015','Ander22b_run_2015']#
    variable='NOx'
    years=['201501','201502','201507','201508']
    cs=['#1b9e77','#e7298a','#d95f02','#7570b3','#e7298a','#66a61e','#e6ab02']
    stde_=[] ; diur_=[]
    _25=[] ; _75=[]
    temp_ds=[] ; temp_25=[] ; temp_75=[]
    
    for rundir in rundirs:
        print( rundir )
        for year in years:
            #if "j25" in rundir and (year == '201507' or year== '201508'):
            #    continue
            print( year )
            path=f'/mnt/lustre/users/mjr583/GC/13.1.2/rundirs/{rundir}/OutputDir/'
            ds0 = xr.open_mfdataset( find_file_list(path, [f'SpeciesConc.{year}']), combine='by_coords' )

            cosw = get_cosweights( ds0.lat )
            ds0 = ( ds0['SpeciesConc_NO']+ds0['SpeciesConc_NO2'] ) * 1e12 * cosw

            surf = site_data( ds0,lon=site['longitude'], lat=site['latitude'] )
            ds0 = pd.DataFrame({rundir:surf.values}, index=surf.time)
            ds0.index = pd.to_datetime( ds0.index, format="%Y-%m-%d %H:%M:%S")
            ds0.index=ds0.index.tz_localize(site['tz']).tz_convert(None)
            ds1 = np.squeeze(ds0.groupby(ds0.index.hour).mean())
            std = np.squeeze(ds0.groupby(ds0.index.hour).std())
            q25 = np.squeeze(ds0.groupby(ds0.index.hour).quantile(.25))
            q75 = np.squeeze(ds0.groupby(ds0.index.hour).quantile(.75))
            temp_ds.append(ds1)
            temp_25.append(q25)
            temp_75.append(q75)

        diur_.append( np.mean( np.array(temp_ds), 0 ) )
        _25.append( np.mean( np.array(temp_25), 0 ) )
        _75.append( np.mean( np.array(temp_75), 0 ) )
    
    df = load_observations( site, 'NOx' ).dropna()
    df = df[df.index.month.isin([1,2,7,8])]
    surf_mean = df.mean()
    df_ = df.groupby(df.index.hour).mean()
    se_ = df.groupby(df.index.hour).std()

    lev=pd.read_csv('/users/mjr583/GC/info_files/GC_72_vertical_levels.csv')['Altitude (km)']
    cs = ['grey','#377eb8','#984ea3','#ff7f00','#e41a1c']

    labels = ['Base','Kasibhatla et al. 2018','Shah et al. 2022','Andersen et al. 2022']
    
    f, ax = plt.subplots(1,1, figsize=(5,5))
    
    ax.plot(    range(0,24), df_, c='k', alpha=1., label='Observations')
    ax.errorbar(range(0,24), df_, yerr=se_, c='k', label='__Observations', zorder=1,
                            capsize=3, elinewidth=1, markeredgewidth=2)

    for n in range(len(rundirs)):
        ax.plot( range( 0,24 ), diur_[n], c=cs[n+1], label=labels[n], zorder=3)
        ax.fill_between( range(0,24), _25[n], _75[n], color=cs[n+1], alpha=.2, zorder=3 )

    ax.set_ylabel( 'NOx (ppt)')
    ax.set_xlabel( 'Hour')
    ax.legend()
    
    plt.tight_layout()
    plt.savefig(f'plots/NOx_diurnal_{site["save_name"]}.png', dpi=200)
    plt.close()

if __name__ == "__main__":
    main()
