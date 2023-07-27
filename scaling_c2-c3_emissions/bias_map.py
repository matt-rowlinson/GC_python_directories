#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#SBATCH --job-name=flasks_noaa
#SBATCH --ntasks=1
#SBATCH --mem=4Gb
#SBATCH --partition=interactive
#SBATCH --time=00:03:00
#SBATCH --output=Logs/flasks_noaa_%A.log
#SBATCH --open-mode=appendltruncate
import pandas as pd
import glob
import numpy as np
import sys
import matplotlib.lines as mlines
import re
import random
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
sys.path.append('/users/mjr583/python_lib')
import GC_tools as GC
import RowPy as rp
import matplotlib
from matplotlib.colors import BoundaryNorm
from CVAO_dict import CVAO_dict as d
#plt.style.use('seaborn-darkgrid')

import warnings
warnings.simplefilter(action='ignore', category=pd.errors.PerformanceWarning)

def get_site_altitude(df):
    altitude = df['altitude'][0]
    lev=pd.read_csv('/users/mjr583/GC/info_files/GC_vertical_levels.csv')['Altitude (km)'] * 1e3
    idx = rp.find_nearest( np.array(lev.values), float(altitude) )
    return idx

def open_site_dataset(infile):
    with open(infile) as thefile:
        txt = thefile.readline()
        n = [int(s) for s in txt.split() if s.isdigit()][0]
        txt=thefile.readlines()
        cols=txt[n-2].split()[2:]
        cols = [c.replace("sample_", "") for c in cols]

        grab_site = [s for s in txt[:n] if "site-code" in s][0]
        site = grab_site.split("code:",1)[1].replace(" ","")[:3]

        data=txt[n:]
        data = [w.replace(" ", ";") for w in data]
        data = [re.sub("\;+", ";", w) for w in data]

    data = [n.split(';') for n in data]
    df = pd.DataFrame(data)
    df.columns = cols
    df.index = pd.to_datetime(df[['year', 'month', 'day', 'hour', 'minute']])
    return df 

def get_site_xy_index(df, gclon, gclat):
    x = np.round(float(df.longitude[0]),2)
    y = np.round(float(df.latitude[0]),2)
    group=rp.allocate_region(x,y)
    
    lon = rp.find_nearest(gclon, x)
    lat = rp.find_nearest(gclat, y)
    return x, y, lon, lat

def model_data_for_site( d, d_time, lon, lat,  alt_idx, species='C2H6', year=2016):
    d_=[]
    for t in range(len(d_time)):
        d_.append(d[t,alt_idx,lat,lon])
    ref_var = pd.DataFrame({'ref':d_},index=d_time)#[:year]
    ref_var = ref_var.resample('D').mean()
    ref_var.index = ref_var.index + pd.DateOffset(year=year)
    #ref_var = ref_var.resample('M').mean()
    return ref_var

def plot_individual_site(df, year='2016'):
    df = df[df.analysis_flag == '...'][year:]
    values = pd.DataFrame(pd.to_numeric(df.analysis_value))
    return values


def calc_stats(obs, model):
    absError= model - obs
    SE = np.square(absError)
    MSE= np.mean(SE)
    RMSE=np.round( np.sqrt(MSE), 3)
    R2 = np.round( 1. - (np.var(absError) / np.var(obs)), 3)

    return RMSE, R2


def main():
    year = '2016'
    version='13.1.2' ; variable='propane'
    models = ['control']
    ### Get model data
    ref, gclat,gclon,lev, ref_time = GC.get_gc_var(rundir='ceds_only',variable=variable,version=version, verbose=False, year='2016')
    print( 'dev' )
    dev, gclat,gclon,lev, dev_time = GC.get_gc_var(rundir='asia-meic_scale',variable=variable,version=version, verbose=False, year='2016')
   
    ### Get observations for each flask site
    bias_base=[] ; bias_dev=[]
    lons=[] ; lats=[] 
    for infile in sorted(glob.glob(f'/users/mjr583/noaa_flask_data/{variable}_datasets/*txt')):
        print(infile)

        #sys.exit()
        df = open_site_dataset(infile)
        alt_idx = get_site_altitude(df)

        if df.index[-1].year == 2016:  ## If site doesn't have data for 2016 then ignore
            pass
        else:
            continue
        x, y, lon, lat = get_site_xy_index(df, gclon, gclat)
        group=rp.allocate_region(x,y)
        site=df.site_code[0]
        
        ref_var = model_data_for_site(ref, ref_time, lon, lat, alt_idx, species=variable) * 1e12
        dev_var = model_data_for_site(dev, dev_time, lon, lat, alt_idx, species=variable) * 1e12

        values = plot_individual_site( df )
        values = values.resample('D').mean().dropna()
        #print( values )
        #print( str(values.index[0]), type( str(values.index[0])) )
        ref_var = ref_var.reindex( values.index )#.dropna()
        dev_var = dev_var.reindex( values.index )#.dropna()
        #values = values[str(values.index) != '2016-02-29 00:00:00' ]
        #values = values.drop( values['2016-02-29'] )

        bias_base.append( np.nanmean((ref_var.values - values.values)) )
        bias_dev.append(  np.nanmean((dev_var.values - values.values)) )
        lons.append( x )
        lats.append( y )


    fig = plt.figure( figsize=(14,7) )
    ax = fig.add_subplot(111, projection=ccrs.EqualEarth(), aspect='auto')#PlateCarree(), aspect='auto')
    ax.coastlines()
    ax.set_global()
    cmap=matplotlib.cm.bwr
    bounds = np.arange(-1000,1010,10)
    norm = BoundaryNorm(bounds, cmap.N)
    
    im = ax.scatter( x=lons, y=lats, c=bias_base, norm=norm, cmap='bwr', 
                    s=200, edgecolor='k', alpha=1., transform=ccrs.PlateCarree())
    cbar_ax = fig.add_axes([0.2, 0.08, 0.6, 0.02])
    cbar = fig.colorbar(im, cax=cbar_ax, orientation='horizontal')
    cbar.ax.set_xlabel(f'Model mean - Observed mean (pptv)', size=15)
    plt.suptitle( 'Model - Observation bias: Base ', fontsize=22 )
    plt.savefig(f'plots/bias_map_{variable}_base.png' )
    plt.close()


    fig = plt.figure( figsize=(14,7) )
    ax = fig.add_subplot(111, projection=ccrs.EqualEarth(), aspect='auto')
    ax.coastlines()
    ax.set_global()
    cmap=matplotlib.cm.bwr
    bounds = np.arange(-1000,1010,10)
    norm = BoundaryNorm(bounds, cmap.N)
    
    im = ax.scatter( x=lons, y=lats, c=bias_dev, norm=norm, cmap='bwr', 
                    s=200, edgecolor='k', alpha=1., transform=ccrs.PlateCarree())
    cbar_ax = fig.add_axes([0.2, 0.08, 0.6, 0.02])
    cbar = fig.colorbar(im, cax=cbar_ax, orientation='horizontal')
    cbar.ax.set_xlabel(f'Model mean - Observed mean (pptv)', size=15)
    plt.suptitle( 'Ethane model/observation bias: Dev.scaled ', fontsize=22 )
    plt.savefig(f'plots/bias_map_{variable}_dev.png' )
    plt.close()


    fig = plt.figure( figsize=(12,12) )
    ax1 = fig.add_subplot(211, projection=ccrs.EqualEarth(), aspect='auto')
    ax2 = fig.add_subplot(212, projection=ccrs.EqualEarth(), aspect='auto')
    ax1.coastlines() ;  ax1.set_global()
    ax2.coastlines() ;  ax2.set_global()

    cmap=matplotlib.cm.bwr
    bounds = np.arange(-1000,1010,10)
    norm = BoundaryNorm(bounds, cmap.N)
    
    im1 = ax1.scatter( x=lons, y=lats, c=bias_base, norm=norm, cmap='bwr', 
                    s=200, edgecolor='k', alpha=1., transform=ccrs.PlateCarree())
    im2 = ax2.scatter( x=lons, y=lats, c=bias_dev, norm=norm, cmap='bwr', 
                    s=200, edgecolor='k', alpha=1., transform=ccrs.PlateCarree())

    NMB_base = np.round(np.sum( bias_base ) / len(bias_base),2)
    NMB_dev  = np.round(np.sum( bias_dev  ) / len(bias_dev ),2)
    
    ax1.set_title(f'Base (NMB={NMB_base} ppt)', fontsize=22)
    ax2.set_title(f'Scaled emissions (NMB={NMB_dev} ppt)', fontsize=22)
    
    plt.tight_layout()
    plt.subplots_adjust(top=.95, bottom=0.1)
    cbar_ax = fig.add_axes([0.2, 0.065, 0.6, 0.02])
    cbar = fig.colorbar(im, cax=cbar_ax, orientation='horizontal')
    cbar.ax.set_xlabel(f'Model mean - Observed mean (pptv)', size=15)
    
    plt.savefig(f'plots/bias_map_{variable}_both.png' )
    plt.close()

    sys.exit()

if __name__=="__main__":
    main()

