#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#SBATCH --job-name=VOC_fig3
#SBATCH --ntasks=1
#SBATCH --mem=2Gb
#SBATCH --partition=test
#SBATCH --time=00:05:00
#SBATCH --output=Logs/Fig.3_%A.log

import xarray as xr
import os
import glob
import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
plt.style.use('seaborn-darkgrid')

### Temporary imports
import warnings
warnings.simplefilter(action='ignore', category=pd.errors.PerformanceWarning)
warnings.simplefilter(action='ignore', category=FutureWarning)

def find_file_list(path, substrs):
    file_list =[]
    for root, directory, files in os.walk(path):
        for f in files:
            for s in substrs:
                if s in f:
                    file_list.append(os.path.join(root, f))
    file_list.sort()
    return file_list

def get_data_as_xr(rundir, year='', lev=0, version='14.0.1', variable=False):
    path=f'/mnt/lustre/users/mjr583/GC/{version}/rundirs/{rundir}'
    ds = xr.open_mfdataset( find_file_list(path, [f'SpeciesConc.{year}']), combine='by_coords' )
    if variable:
        ds = ds[f'SpeciesConc_{d[variable]["GC_name"]}'].isel(lev=lev) * float(d[variable]['scale'])
    return ds

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

def open_site_dataset_ALK4(infile):
    ALK4_df=pd.DataFrame()
    site = infile.split('_')[4]
    alk4_species = ['ic4h10','nc4h10','ic5h12','nc5h12']
    for s in alk4_species:
        f = f'/users/mjr583/noaa_flask_data/ALK4_datasets/{s}_{site}_surface-flask_1_arl_event.txt'
        with open( f ) as thefile:
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
        if ALK4_df.empty:
            ALK4_df = pd.DataFrame(data)
            ALK4_df.columns = cols
            ALK4_df['analysis_value'] = pd.to_numeric(ALK4_df['analysis_value'])
            ALK4_df.index = pd.to_datetime(ALK4_df[['year', 'month', 'day', 'hour', 'minute']])
        else:
            df = pd.DataFrame(data)
            df.columns = cols
            df.index = pd.to_datetime(df[['year', 'month', 'day', 'hour', 'minute']])
            try:
                ALK4_df['analysis_value'] = ALK4_df['analysis_value'].values + pd.to_numeric(df['analysis_value'].values)
            except:
                pass
    return ALK4_df 

def find_nearest_3d(ds,lon_value, lat_value, alt_value):
    gc_alts = pd.read_csv( '/users/mjr583/GC/info_files/GC_72_vertical_levels.csv')['Altitude (km)'] * 1e3
    x_idx=(np.abs( ds.lon - float(lon_value) )).argmin()
    y_idx=(np.abs( ds.lat - float(lat_value) )).argmin()
    z_idx=(np.abs( gc_alts.values - float(alt_value) )).argmin()
    return x_idx, y_idx, z_idx

def process_obs_and_model_data( rundirs, ds_, var ):
    first=True
    for infile in sorted( glob.glob( f'/users/mjr583/noaa_flask_data/{var}_datasets/*txt')):
        ## only process observations if data for 2016 exists
        try:
            if var=='ALK4':
                df = open_site_dataset_ALK4( infile ).loc['2016']                
            else:
                df = open_site_dataset( infile ).loc['2016']
        except:
            continue

        ## process obs
        obs = df[df.analysis_flag == '...']
        obs = pd.DataFrame( pd.to_numeric(obs.analysis_value) )
        obs = obs['analysis_value'].resample('D').mean().dropna().to_frame()
        
        ## process model data
        for ds, rundir in zip(ds_, rundirs):
            d = ds[f'SpeciesConc_{var}']
            x_idx, y_idx, z_idx = find_nearest_3d( d, df.longitude[0], df.latitude[0], df.altitude[0] )
            d = d.isel( lon=x_idx, lat=y_idx, lev=z_idx )
            d = pd.DataFrame( {f'{rundir}': d.values * 1e12}, index=d.time ).resample('D').mean()


            d = d.reindex( obs.index ).dropna()
            obs = pd.concat( [obs, d], axis=1 )
        
        if first:
            out=obs
            first=False
        else: 
            out = pd.concat( [out, obs ], axis=0 )
    return out 

def calc_stats(obs, model):
    absError= model - obs
    SE = np.square(absError)
    MSE= np.mean(SE)
    RMSE=np.round( np.sqrt(MSE), 3)
    R2 = np.round( 1. - (np.var(absError) / np.var(obs)), 3)
    nmb = 100 * (np.sum( absError - obs ) / np.sum( obs ))
    #print( 'SE', SE )
    print( 'MSE', MSE )
    print( 'nmb',  nmb )
    print( 'RMSE', RMSE)
    print( 'r2', R2 )
    return RMSE, R2

######---------------------------MAIN-SCRIPT------------------------------------------------######
def main():
    ## Set plot input criteria
    rundirs = ['geo_2x25','all_2x25']
    year    = '2016'
    
    ## Read model output as xarray datasets
    ds_=[]
    for r in rundirs:
        ds = get_data_as_xr(r, year)
        ds_.append( ds )

    ## Read observations, and process model output to match
    ethane  = process_obs_and_model_data(rundirs, ds_, 'C2H6')
    propane = process_obs_and_model_data(rundirs, ds_, 'C3H8')
    alk4    = process_obs_and_model_data(rundirs, ds_, 'ALK4')

    ## Filter 99th percentile of data
    ethane  = ethane[  ethane.analysis_value  <= np.percentile( ethane.analysis_value,  99 )]
    propane = propane[ propane.analysis_value <= np.percentile( propane.analysis_value, 99 )]
    alk4    = alk4[    alk4.analysis_value    <= np.percentile( alk4.analysis_value,    99 )]

    ethane  = ethane[  ethane[rundirs[0]]  <= np.percentile( ethane[rundirs[0]],  99 )]
    propane = propane[ propane[rundirs[0]] <= np.percentile( propane[rundirs[0]], 99 )]
    alk4    = alk4[    alk4[rundirs[0]]    <= np.percentile( alk4[rundirs[0]],    99 )]

    v = calc_stats(ethane.analysis_value, ethane['geo_2x25'] )
    v = calc_stats(ethane.analysis_value, ethane['all_2x25'] )
           
if __name__=="__main__":
    main()
