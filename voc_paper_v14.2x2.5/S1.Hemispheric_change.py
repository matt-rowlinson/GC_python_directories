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
import datetime
import sys

#import matplotlib.pyplot as plt
#import matplotlib.dates as mdates
#import matplotlib
#matplotlib.use('agg')
#sys.path.append('/users/mjr583/python_lib')
#import GC_tools as GC
#import RowPy as rp
#from CVAO_dict import CVAO_dict as d
#import CVAO_tools as CV
#plt.style.use('seaborn-darkgrid')

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

def process_obs_and_model_data( rundirs, ds_, var, hemis ):
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
        
        ## Hemisphere split
        if hemis=='NH':
            if float(df.latitude[0]) < 0.:
                continue
        elif hemis=='SH':
            if float(df.latitude[0]) > 0.:
                continue

        ## process model data
        for ds, rundir in zip(ds_, rundirs):
            d = ds[f'SpeciesConc_{var}']
            x_idx, y_idx, z_idx = find_nearest_3d( d, df.longitude[0], df.latitude[0], df.altitude[0] )
            d = d.isel( lon=x_idx, lat=y_idx, lev=z_idx )
            d = pd.DataFrame( {f'{rundir}': d.values * 1e12}, index=d.time ).resample('D').mean()

            ## temporarily overwrite 2017 data as 2016 (including changing leap day to 28/02)
            d.index = d.index + pd.DateOffset(year=2016)
            if datetime.datetime( 2016, 2, 29) in obs.index:
                obs.index = obs.index.to_series().replace({pd.Timestamp('2016-02-29'): pd.Timestamp('2016-02-28')})
            ## 

            d = d.reindex( obs.index ).dropna()
            obs = pd.concat( [obs, d], axis=1 )
        
        if first:
            out=obs
            first=False
        else: 
            out = pd.concat( [out, obs ], axis=0 )
    return out 

######---------------------------MAIN-SCRIPT------------------------------------------------######
def main():
    ## Set plot input criteria
    rundirs = ['geo_2x25', 'all_2x25']
    year    = '2016'
    
    ## Read model output as xarray datasets
    ds_=[]
    for r in rundirs:
        ds = get_data_as_xr(r, year)
        ds_.append( ds )

    ## Read observations, and process model output to match
    NH_ethane  = process_obs_and_model_data(rundirs, ds_, 'C2H6', 'NH')
    SH_ethane  = process_obs_and_model_data(rundirs, ds_, 'C2H6', 'SH')
    NH_propane = process_obs_and_model_data(rundirs, ds_, 'C3H8', 'NH')
    SH_propane = process_obs_and_model_data(rundirs, ds_, 'C3H8', 'SH')
    NH_alk4    = process_obs_and_model_data(rundirs, ds_, 'ALK4', 'NH')
    SH_alk4    = process_obs_and_model_data(rundirs, ds_, 'ALK4', 'SH')

    ## Filter 99th percentile of data
    NH_ethane  = NH_ethane[  NH_ethane.analysis_value  <= np.percentile( NH_ethane.analysis_value,  99 )]
    NH_propane = NH_propane[ NH_propane.analysis_value <= np.percentile( NH_propane.analysis_value, 99 )]
    NH_alk4    = NH_alk4[    NH_alk4.analysis_value    <= np.percentile( NH_alk4.analysis_value,    99 )]
    
    SH_ethane  = SH_ethane[  SH_ethane.analysis_value  <= np.percentile( SH_ethane.analysis_value,  99 )]
    SH_propane = SH_propane[ SH_propane.analysis_value <= np.percentile( SH_propane.analysis_value, 99 )]
    SH_alk4    = SH_alk4[    SH_alk4.analysis_value    <= np.percentile( SH_alk4.analysis_value,    99 )]
    
    NH_ethane  = NH_ethane[  NH_ethane[rundirs[0]]  <= np.percentile( NH_ethane[rundirs[0]],  99 )]
    NH_propane = NH_propane[ NH_propane[rundirs[0]] <= np.percentile( NH_propane[rundirs[0]], 99 )]
    NH_alk4    = NH_alk4[    NH_alk4[rundirs[0]]    <= np.percentile( NH_alk4[rundirs[0]],    99 )]

    SH_ethane  = SH_ethane[  SH_ethane[rundirs[0]]  <= np.percentile( SH_ethane[rundirs[0]],  99 )]
    SH_propane = SH_propane[ SH_propane[rundirs[0]] <= np.percentile( SH_propane[rundirs[0]], 99 )]
    SH_alk4    = SH_alk4[    SH_alk4[rundirs[0]]    <= np.percentile( SH_alk4[rundirs[0]],    99 )]
    
    print( NH_ethane.shape, SH_ethane.shape )
    print( NH_propane.shape, SH_propane.shape )
    print( NH_alk4.shape, SH_alk4.shape )

    ## Plot on 3-panel scatter plot
    scat_colours = ['forestgreen','orange']
    line_colours = ['darkgreen'  ,'darkorange']
    labels=['Base','Scaled CEDS Emissions']
    alphas=[.3,.3,.05,.3,.3,.05]

    f, ((ax1,ax2,ax3),(ax4,ax5,ax6)) = plt.subplots( 2,3, figsize=(15,10) )
    axes=[ax1,ax2,ax3,ax4,ax5,ax6]
    for n, ax, var,alpha in zip(range(6), axes, [NH_ethane,NH_propane,NH_alk4,SH_ethane,SH_propane,SH_alk4],alphas):
        for n, rundir in enumerate(rundirs):
            ax.scatter( var.analysis_value, var[rundir], c=scat_colours[n], alpha=alpha, zorder=1)
            a2,_,_,_ = np.linalg.lstsq( var.analysis_value[:,np.newaxis], var[rundir], rcond=None)
            ax.plot( var.analysis_value, a2 * var.analysis_value, c=line_colours[n], label=labels[n]+" (z=%.2f)" %a2 )
            
        line  = mlines.Line2D([0,1],[0,1], color='darkgrey')
        line.set_transform(ax.transAxes)
        ax.add_line(line)
        ax.legend()
        if n >0:
            ax.set_xlim( 0., np.max([ var.analysis_value]))
            ax.set_ylim( 0., np.max([ var.analysis_value]))
    #ax2.set_xlim( 0., np.max([ var.analysis_value]))
    #ax2.set_ylim( 0., np.max([ var.analysis_value]))
    #ax3.set_xlim( 0., np.max([ var.analysis_value]))
    #ax3.set_ylim( 0., np.max([ var.analysis_value]))
    #ax4.set_xlim( 0., np.max([ var.analysis_value]))
    #ax4.set_ylim( 0., np.max([ var.analysis_value]))

    #ax5.set_xscale('log') ; ax5.set_yscale('log') 
    #ax6.set_xscale('log') ; ax6.set_yscale('log')

    ax1.set_ylabel(f'NH Model $C_2H_6$ (ppt)')
    ax1.set_xlabel(f'NH Observed $C_2H_6$ (ppt)')
    ax2.set_ylabel(f'NH Model $C_3H_8$ (ppt)')
    ax2.set_xlabel(f'NH Observed $C_3H_8$ (ppt)')
    ax3.set_ylabel(f'NH Model ALK4 (ppt)')
    ax3.set_xlabel(f'NH Observed ALK4 (ppt)')
    ax4.set_ylabel(f'SH Model $C_2H_6$ (ppt)')
    ax4.set_xlabel(f'SH Observed $C_2H_6$ (ppt)')
    ax5.set_ylabel(f'SH Model $C_3H_8$ (ppt)')
    ax5.set_xlabel(f'SH Observed $C_3H_8$ (ppt)')
    ax6.set_ylabel(f'SH Model ALK4 (ppt)')
    ax6.set_xlabel(f'SH Observed ALK4 (ppt)')



    plt.tight_layout()
    plt.savefig('plots/Sfigure_S1.png')
    plt.close()


if __name__=="__main__":
    main()
