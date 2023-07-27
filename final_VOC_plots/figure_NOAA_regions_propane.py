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
import warnings
warnings.filterwarnings("ignore")

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

def process_obs_and_model_data( rundirs, ds_, var, region ):
    first=True
    xy=[]
    for n, infile in enumerate(sorted( glob.glob( f'/users/mjr583/noaa_flask_data/{var}_datasets/*txt'))):
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
        if region=='n-america':
            if 10. < float(df.latitude[0]) < 90. and -180. < float(df.longitude[0]) < -85.:
                pass
            elif 55. < float(df.latitude[0]) < 90. and -100. < float(df.longitude[0]) < -30.: 
                pass
            else:
                continue
        if region=='europe':
            if 36.5 < float(df.latitude[0]) < 70. and -27. < float(df.longitude[0]) < 50.:
                pass
            elif 20. < float(df.latitude[0]) < 40. and 0. < float(df.longitude[0]) < 20.:
                pass
            else:
                continue
        if region=='asia':
            if 5. < float(df.latitude[0]) < 80. and 50. < float(df.longitude[0]) < 150.:
                pass
            else:
                continue
        if region=='pacific':
            if -60. < float(df.latitude[0]) < 70. and (-140. > float(df.longitude[0]) or float(df.longitude[0]) > 130.):
                pass
            else:
                continue
        if region=='atlantic':
            if 0. < float(df.latitude[0]) < 50. and -70. < float(df.longitude[0]) < -10.:
                pass
            else:
                continue
        elif region=='SH':
            if float(df.latitude[0]) < 0.:
                pass
            else:
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
            ll = pd.DataFrame({'lon':df.longitude[0], 'lat':df.latitude[0], 'region':region},
                    index=[0])
        else: 
            out = pd.concat( [out, obs ], axis=0 )

            ll.loc[n] = [df.longitude[0]] + [df.latitude[0]] + [region]
    #ll.to_csv(f'lat-lon_{region}.csv')  
    return out 

def func(p, x):
    a, b = p
    return a * x + b

def get_odr(x,y):
    from scipy import odr
    quad_model = odr.Model(func) 

    data = odr.RealData( np.log10(x),np.log10(y) )
    odr  = odr.ODR( data, quad_model, beta0=[0., 1.,] ) 
    out  = odr.run()
    popt = out.beta
    perr = out.sd_beta
    #print( popt )

    for  i in range(len(popt)):
        nstd=10
        popt_up = popt + nstd * perr
        popt_dw = popt - nstd * perr
        x_fit = np.linspace( np.log10(x).min(), np.log10(x).max(), len(x) )
        x_fit = np.linspace( 0, 4400, len(x) )

        fit = func( popt, x_fit )
        return x_fit, fit, popt[0] 

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
    EU_propane  = process_obs_and_model_data(rundirs, ds_, 'C3H8', 'europe')
    NA_propane  = process_obs_and_model_data(rundirs, ds_, 'C3H8', 'n-america')
    AS_propane  = process_obs_and_model_data(rundirs, ds_, 'C3H8', 'asia')
    PA_propane  = process_obs_and_model_data(rundirs, ds_, 'C3H8', 'pacific')
    AT_propane  = process_obs_and_model_data(rundirs, ds_, 'C3H8', 'atlantic')
    SH_propane  = process_obs_and_model_data(rundirs, ds_, 'C3H8', 'SH')
    
    ## Filter 99th percentile of data
    SH_propane  = SH_propane[  SH_propane.analysis_value  <= np.percentile( SH_propane.analysis_value,  99 )]
    NA_propane  = NA_propane[  NA_propane.analysis_value  <= np.percentile( NA_propane.analysis_value,  95 )]
    EU_propane  = EU_propane[  EU_propane.analysis_value  <= np.percentile( EU_propane.analysis_value,  99 )]
    AS_propane  = AS_propane[  AS_propane.analysis_value  <= np.percentile( AS_propane.analysis_value,  99 )]
    PA_propane  = PA_propane[  PA_propane.analysis_value  <= np.percentile( PA_propane.analysis_value,  99 )]
    AT_propane  = AT_propane[  AT_propane.analysis_value  <= np.percentile( AT_propane.analysis_value,  99 )]
    
    NA_propane  = NA_propane[  NA_propane[rundirs[0]]  <= np.percentile( NA_propane[rundirs[0]],  99 )]
    EU_propane  = EU_propane[  EU_propane[rundirs[0]]  <= np.percentile( EU_propane[rundirs[0]],  95 )]
    AS_propane  = AS_propane[  AS_propane[rundirs[0]]  <= np.percentile( AS_propane[rundirs[0]],  99 )]
    PA_propane  = PA_propane[  PA_propane[rundirs[0]]  <= np.percentile( PA_propane[rundirs[0]],  99 )]
    AT_propane  = AT_propane[  AT_propane[rundirs[0]]  <= np.percentile( AT_propane[rundirs[0]],  99 )]
    SH_propane  = SH_propane[  SH_propane[rundirs[0]]  <= np.percentile( SH_propane[rundirs[0]],  99 )]
    
    #print( NA_propane.shape, EU_propane.shape,AS_propane.shape, PA_propane.shape,AT_propane.shape, SH_propane.shape )
    ## Plot on 3-panel scatter plot
    scat_colours = ['blue','orange'][::-1]
    line_colours = ['blue','darkorange'][::-1]
    labels=['Base','Scaled Emissions']
    titles=['North America','Europe','Asia','Pacific Ocean','Atlantic Ocean','Southern Hemisphere']
    alphas=[.2,.2,.2,.2,.2,.2]
    
    from sklearn.metrics import mean_squared_error
    
    ## Log-log plot
    f, ((ax1,ax2,ax3),(ax4,ax5,ax6)) = plt.subplots( 2,3, figsize=(15,10),sharey=True, sharex=True )
    axes=[ax1,ax2,ax3,ax4,ax5,ax6]
    ps = [NA_propane, EU_propane,AS_propane, PA_propane, AT_propane, SH_propane]
    for nn, ax, var,alpha in zip(range(6), axes, ps,alphas):
        for n, rundir in enumerate(rundirs[:1]):
            ax.scatter( var.analysis_value, var[rundir], c=scat_colours[n], alpha=alpha, zorder=1)
            rms = mean_squared_error(var.analysis_value, var[rundir], squared=False)
            nmb = 100 * np.sum(var[rundir] - var.analysis_value) / np.sum(var.analysis_value)
            print( titles[nn], rundir, nmb.round(2), '%'  )

            x_fit, fit, a2 = get_odr( var.analysis_value, var[rundir] )
            ax.plot( 10**x_fit, 10**fit, c=line_colours[n], zorder=200,
                        label=labels[n]+" (NMB=%.0f%%)" %nmb)
            rms = mean_squared_error(var.analysis_value, var[rundir], squared=False)

        line  = mlines.Line2D([0,1],[0,1], color='darkgrey', zorder=1)
        line.set_transform(ax.transAxes)
        ax.add_line(line)
        ax.legend(fontsize=12)
        ax.set_xscale('log') ; ax.set_yscale('log') 
        ax.set_title(titles[nn],weight='bold', fontsize=16)

        ax.set_xlim(0.45, 4000) ; ax.set_ylim(0.45, 4000)  
        print( nn ) 
        if nn==0 or nn==3:
            ax.set_ylabel(f'Model $C_2H_6$ (ppt)', fontsize=16)
        if nn>2:
            ax.set_xlabel(f'Observed $C_2H_6$ (ppt)', fontsize=16)
        ax.tick_params(axis='both', which='major', labelsize=14)

    plt.tight_layout()
    plt.savefig('plots/NOAA_regions_propane_base.png')
    plt.close()

    ## Log-log plot
    f, ((ax1,ax2,ax3),(ax4,ax5,ax6)) = plt.subplots( 2,3, figsize=(15,10),sharey=True, sharex=True )
    axes=[ax1,ax2,ax3,ax4,ax5,ax6]
    ps = [NA_propane, EU_propane,AS_propane, PA_propane, AT_propane, SH_propane]
    for nn, ax, var,alpha in zip(range(6), axes, ps,alphas):
        for n, rundir in enumerate(rundirs):
            ax.scatter( var.analysis_value, var[rundir], c=scat_colours[n], alpha=alpha, zorder=1)
            rms = mean_squared_error(var.analysis_value, var[rundir], squared=False)
            nmb = 100 * np.sum(var[rundir] - var.analysis_value) / np.sum(var.analysis_value)
            print( titles[nn], rundir, nmb.round(2), '%'  )

            x_fit, fit, a2 = get_odr( var.analysis_value, var[rundir] )
            ax.plot( 10**x_fit, 10**fit, c=line_colours[n], zorder=200,
                        label=labels[n]+" (NMB=%.0f%%)" %nmb)
            rms = mean_squared_error(var.analysis_value, var[rundir], squared=False)

        line  = mlines.Line2D([0,1],[0,1], color='darkgrey', zorder=1)
        line.set_transform(ax.transAxes)
        ax.add_line(line)
        ax.legend(fontsize=12)
        ax.set_xscale('log') ; ax.set_yscale('log') 
        ax.set_title(titles[nn],weight='bold', fontsize=16)

        ax.set_xlim(0.45, 4000) ; ax.set_ylim(0.45, 4000)  
        print( nn ) 
        if nn==0 or nn==3:
            ax.set_ylabel(f'Model $C_2H_6$ (ppt)', fontsize=16)
        if nn>2:
            ax.set_xlabel(f'Observed $C_2H_6$ (ppt)', fontsize=16)
        ax.tick_params(axis='both', which='major', labelsize=14)

    plt.tight_layout()
    plt.savefig('plots/NOAA_regions_propane.png')
    plt.close()



if __name__=="__main__":
    main()
