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
import warnings
warnings.simplefilter(action='ignore', category=pd.errors.PerformanceWarning)
import glob
import numpy as np
import sys
import matplotlib.lines as mlines
import re
import random
import matplotlib.pyplot as plt
sys.path.append('/users/mjr583/python_lib')
import GC_tools as GC
import RowPy as rp
from CVAO_dict import CVAO_dict as d
plt.style.use('seaborn-darkgrid')

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

def get_site_xy_index(df):
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

def allocate_to_group(group, variables):
    if group=='SH':
        shs.append([variables[0],variables[1],variables[2],variables[3],variables[4],variables[5],variables[6],variables[7] ])
    elif group=='Polar':
        polars.append([variables[0],variables[1],variables[2],variables[3],variables[4],variables[5],variables[6],variables[7] ])
    elif group=='NAM':
        nams.append([variables[0],variables[1],variables[2],variables[3],variables[4],variables[5],variables[6],variables[7] ])
    elif group=='EUR':
        eurs.append([variables[0],variables[1],variables[2],variables[3],variables[4],variables[5],variables[6],variables[7] ])
    elif group=='AFR':
        afrs.append([variables[0],variables[1],variables[2],variables[3],variables[4],variables[5],variables[6],variables[7] ])
    elif group=='Asia':
        asias.append([variables[0],variables[1],variables[2],variables[3],variables[4],variables[5],variables[6],variables[7] ])
    mapped.append([variables[5],variables[6],variables[7]])
    return shs, polars, nams, eurs, asias, afrs, mapped

def calc_stats(obs, model):
    absError= model - obs
    SE = np.square(absError)
    MSE= np.mean(SE)
    RMSE=np.round( np.sqrt(MSE), 3)
    R2 = np.round( 1. - (np.var(absError) / np.var(obs)), 3)

    return RMSE, R2

def main():
    global gclon, gclat
    obs_=[]
    model1_=[]
    model2_=[] 
    model3_=[]
    variables=['ethane','propane']
    a=['','NH','SH']
    for aa in a:
        print( aa )
        obs_=[]
        dev_=[]
        dev2_=[] 
        dev3_=[]
        dev4_=[]
        _obs=np.arange(1) ; _dev=np.arange(1) ; _dev2=np.arange(1) ; _dev3=np.arange(1) ; _dev4=np.arange(1)
        for variable in variables:
            print( variable )
            year = '2017'
            version='13.4.0' #; variable='ethane'
            models = ['control']
            ### Get model data
            dev, gclat,gclon,lev, dev_time = GC.get_gc_var(rundir='ceds_only',variable=variable,version=version, verbose=False, year=year)
            dev2, gclat,gclon,lev, dev2_time = GC.get_gc_var(rundir='all_scaled',variable=variable,version=version, verbose=False, year=year)
            dev3, gclat,gclon,lev, dev3_time = GC.get_gc_var(rundir='geo_only',variable=variable,version=version, verbose=False, year=year)
            dev4, gclat,gclon,lev, dev4_time = GC.get_gc_var(rundir='scale_only',variable=variable,version=version, verbose=False, year=year)
            #print( dev.shape)
            ### Get observations for each flask site
            #obs=np.arange(1) ; model1=np.arange(1) ; model2=np.arange(1) ; model3=np.arange(1)
            for infile in sorted(glob.glob(f'/users/mjr583/noaa_flask_data/{variable}_datasets/*txt')):
                #print(infile)
                df = open_site_dataset(infile)
                alt_idx = get_site_altitude(df)

                if df.index[-1].year == 2016:  ## If site doesn't have data for 2016 then ignore
                    pass
                else:
                    continue
                x, y, lon, lat = get_site_xy_index(df)
                if aa=='':
                    pass
                elif aa=='NH':
                    if y < 0.:
                        continue
                elif aa=='SH':
                    if y > 0.:
                        continue

                group=rp.allocate_region(x,y)
                site=df.site_code[0]
                df = df[df.analysis_flag == '...']

                dev_var  = model_data_for_site(dev, dev_time, lon, lat, alt_idx, species=variable)
                dev2_var = model_data_for_site(dev2, dev_time, lon, lat, alt_idx, species=variable)
                dev3_var = model_data_for_site(dev3, dev_time, lon, lat, alt_idx, species=variable)
                dev4_var = model_data_for_site(dev4, dev_time, lon, lat, alt_idx, species=variable)
                
                ## Plot comparison at each individual site
                values = plot_individual_site( df, year )
                values = values.resample('D').mean().dropna()
                import datetime
                if datetime.datetime( 2016, 2, 29) in values.index:
                    values.index = values.index.to_series().replace({pd.Timestamp('2016-02-29'): pd.Timestamp('2016-02-28')})
               
                dev_var = dev_var.reindex( values.index)
                dev2_var = dev2_var.reindex( values.index)
                dev3_var = dev3_var.reindex( values.index)
                dev4_var = dev4_var.reindex( values.index)
                
                check = dev_var.values.mean() 
                if np.isfinite(check) == False:
                    continue

                _obs  = np.append( _obs, values.values )
                _dev  = np.append( _dev, dev_var.values )
                _dev2 = np.append( _dev2, dev2_var.values )
                _dev3 = np.append( _dev3, dev3_var.values )
                _dev4 = np.append( _dev4, dev4_var.values )
            

            obs_.append(  _obs )
            dev_.append(  _dev  * 1e12)
            dev2_.append( _dev2 * 1e12 )
            dev3_.append( _dev3 * 1e12)
            dev4_.append( _dev4 * 1e12 )

        # save all to csv files
        np.savetxt(f"v{version}_csvs/{aa}noaa_ethane_observations_{year}.csv", obs_[0], delimiter=",")
        np.savetxt(f"v{version}_csvs/{aa}noaa_propane_observations_{year}.csv", obs_[1], delimiter=",")

        np.savetxt(f"v{version}_csvs/{aa}ethane_ceds_only_{version}_{year}.csv", dev_[0], delimiter=",")
        np.savetxt(f"v{version}_csvs/{aa}propane_ceds_only_{version}_{year}.csv", dev_[1], delimiter=",")

        np.savetxt(f"v{version}_csvs/{aa}ethane_all_scale_{version}_{year}.csv", dev2_[0], delimiter=",")
        np.savetxt(f"v{version}_csvs/{aa}propane_all_scale_{version}_{year}.csv", dev2_[1], delimiter=",")

        np.savetxt(f"v{version}_csvs/{aa}ethane_geo_only_{version}_{year}.csv", dev3_[0], delimiter=",")
        np.savetxt(f"v{version}_csvs/{aa}propane_geo_only_{version}_{year}.csv", dev3_[1], delimiter=",")

        np.savetxt(f"v{version}_csvs/{aa}ethane_scale_only_{version}_{year}.csv", dev4_[0], delimiter=",")
        np.savetxt(f"v{version}_csvs/{aa}propane_scale_only_{version}_{year}.csv", dev4_[1], delimiter=",")
    

if __name__=="__main__":
    main()

