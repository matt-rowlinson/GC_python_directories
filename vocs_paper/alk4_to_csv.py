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

def open_alk4_datasets(infile):
    infile = infile.replace('/users/mjr583/noaa_flask_data/ALK4_datasets/','')
    infile = infile.replace('_surface-flask_1_arl_event.txt','')
    site = infile[-3:]
    spec = ['ic4h10','nc4h10','ic5h12','nc5h12']
    alk_files = [f'/users/mjr583/noaa_flask_data/ALK4_datasets/ic4h10_{site}_surface-flask_1_arl_event.txt',
                 f'/users/mjr583/noaa_flask_data/ALK4_datasets/nc4h10_{site}_surface-flask_1_arl_event.txt',
                 f'/users/mjr583/noaa_flask_data/ALK4_datasets/ic5h12_{site}_surface-flask_1_arl_event.txt',
                 f'/users/mjr583/noaa_flask_data/ALK4_datasets/nc5h12_{site}_surface-flask_1_arl_event.txt',
                 ]
    #print( site )
    ALK4_df = pd.DataFrame()
    for n,f in enumerate(alk_files):
        #print( f )
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

    #print( ALK4_df['analysis_value'] )
    #sys.exit()
    return ALK4_df


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

if __name__=="__main__":
    year = '2016'
    version='13.1.2' ; variable='ALK4'
    models = ['control']
    ### Get model data
    dev, gclat,gclon,lev, dev_time = GC.get_gc_var(rundir='ceds_only',variable=variable,version=version, verbose=False, year=year)
    dev2, gclat,gclon,lev, dev2_time = GC.get_gc_var(rundir='asia-meic_scale',variable=variable,version=version, verbose=False, year=year)
    #dev3, gclat,gclon,lev, dev3_time = GC.get_gc_var(rundir='scale_only',variable=variable,version=version, verbose=False, year='2016')
    #dev4, gclat,gclon,lev, dev4_time = GC.get_gc_var(rundir='all_scaled',variable=variable,version=version, verbose=False, year='2016')
    aa = ['','NH','SH']
    for a in aa:
        print( a )
        ### Get observations for each flask site
        obs=np.arange(1) ; model1=np.arange(1) ; model2=np.arange(1) ; model3=np.arange(1)
        model4=np.arange(1) ; model5=np.arange(1)
        for infile in sorted(glob.glob(f'/users/mjr583/noaa_flask_data/{variable}_datasets/*txt')):
            #print( infile )
            if variable=='ALK4':
                df = open_alk4_datasets(infile)            
            else:
                #print(infile)
                df = open_site_dataset(infile)
            alt_idx = get_site_altitude(df)

            if df.index[-1].year == 2016:  ## If site doesn't have data for 2016 then ignore
                pass
            else:
                continue
            
            x, y, lon, lat = get_site_xy_index(df)
            if a == '':
                pass
            elif a=='NH':
                if y < 0.:
                    continue
            elif a=='SH':
                if y > 0.:
                    continue
            group=rp.allocate_region(x,y)
            site=df.site_code[0]
            
            dev_var = model_data_for_site( dev,  dev_time, lon, lat, alt_idx, species=variable) *1e12
            dev2_var = model_data_for_site(dev2, dev_time, lon, lat, alt_idx, species=variable) *1e12
            #dev3_var = model_data_for_site(dev3, dev_time, lon, lat, alt_idx, species=variable) *1e12
            #dev4_var = model_data_for_site(dev4, dev_time, lon, lat, alt_idx, species=variable) *1e12

            df = df[df.analysis_flag == '...']

            ## Plot comparison at each individual site
            values = plot_individual_site( df )
            values = values.resample('D').mean().dropna()
            import datetime
            if datetime.datetime( 2016, 2, 29) in values.index:
                values.index = values.index.to_series().replace({pd.Timestamp('2016-02-29'): pd.Timestamp('2016-02-28')})
            dev_var  = dev_var.reindex(  values.index) 
            dev2_var = dev2_var.reindex( values.index)
            #dev3_var = dev3_var.reindex( values.index) 
            #dev4_var = dev4_var.reindex( values.index)
            
            check = dev_var.values.mean() 
            if np.isfinite(check) == False:
                continue
        
            obs = np.append( obs, values.values )
            model2 = np.append( model2, dev_var.values )
            model3 = np.append( model3, dev2_var.values )
            #model4 = np.append( model4, dev3_var.values )
            #model5 = np.append( model5, dev4_var.values )
        
        np.savetxt(f"v{version}_csvs/{a}ALK4_obs_{version}_DLim.csv", obs, delimiter=",")
        np.savetxt(f"v{version}_csvs/{a}ALK4_ceds_only_{version}_DLim.csv", model2, delimiter=",")
        np.savetxt(f"v{version}_csvs/{a}ALK4_all_scale_{version}_DLim.csv", model3, delimiter=",")
        #np.savetxt(f"v{version}_csvs/{a}ALK4_scale_only_{version}_DLim.csv", model4, delimiter=",")
        #np.savetxt(f"v{version}_csvs/{a}ALK4_all_scale_{version}_DLim.csv", model5, delimiter=",")
