#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#SBATCH --job-name=flasks_noaa
#SBATCH --ntasks=1
#SBATCH --mem=4Gb
#SBATCH --partition=test 
#SBATCH --time=00:13:00
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

def plot_individual_site(df, GC_data, site, group, year='2016', species='C2H6'):
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
    for variable in variables:
        year = '2017'
        version='14.0.1' #; variable='ethane'
        models = ['control']
        ### Get model data
        #ref, gclat,gclon,lev, ref_time = GC.get_gc_var(rundir='nitrate_photol_control',variable=variable,version=version, verbose=False, year='2017')
        dev, gclat,gclon,lev, dev_time = GC.get_gc_var(rundir='geo_2x25',variable=variable,version=version, verbose=False, year=year)
        dev2, gclat,gclon,lev, dev2_time = GC.get_gc_var(rundir='all_2x25',variable=variable,version=version, verbose=False, year=year)
        dev3, gclat,gclon,lev, dev_time = GC.get_gc_var(rundir='geo',variable=variable,version=version, verbose=False, year=year)
        dev4, gclat,gclon,lev, dev2_time = GC.get_gc_var(rundir='scale_all',variable=variable,version=version, verbose=False, year=year)
        ref = dev ; ref_time=dev_time
        ### Get observations for each flask site
        obs=np.arange(1) ; model1=np.arange(1) ; model2=np.arange(1) ; model3=np.arange(1)
        for infile in sorted(glob.glob(f'/users/mjr583/noaa_flask_data/{variable}_datasets/*txt')):
            print(infile)
            df = open_site_dataset(infile)
            alt_idx = get_site_altitude(df)

            if df.index[-1].year == 2016:  ## If site doesn't have data for 2016 then ignore
                pass
            else:
                continue
            x, y, lon, lat = get_site_xy_index(df)
            group=rp.allocate_region(x,y)
            site=df.site_code[0]
            
            ref_var = model_data_for_site(ref, ref_time, lon, lat, alt_idx, species=variable)
            dev_var = model_data_for_site(dev, dev_time, lon, lat, alt_idx, species=variable)
            dev2_var = model_data_for_site(dev2, dev2_time, lon, lat, alt_idx, species=variable)
            dev3_var = model_data_for_site(dev3, dev3_time, lon, lat, alt_idx, species=variable)
            dev4_var = model_data_for_site(dev4, dev4_time, lon, lat, alt_idx, species=variable)

            
            ## Plot comparison at each individual site
            values = plot_individual_site( df, [ref_var, dev_var, dev2_var], site, group, species=variable)
            values = values.resample('D').mean().dropna()
            import datetime
            if datetime.datetime( 2016, 2, 29) in values.index:
                values.index = values.index.to_series().replace({pd.Timestamp('2016-02-29'): pd.Timestamp('2016-02-28')})
           
            ref_var = ref_var.reindex( values.index).dropna() #* 1e12 
            dev_var = dev_var.reindex( values.index).dropna() #* 1e12
            dev2_var = dev2_var.reindex( values.index).dropna() #* 1e12
            dev3_var = dev3_var.reindex( values.index).dropna() #* 1e12
            dev4_var = dev4_var.reindex( values.index).dropna() #* 1e12

            
            values = values.reindex( dev_var.index ).dropna()

            check = dev_var.values.mean() 
            if np.isfinite(check) == False:
                continue
                #sys.exit('NaN value getting through')

            obs = np.append( obs, values.values )
            model1 = np.append( model1, ref_var.values )
            model2 = np.append( model2, dev_var.values )
            model3 = np.append( model3, dev2_var.values )
            model4 = np.append( model4, dev3_var.values )
            model5 = np.append( model5, dev4_var.values )

        
        from numpy.polynomial.polynomial import polyfit
        obs=obs[1:]
        model1=model1[1:]
        model2=model2[1:]
        model3=model3[1:]
        model4=model4[1:]
        model5=model5[1:]

        print( obs )
        print( model2)
        fltr = np.where( obs <= np.percentile( obs,99) )
        obs = obs[fltr] 
        model1 = model1[fltr] 
        model2 = model2[fltr]
        model3 = model3[fltr]
        model4 = model4[fltr]
        model5 = model5[fltr]

        model1=model2

        RMSE_1, r2_1 = calc_stats( obs, model1 )
        RMSE_2, r2_2 = calc_stats( obs, model2 )
        RMSE_3, r2_3 = calc_stats( obs, model3 )
        print( obs.shape )
        print( model2.shape )
        
        obs_.append( obs )
        model1_.append( model1 )
        model2_.append( model2 * 1e12)
        model3_.append( model3 * 1e12)
        model4_.append( model4 * 1e12)
        model5_.append( model5 * 1e12)

    
    alk4_obs = np.genfromtxt('csvs/ALK4_obs.csv', delimiter=',')
    alk4_model2 = np.genfromtxt('csvs/ALK4_GeoOnly_v14.0.csv', delimiter=',') *1e12
    alk4_model3= np.genfromtxt('csvs/ALK4_AllScale_v14.0.csv', delimiter=',') *1e12

    fig, ((ax1,ax2,ax3)) = plt.subplots(1,3,figsize=(15,5))  

    # ethane
    print( 'Ethane \n', calc_stats( obs_[0], model2_[0] ) )
    print( 'Propane \n', calc_stats( obs_[1], model2_[1] ) )
    print( 'ALK4 \n', calc_stats( alk4_obs, alk4_model2 ) )
    #sys.exit()

    Pdev = ax1.scatter( obs_[0], model2_[0], c='forestgreen', alpha=.3, zorder=1)
    a2, _, _, _ = np.linalg.lstsq(obs_[0][:,np.newaxis], model2_[0], rcond=None)
    print( 'Base a2', a2 )
    ax1.plot(obs_[0], a2 * obs_[0], '-', c='darkgreen', zorder=10, label=f'2x25 Base (z={a2[0].round(2)})')

    Pdev2 = ax1.scatter( obs_[0], model3_[0], c='orange', alpha=.3, zorder=2 )
    a2, _, _, _ = np.linalg.lstsq(obs_[0][:,np.newaxis], model3_[0], rcond=None)
    print( 'Dev a2', a2 )
    ax1.plot( obs_[0], a2 * obs_[0], '-', c='darkorange', zorder=10, label=f'2x25 Scaled CEDS emissions (z={a2[0].round(2)})')

    Pdev = ax1.scatter( obs_[0], model4_[0], c='lightgreen', alpha=.3, zorder=1)
    a2, _, _, _ = np.linalg.lstsq(obs_[0][:,np.newaxis], model4_[0], rcond=None)
    print( 'Base a2', a2 )
    ax1.plot(obs_[0], a2 * obs_[0], '--', c='darkgreen', zorder=10, label=f'4x5 Base (z={a2[0].round(2)})')

    Pdev2 = ax1.scatter( obs_[0], model5_[0], c='lightorange', alpha=.3, zorder=2 )
    a2, _, _, _ = np.linalg.lstsq(obs_[0][:,np.newaxis], model5_[0], rcond=None)
    print( 'Dev a2', a2 )
    ax1.plot( obs_[0], a2 * obs_[0], '--', c='darkorange', zorder=10, label=f'4x5 Scaled CEDS emissions (z={a2[0].round(2)})')


    ax1.set_ylabel(f'Model {d[variables[0]]["abbr"]} (ppt)')
    ax1.set_xlabel(f'Observed {d[variables[0]]["abbr"]} (ppt)')

    ax1.set_xlim(0. , np.max([obs_[0].max(), model2_[0].max()]))
    ax1.set_ylim(0. , np.max([obs_[0].max(), model2_[0].max()]))
    ax1.legend() 

    line = mlines.Line2D([0, 1], [0, 1], color='darkgrey', alpha=1.)
    transform = ax1.transAxes
    line.set_transform(transform)
    ax1.add_line(line)#, zorder=11)

    # propane
    Pdev = ax2.scatter( obs_[1], model2_[1], c='forestgreen', alpha=.3, zorder=1)
    a2, _, _, _ = np.linalg.lstsq(obs_[1][:,np.newaxis], model2_[1], rcond=None)
    ax2.plot(obs_[1], a2 * obs_[1], '-', c='darkgreen', zorder=10, label=f'Base (z={a2[0].round(2)})')

    Pdev2 = ax2.scatter( obs_[1], model3_[1], c='orange', alpha=.3, zorder=1)
    a2, _, _, _ = np.linalg.lstsq(obs_[1][:,np.newaxis], model3_[1], rcond=None)
    ax2.plot(obs_[1], a2 * obs_[1], '-', c='darkorange', zorder=10, label=f'Scaled CEDS emissions (z={a2[0].round(2)})')

    ax2.set_ylabel(f'Model {d[variables[1]]["abbr"]} (ppt)')
    ax2.set_xlabel(f'Observed {d[variables[1]]["abbr"]} (ppt)')
    ax2.legend()

    ax2.set_xlim(0. , np.max([obs_[1].max(), model2_[1].max()]))
    ax2.set_ylim(0. , np.max([obs_[1].max(), model2_[1].max()]))

    line = mlines.Line2D([0, 1], [0, 1], color='darkgrey', alpha=1.)
    transform = ax2.transAxes
    line.set_transform(transform)
    ax2.add_line(line)#, zorder=11)

    # ALK4
    Pdev = ax3.scatter( alk4_obs, alk4_model2, c='forestgreen', alpha=.05, zorder=1)
    a2, _, _, _ = np.linalg.lstsq(alk4_obs[:,np.newaxis], alk4_model2, rcond=None)
    ax3.plot(alk4_obs, a2 * alk4_obs, '-', c='darkgreen', zorder=10, label=f'Base (z={a2[0].round(2)})')

    Pdev = ax3.scatter( alk4_obs, alk4_model3, c='orange', alpha=.05, zorder=1)
    a2, _, _, _ = np.linalg.lstsq(alk4_obs[:,np.newaxis], alk4_model3, rcond=None)
    ax3.plot(alk4_obs, a2 * alk4_obs, '-', c='darkorange', zorder=10, label=f'Scaled CEDS emissions (z={a2[0].round(2)})')

    ax3.set_ylabel(f'Model ALK4 (ppt)')
    ax3.set_xlabel(f'Observed ALK4 (ppt)')
    ax3.legend()

    ax3.set_xlim(0. , np.max([np.quantile(alk4_obs,.95), np.quantile(alk4_model2, .95)]))#, np.quantile(model2, .995)]))
    ax3.set_ylim(0. , np.max([np.quantile(alk4_obs,.95), np.quantile(alk4_model2, .95)]))#, np.quantile(model2, .995)]))


    line = mlines.Line2D([0, 1], [0, 1], color='darkgrey', alpha=1.)
    transform = ax3.transAxes
    line.set_transform(transform)
    ax3.add_line(line)#, zorder=11)


    #plt.delaxes(ax4)
    plt.tight_layout()
    plt.savefig('plots/figure_5.png')
    plt.close()

if __name__=="__main__":
    main()

