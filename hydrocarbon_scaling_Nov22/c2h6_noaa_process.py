#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#SBATCH --job-name=c2h6_noaa
#SBATCH --ntasks=1
#SBATCH --mem=134Gb
#SBATCH --partition=nodes
#SBATCH --time=01:45:00
#SBATCH --output=Logs/c2h6_noaa_%A.log
#SBATCH --open-mode=appendltruncate
import pandas as pd
import glob
import numpy as np
import sys
from mpl_toolkits.basemap import Basemap
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

def get_site_xy_index(df):
    x = np.round(float(df.longitude[0]),2)
    y = np.round(float(df.latitude[0]),2)
    group=rp.allocate_region(x,y)
    
    lon = rp.find_nearest(gclon, x)
    lat = rp.find_nearest(gclat, y)
    return x, y, lon, lat

def model_data_for_site( d, d_time, lon, lat,  alt_idx, species='C2H6', year=2016, label='ref'):
    d_=[]
    for t in range(len(d_time)):
        d_.append(d[t,alt_idx,lat,lon])
    ref_var = pd.DataFrame({label:d_},index=d_time)#[:year]
    ref_var = ref_var.resample('D').mean()
    ref_var.index = ref_var.index + pd.DateOffset(year=year)
    ref_var = ref_var.resample('MS').mean() * 1e12
    ref_var.index += pd.Timedelta(14, 'd')
    return ref_var

def plot_individual_site(df, GC_data, site, group, year='2016', species='C2H6'):
    df = df[df.analysis_flag == '...'][year:]
    values = pd.DataFrame(pd.to_numeric(df.analysis_value))

    f, ax = plt.subplots(figsize=(12,4))
    ax.scatter(values.index, values.analysis_value, label='_Dont', alpha=.5, color='grey')
    mmeans = values.resample('MS').mean()
    mmeans.index += pd.Timedelta(14, 'd')
    ax.plot(mmeans, 'k-', label=f'{site} ({y}N, {x}E)')#, ax=ax)

    im2 = GC_data[0].plot(ax=ax,linestyle='-', label='Base', color='#1b9e77', x_compat=True)
    im3 = GC_data[1].plot(ax=ax,linestyle='-', label='Scaled CEDS', color='#d95f02', x_compat=True)

    plt.legend()
    plt.savefig(f'plots/individual_sites/{species}_noaa_{site}_{group}_{year}.png') 
    plt.close()
    return values

def allocate_to_group(group, variables):
    if group=='SH':
        shs.append([variables[0],variables[1],variables[2],variables[3],variables[4],variables[5]])#,variables[6],variables[7] ])
    elif group=='Polar':
        polars.append([variables[0],variables[1],variables[2],variables[3],variables[4],variables[5]])#,variables[6],variables[7] ])
    elif group=='NAM':
        nams.append([variables[0],variables[1],variables[2],variables[3],variables[4],variables[5]])#,variables[6],variables[7] ])
    elif group=='EUR':
        eurs.append([variables[0],variables[1],variables[2],variables[3],variables[4],variables[5]])#,variables[6],variables[7] ])
    elif group=='AFR':
        afrs.append([variables[0],variables[1],variables[2],variables[3],variables[4],variables[5]])#,variables[6],variables[7] ])
    elif group=='Asia':
        asias.append([variables[0],variables[1],variables[2],variables[3],variables[4],variables[5]])#,variables[6],variables[7] ])
    mapped.append([variables[3],variables[4],variables[5]])
    return shs, polars, nams, eurs, asias, afrs, mapped

def plot_by_region(all_regions, region_names, region_longnames, species='C2H6'):
    for y, region in enumerate(all_regions):
        if len(region) == 0:
            continue
        fig, axes = rp.create_figure(len(region),vertical=False)
        for n in range(len(region)):
            #print(region[n])

            axes[n].scatter(region[n][0].index, region[n][0].analysis_value, label='__', alpha=.5,color='grey')
            axes[n].plot(region[n][0].resample('M').mean(), 'k-', label='%s (%sN, %sE)' 
                    %(region[n][3], region[n][5],region[n][4]) )
            axes[n].legend()
            cs = ['#377eb8','#4daf4a','#984ea3','#ff7f00']

            region[n][1].plot(ax=axes[n], linestyle='-', color=cs[1], x_compat=True, legend=False)
            region[n][2].plot(ax=axes[n], linestyle='-', color=cs[3], x_compat=True, legend=False)

        axes[0].legend()#['%s (%sN, %sE)' %(region[0][3], region[0][5], region[0][4])
                        #,'Base', 'Scaled CEDS' ])
        plt.suptitle(f'NOAA {species}: {region_longnames[y]}')
        plt.savefig(f'plots/{species}_noaa_{region_names[y]}-sites.png')
        plt.close()

def plot_sites_on_map(mapped):
    f,ax = plt.subplots(figsize=(12,6))
    m = Basemap()
    m.bluemarble()
    for i in range(len(mapped)):
        m.scatter(mapped[i][1],mapped[i][2],marker='*',color='y',s=150,zorder=10, label=mapped[i][0])
        plt.text(mapped[i][1]+random.uniform(-2,2), mapped[i][2]+random.uniform(-.5,.5),
                 mapped[i][0], weight='bold', c='y')
    plt.savefig('plots/noaa_sitemap.png')
    plt.close()
    
if __name__=="__main__":
    year = '2016'
    version='14.0.1' ; variable='ethane'
    ### Get model data
    dev, gclat,gclon,lev, dev_time = GC.get_gc_var(rundir='geo',variable=variable,version=version, verbose=True, year='2017')
    dev3, gclat,gclon,lev, dev3_time = GC.get_gc_var(rundir='scale_all',variable=variable,version=version,verbose=True, year='2017')
    
    ### Get observations for each flask site
    mapped=[]
    shs=[] ; polars=[] ; nams=[] ; eurs=[] ; asias=[] ; afrs=[]
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
        
        dev_var = model_data_for_site(dev, dev_time, lon, lat, alt_idx, label='Base')
        dev3_var = model_data_for_site( dev3, dev_time, lon, lat, alt_idx, label='Scaled CEDS')
        
        ## Plot comparison at each individual site
        values = plot_individual_site( df, [ dev_var, dev3_var], site, group, species=variable)
        shs, polars, nams, eurs, asias, afrs, mapped = allocate_to_group(group, [values,dev_var,dev3_var,site,x,y] )
    ## Now create regionally grouped plots
    all_regions=[shs,polars,nams,eurs,afrs,asias]
    region_names=['SH','Polar','NAM','EUR','AFR','Asia']
    region_longnames=['SH','Polar','N. America','Europe','North Africa','Asia']
    plot_by_region(all_regions, region_names, region_longnames, species=variable)

    plot_sites_on_map(mapped)  ### Uncomment to create map showing site locations
