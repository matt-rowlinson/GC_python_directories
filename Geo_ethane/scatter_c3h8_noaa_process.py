#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#SBATCH --job-name=c3h8_noaa
#SBATCH --ntasks=1
#SBATCH --mem=140Gb
#SBATCH --partition=nodes
#SBATCH --time=01:45:00
#SBATCH --output=Logs/c3h8_noaa_%A.log
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
import matplotlib.lines as mlines
plt.style.use('seaborn-darkgrid')

from sklearn import datasets, linear_model
from sklearn.metrics import mean_squared_error, r2_score

print('Load reference GC')
ref, gclat,gclon,lev, ref_time = GC.get_gc_var(rundir='fullchem_4x5_LVOCfalse',variable='propane',version='GEOS-Chem', verbose=False)
print('Load development GC')
dev, gclat,gclon,lev, dev_time = GC.get_gc_var(rundir='geo_ethane_propane',variable='propane',version='GEOS-Chem',verbose=False)

mapped=[]
shs=[] ; polars=[] ; nams=[] ; eurs=[] ; asias=[] ; afrs=[]
for infile in sorted(glob.glob('propane_datasets/*txt')):
    print(infile)

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

    if df.index[-1].year == 2016:
        pass
    else:
        continue
    
    site=df.site_code[0]
    x = np.round(float(df.longitude[0]),2)
    y = np.round(float(df.latitude[0]),2)
    group=rp.allocate_region(x,y)
    
    lon = rp.find_nearest(gclon, x)
    lat = rp.find_nearest(gclat, y)
    
    ref_=[] ; dev_=[]
    for t in range(len(ref_time)):
        ref_.append(ref[t,0,lat,lon])
    ref_var = pd.DataFrame({'ref':ref_},index=ref_time)[:'2016']
    for t in range(len(dev_time)):
        dev_.append(dev[t,0,lat,lon])
    dev_var = pd.DataFrame({'dev':dev_}, index=dev_time)[:'2016']

    df = df[df.analysis_flag == '...']['2015':]
    values = pd.DataFrame(pd.to_numeric(df.analysis_value))
    
    if group=='SH':
        shs.append([values,ref_var, dev_var, site, x, y ])
    elif group=='Polar':
        polars.append([values,ref_var, dev_var, site, x, y ])
    elif group=='NAM':
        nams.append([values,ref_var,dev_var, site, x, y])
    elif group=='EUR':
        eurs.append([values,ref_var,dev_var, site, x, y])
    elif group=='AFR':
        afrs.append([values,ref_var, dev_var, site, x, y])
    elif group=='Asia':
        asias.append([values,ref_var, dev_var, site, x, y])
    mapped.append([site,x,y])

all_regions=[shs,polars,nams,eurs,afrs,asias]
region_names=['SH','Polar','NAM','EUR','AFR','Asia']
region_longnames=['SH','Polar','N. America','Europe','North Africa','Asia']
for y, region in enumerate(all_regions):
    fig, axes = rp.create_figure(len(region),vertical=False)
    for n in range(len(region)):
        obs = region[n][0].resample('D').mean()
        ref = region[n][1].resample('D').mean().reindex(obs.index)
        dev = region[n][2].resample('D').mean().reindex(obs.index)

        axes[n].scatter( obs , ref , color='#1b9e77')
        axes[n].scatter( obs , dev , color='#d95f02')

        Max=np.nanmax([np.nanquantile(obs, 0.9), np.nanquantile(ref, .9), np.nanquantile(dev, .9)])
        Min=np.nanmin([np.nanquantile(obs, 0.1), np.nanquantile(ref, .1), np.nanquantile(dev, .1)])
        axes[n].set_xlim(0.,Max)
        axes[n].set_ylim(0.,Max)
        axes[n].set_xlabel('Observed C3H8 (ppb)')
        axes[n].set_ylabel('Modelled C3H8 (ppb)')

        line = mlines.Line2D([0, 1], [0, 1], color='red', alpha=0.2)
        transform = axes[n].transAxes
        line.set_transform(transform)
        axes[n].add_line(line)

        XX = np.arange(len(obs[obs.columns[0]]))
        idx = np.isfinite(obs[obs.columns[0]])
        X = obs[obs.columns[0]][idx]
        Y = ref[ref.columns[0]][idx]
        Z = dev[dev.columns[0]][idx]
       
        z = np.polyfit(X, Y, 1)
        p = np.poly1d(z)
        axes[n].plot(obs[idx],p(X),"#1b9e77")

        z = np.polyfit(X, Z, 1)
        p = np.poly1d(z)
        axes[n].plot(obs[idx],p(X),"#d95f02")

        absError= X - Y
        SE = np.square(absError)
        MSE= np.mean(SE)
        REFRMSE=np.round( np.sqrt(MSE), 3)
        REFR2 = np.round( 1. - (np.var(absError) / np.var(X)), 3)

        absError= X - Z
        SE = np.square(absError)
        MSE= np.mean(SE)
        DEVRMSE=np.round( np.sqrt(MSE), 3)
        DEVR2 = np.round( 1. - (np.var(absError) / np.var(X)), 3)
        print(region[n][3],'\n RMSE=', REFRMSE, DEVRMSE, 'ppb \n R2=', REFR2, DEVR2, '%')

        axes[n].set_title('%s (%sN, %sE)' %(region[n][3], region[n][4], region[n][5]) )# , pad=-20.6  )

    axes[0].legend(['Tzompa propane', 'CEDS + geological propane' ], bbox_to_anchor=(2., 1.11), ncol=3)
    plt.subplots_adjust(hspace=0.3) 
    plt.suptitle('NOAA C3H8: %s ' %region_longnames[y])
    plt.savefig('working_plots/c3h8_noaa_%s-sites.png' %region_names[y])
    plt.close()
