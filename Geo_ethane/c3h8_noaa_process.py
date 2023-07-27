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
plt.style.use('seaborn-darkgrid')

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
    #ref_var = ref_var.resample('M').mean()
    for t in range(len(dev_time)):
        dev_.append(dev[t,0,lat,lon])
    dev_var = pd.DataFrame({'dev':dev_}, index=dev_time)[:'2016']
    #dev_var = dev_var.resample('M').mean()

    df = df[df.analysis_flag == '...']['2015':]
    values = pd.DataFrame(pd.to_numeric(df.analysis_value))
    
    f, ax = plt.subplots(figsize=(12,4))
    ax.scatter(values.index, values.analysis_value, alpha=.5, color='grey')
    mmeans = values.resample('M').mean()
    ax.plot(mmeans, 'k-', label=site)

    ref_var.plot(ax=ax,linestyle=':', label='Tzompa propane', color='#1b9e77', x_compat=True)
    dev_var.plot(ax=ax,linestyle='--', label='Geological propane', color='#d95f02', x_compat=True)
    ax.legend(['%s (%sN, %sE)' %(site,y,x), 'Tzompa propane', 'CEDS + geological propane'])
    plt.savefig('plots/individual_sites/c3h8_noaa_%s_%s.png' %(site, group))
    plt.close()
    
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
        axes[n].scatter(region[n][0].index, region[n][0].analysis_value, alpha=.5,color='grey')
        axes[n].plot(region[n][0].resample('M').mean(), 'k-', label='%s (%sN, %sE)' 
                %(region[n][3], region[n][4],region[n][5]) )
        axes[n].legend()

        region[n][1].plot(ax=axes[n], linestyle=':', color='#1b9e77', x_compat=True, legend=False)
        region[n][2].plot(ax=axes[n], linestyle='--', color='#d95f02', x_compat=True, legend=False)

    axes[0].legend(['%s (%sN, %sE)' %(region[0][3], region[0][4], region[0][5])
                    ,'Tzompa propane', 'CEDS + geological propane' ])
    plt.suptitle('NOAA C3H8: %s ' %region_longnames[y])
    plt.savefig('plots/c3h8_noaa_%s-sites.png' %region_names[y])
    plt.close()

f,ax = plt.subplots(figsize=(12,6))
m = Basemap()
m.bluemarble()
for i in range(len(mapped)):
    m.scatter(mapped[i][1],mapped[i][2],marker='*',color='y',s=150,zorder=10, label=mapped[i][0])
    plt.text(mapped[i][1]+random.uniform(-2,2), mapped[i][2]+random.uniform(-.5,.5),
             mapped[i][0], weight='bold', c='y')
plt.savefig('plots/noaa_sitemap.png')
plt.close()
