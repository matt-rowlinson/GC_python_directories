#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#SBATCH --job-name=timeseries
#SBATCH --ntasks=1
#SBATCH --mem=10gb
#SBATCH --partition=interactive
#SBATCH --time=00:10:00
#SBATCH --output=LOGS/timeseries.log
import sys
sys.path.append('/users/mjr583/python_lib')
import RowPy as rp
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib
import read as R
from matplotlib.colors import LogNorm
import numpy as np
from netCDF4 import Dataset
matplotlib.use('agg')
plt.style.use('seaborn-darkgrid')

df=pd.read_csv('NEW_EMEP_NMVOCs_2000-2018_by_country.csv', index_col=0)

years = range(2000, 2020, 1)
years = [ str(xx) for xx in years ]
x = pd.to_datetime( years, format="%Y" )

inpath = '/users/mjr583/GC/c2_c3_emissions/temp_files/CO_2019_GRID_1990_to_2017.nc'
fh = Dataset(inpath)
time=fh.variables['time'][:]
Elat = fh.variables['lat'][:]
Elon = fh.variables['lon'][:]
emep_em = fh.variables['sumallsectors'][:]
emep_em = emep_em[10:]#.sum()

v21_voc = np.load('temp_files/v21_voc.npy')
fh = Dataset('/mnt/lustre/groups/chem-acm-2018/earth0_data/GEOS/HEMCO/CEDS/v2021-06/2017/CO-em-anthro_CMIP_CEDS_2017.nc')
lat = fh.variables['lat'][:]
lon = fh.variables['lon'][:]

ceds_em=[]
for i in range(len(v21_voc)):
    ceds_em.append( R.regional( v21_voc[i,:,:], lat, lon, region='Europe') )
Max = np.nanmax(emep_em)
maX = np.nanmax(ceds_em)
Max = np.max((Max, maX))
for n, y in enumerate(years):
    try:
        ax = np.nanmax(emep_em[n])
        maX = np.nanmax(ceds_em[n])
        #print(y, ':', np.round(np.max((ax, maX)), 2) )
        f,(ax,ax1) =plt.subplots(2,1, figsize=(8,8))
        m=rp.get_basemap(lllat=Elat.min()-5, urlat=Elat.max()+5, lllon=Elon.min()-5, urlon=Elon.max()+5, ax=ax)
        X,Y=np.meshgrid(Elon,Elat)
        m.pcolormesh(X,Y, emep_em[n], cmap='viridis', vmax=Max, norm=LogNorm())

        m=rp.get_basemap(lllat=Elat.min()-5, urlat=Elat.max()+5, lllon=Elon.min()-5, urlon=Elon.max()+5, ax=ax1)
        X,Y=np.meshgrid(lon,lat)
        m.pcolormesh(X,Y, ceds_em[n], cmap='viridis', vmax=Max, norm=LogNorm())

        f.text(0.16, .845, f'EMEP {y}: %s Tg' %np.round(np.nansum(emep_em[n]),2),
                                    fontsize=14, bbox=dict(facecolor='w', alpha=1.) )
        f.text(0.16, .425, f'CEDS {y}: %s Tg' %np.round(np.nansum(ceds_em[n]),2),
                                    fontsize=14, bbox=dict(facecolor='w', alpha=1.) )

        plt.savefig(f'plots/CO_emep_ceds_map_{y}.png')
        plt.close()
    except:
        print("No data for {y}")
        pass


ceds_em=[]
for i in range(len(v21_voc)):
    ceds_em.append( R.regional_totals( v21_voc[i,:,:], lat, lon, region='Europe') )

emep_em = np.sum(np.sum(emep_em, 1),1)
print(emep_em)
print(ceds_em)
f,ax=plt.subplots()
ax.plot( x[:len(ceds_em)], ceds_em, label='CEDS')
ax.plot( x[:len(emep_em)], emep_em, label='EMEP')
#ax.plot( x[:len(df['Total (Tg)'])], df['Total (Tg)'].values, label='EMEP by country')
plt.legend()
plt.savefig('plots/CO_emep.png')
plt.close()
