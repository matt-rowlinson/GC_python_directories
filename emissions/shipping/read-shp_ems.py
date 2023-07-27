#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#SBATCH --job-name=timeseries
#SBATCH --ntasks=1
#SBATCH --mem=10gb
#SBATCH --partition=interactive
#SBATCH --time=00:10:00
#SBATCH --output=LOGS/timeseries.log
import sys
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib
import glob
import re
import numpy as np
from netCDF4 import Dataset
matplotlib.use('agg')
sys.path.append('/users/mjr583/python_lib')
import GC_tools as GC
import RowPy as rp
from CVAO_dict import CVAO_dict as d
import CVAO_tools as CV
plt.style.use('seaborn-darkgrid')

def date_to_string(stamp):
    new_date = datetime.datetime(1950,1,1,0,0) + datetime.timedelta(float(stamp ))
    return new_date.strftime("%Y-%m-%d")

path = '/users/mjr583/GC/emissions/shipping/'
macc = Dataset('%s/MACCity-anthro_Glb_0.5x0.5_anthro_NO__monthly.nc' %path)
macc_shp = macc.variables['emiss_shp']

time=macc.variables['time']
match = re.search('\d{4}', time.units)
year = match.group(0)
macc_time = rp.days_to_datetime(time[:], start=int(year))

lon=macc.variables['lon'][:]
lat=macc.variables['lat'][:]
x = rp.find_nearest(lon, 360-24.4)
y = rp.find_nearest(lat,16.9)

cams = Dataset('%s/CAMS-GLOB-ANT_Glb_0.1x0.1_anthro_nox_v4.2_monthly.nc' %path)
cams_shp = cams.variables['shp']
time=cams.variables['time']
match = re.search('\d{4}', time.units)
year = match.group(0) if match else '2013'
cams_time = rp.hours_to_datetime(time[:], start=int(year))

lon=cams.variables['lon'][:]
lat=cams.variables['lat'][:]
x = rp.find_nearest(lon, 360-24.4)
y = rp.find_nearest(lat,16.9)

sae = np.swapaxes(rp.surface_area_earth(cams_shp.shape[2], cams_shp.shape[1]),0,1)
cams=[]
for i in range(len(cams_shp)):
    cams.append( cams_shp[i] * sae * (3600 * 24 ) * 1e-9 )
cams_shp=np.array(cams)

sae = np.swapaxes(rp.surface_area_earth(macc_shp.shape[2], macc_shp.shape[1]),0,1)
macc=[]
for i in range(len(macc_shp)):
    macc.append( macc_shp[i] * sae * (3600 * 24 ) * 1e-9 )
macc_shp=np.array(macc)

glo = np.sum(np.sum(macc_shp,1),1)
nh = np.sum(np.sum(macc_shp[:,180:,:],1),1)
nat = np.sum(np.sum(macc_shp[:, 180:, 90:180],1),1)
macc = pd.DataFrame( {'Global':glo, 'NH':nh, 'N. Atlantic': nat  }, index=macc_time)

glo = np.sum(np.sum(cams_shp,1),1)
nh = np.sum(np.sum(cams_shp[:,900:,:],1),1)
nat = np.sum(np.sum(cams_shp[:, 900:, 900:1800],1),1)
cams = pd.DataFrame( {'Global':glo, 'NH':nh, 'N. Atlantic': nat  }, index=cams_time)

plt.plot(cams.index, cams['Global'], label='CAMS')
plt.plot(macc.index, macc['Global'], label='MACCity')
plt.legend()
plt.savefig('shp_global.png')
plt.close()

plt.plot(cams.index, cams['NH'], label='CAMS')
plt.plot(macc.index, macc['NH'], label='MACCity')
plt.legend()
plt.savefig('shp_nh.png')
plt.close()

plt.plot(cams.index, cams['N. Atlantic'], label='CAMS')
plt.plot(macc.index, macc['N. Atlantic'], label='MACCity')
plt.legend()
plt.savefig('shp_NAt.png')
plt.close()
