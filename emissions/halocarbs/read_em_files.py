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
import numpy as np
from netCDF4 import Dataset
matplotlib.use('agg')
sys.path.append('/users/mjr583/python_lib')
import GC_tools as GC
import RowPy as rp
from CVAO_dict import CVAO_dict as d
import CVAO_tools as CV
#plt.style.use('seaborn-darkgrid')

inpath = '/mnt/lustre/groups/chem-acm-2018/earth0_data/GEOS/HEMCO/CEDS/v2020-08/2017/'

for infile in sorted(glob.glob('%s/C2H6-*anthro*' %inpath)):
    print(infile)
    fh = Dataset(infile)
    lat=fh.variables['lat'][:]
    lon=fh.variables['lon'][:]
    keys=list(fh.variables.keys())
    em=0
    for key in keys[3:]:
        em += fh.variables[key][:]
em = np.sum(em, 0)
sae = np.swapaxes(rp.surface_area_earth(720, 360),0,1)
em = em * sae * (3600*24)*365/12 * 1e-9
for ilon in range(len(lon)):
    for ilat in range(len(lat)):
        if -10.2 < lon[ilon] < 2.:
            if 49.7 < lat[ilat] < 60.6:
                pass
            else:
                em[ilat, ilon] = np.nan
        else:
            em[ilat, ilon] = np.nan
print(np.nansum(em))
CEDS=em


inpath = '/mnt/lustre/groups/chem-acm-2018/earth0_data/GEOS/HEMCO/C2H6_2010/v2019-06/'
hold=np.zeros((12,91,144))
for n, infile in enumerate(sorted(glob.glob('%s/C2H6*anth*2010*2x25*' %inpath))):
    print(infile)
    fh = Dataset(infile)
    Tlat=fh.variables['lat'][:]
    Tlon=fh.variables['lon'][:]
    keys=list(fh.variables.keys())
    
    anthro = fh.variables['ANTHR_C2H6'][:]
    biof = fh.variables['BIOFUEL_C2H6'][:]
    total=anthro+biof
    hold[n]=total

factor = 30.07 / (2*12.011)
em = hold * factor
sae = np.swapaxes(rp.surface_area_earth(144, 91),0,1)
new=[]
for n in range(12):
    x = em[n] * sae * (3600 * 24 * 30.1) * 1e-9
    new.append(x)
em = np.sum(np.array(new),0)
print(em.sum())

print(em.shape)
for ilon in range(len(Tlon)):
    for ilat in range(len(Tlat)):
        if -10.2 < Tlon[ilon] < 2.:
            if 49.7 < Tlat[ilat] < 60.6:
                pass
            else:
                em[ilat, ilon] = np.nan
        else:
            em[ilat, ilon] = np.nan
print(np.nansum(em))
Tzompa=em

f, (ax1,ax2) = plt.subplots(1,2, figsize=(12,4))
X,Y = np.meshgrid(lon, lat)
m = rp.get_basemap(lllat=45, urlat=65, lllon=-16, urlon=8, ax=ax1)
m.pcolormesh(X,Y, CEDS)
X,Y = np.meshgrid(Tlon, Tlat)
m = rp.get_basemap(lllat=45, urlat=65, lllon=-16, urlon=8, ax=ax2)
m.pcolormesh(X,Y, Tzompa)

plt.savefig('BOTH_UK_C2H6.png')
plt.close()

sys.exit()

