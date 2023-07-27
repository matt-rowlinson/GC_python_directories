#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#SBATCH --job-name=shipping_CEDS
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
import datetime
import glob
import numpy as np
from netCDF4 import Dataset
matplotlib.use('agg')
sys.path.append('/users/mjr583/python_lib')
import GC_tools as GC
import RowPy as rp
from CVAO_dict import CVAO_dict as d
import CVAO_tools as CV
plt.style.use('seaborn-darkgrid')

inpath = '/mnt/lustre/groups/chem-acm-2018/earth0_data/GEOS/HEMCO/CEDS/v2020-08/'
outpath = '/users/mjr583/scratch/GC/12.9.3/rundirs/irma_025x03125_noAfrBB/Emissions/'

sae = np.swapaxes(rp.surface_area_earth(720,360),0,1) 
sec_per_year = 3600 * 24 * 365

def date_to_string(stamp):
    new_date = datetime.datetime(1950,1,1,0,0) + datetime.timedelta(float(stamp ))
    return new_date.strftime("%Y-%m-%d")

hold=[] ; times=[]
for infile in sorted(glob.glob('%s/*/NO-em-total-anthro*' %inpath)):
    fh=Dataset(infile)
    shp = np.sum(fh.variables['NO_shp'][:],0)
    #shp = shp * sae * sec_per_year * 1e-9 / 12
    hold.append(shp)

    time=fh.variables['time']
    times.append( date_to_string(time[0]) )

shp = np.array(hold)
shp = np.sum(np.sum(shp,1),1) 
df = pd.DataFrame({'shp':shp}, index=times)

f,ax = plt.subplots(figsize=(12,5))
df.plot(ax=ax,x_compat=True)
plt.savefig('plots/total_shipping_1970-2017.png')
plt.close()

f,ax = plt.subplots(figsize=(12,5))
df=df['2005':]
df.plot(ax=ax,x_compat=True)
plt.savefig('plots/total_shipping_2005-2017.png')
plt.close()
