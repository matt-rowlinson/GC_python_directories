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
import glob
import re
import datetime
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib
matplotlib.use('agg')
sys.path.append('/users/mjr583/python_lib')
import GC_tools as GC
import RowPy as rp
from CVAO_dict import CVAO_dict as d
import CVAO_tools as CV

plt.style.use('seaborn-darkgrid')

version='12.9.3' ; year='2016'
rundir='tropchem_merra_4x5'
path='/users/mjr583/scratch/GC/%s/rundirs/%s/OutputDir/' %(version, rundir)

from GC_tools import Cly as D

var, lat, lon, lev, time = GC.get_var_group(D, rundir=rundir, version='12.9.3', year=year)
T91, lat, lon, lev, Ttime = GC.get_var_group(D, rundir=rundir, version='12.9.1', year=year)
nohal, lat, lon, lev, haltime = GC.get_var_group(D, rundir=rundir, version='geos-chem', year=year)

y = rp.find_nearest(lat, 16.9)
x = rp.find_nearest(lon, -24.9)
print(lon[x], lat[y])

var_time=[] ; nohal_time=[] ; T91_time=[] ; var1=[]
for t in range(len(time)):
    v=var[t,0,y,x]
    var_time.append(v)

for t in range(len(Ttime)):
    v=T91[t,0,y,x]
    T91_time.append(v)

for t in range(len(haltime)):
    v=nohal[t,0,y,x]
    nohal_time.append(v)

gc=pd.DataFrame({'Value':var_time}, index=time)
T91=pd.DataFrame({'Value':T91_time}, index=Ttime)
nohal=pd.DataFrame({'Value':nohal_time}, index=haltime)

print(gc)
print(T91)
print(nohal)
gc=gc.resample('M').mean()
T91=T91.resample('M').mean()
nohal=nohal.resample('M').mean()

f,ax= plt.subplots(figsize=(12,4))
ax.plot(gc.index, gc.Value, 'g', label='12.9.3')
ax.plot(gc.index, T91.Value, 'lightgreen', linestyle='--',label='12.9.1')
ax.plot(gc.index, nohal.Value, 'purple',linestyle='--', label='12.9.1_main.NOHAL')
plt.ylabel('%s (%s)'  %( D['abbr'], D['unit'] ) ) 

plt.gca().xaxis.set_major_locator(mdates.MonthLocator(interval=1))

plt.legend()
plt.savefig('./plots/TESTtimeseries_%s.png' %D['name'] )
