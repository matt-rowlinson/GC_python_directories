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
matplotlib.use('agg')
sys.path.append('/users/mjr583/python_lib')
import GC_tools as GC
import RowPy as rp
from CVAO_dict import CVAO_dict as d
import CVAO_tools as CV
plt.style.use('seaborn-darkgrid')

inputs=GC.get_arguments()
try:
    variable=inputs.var
    print('ok')
except:
    variable='O3'
version='12.9.3'

var, lat, lon, lev, time = GC.get_gc_var(rundir='tropchem_merra_4x5', variable=variable,year='2015')
T91, lat,lon,lev,Ttime = GC.get_gc_var(rundir='tropchem_merra_4x5', variable=variable, version='12.9.1',year='2015')
nohal, lat,lon,lev,haltime = GC.get_gc_var(rundir='tropchem_merra_4x5', variable=variable, version='geos-chem',year='2015')

delta,interval=GC.find_timestep(time)
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
gc=gc.resample('MS').mean()
T91=T91.resample('MS').mean()
nohal=nohal.resample('MS').mean()

print(gc.index)

## Get Merge observations
try:
    print('here')
    df=CV.get_from_merge(d[variable])
    df=df[gc.index[0]:gc.index[-1]]
    df=df.resample('MS').mean()
    #df=df['2013']
    f,ax= plt.subplots(figsize=(12,4))
    ax.plot(df.index, df.Value, 'k', label='CVAO')
    ax.plot(gc.index, gc.Value, 'g', label='12.9.3')
    ax.plot(gc.index, T91.Value, 'lightgreen',linestyle='--', label='12.9.1')
    ax.plot(gc.index, nohal.Value, 'purple',linestyle='--', label='12.9.1_main.NOHAL')
    plt.ylabel('%s (%s)' % (d[variable]['abbr'], d[variable]['unit']) )

    plt.gca().xaxis.set_major_locator(mdates.MonthLocator(interval=1))
    print('here')
    plt.legend()
    plt.savefig('./plots/timeseries_%s.png' %variable )

except:
    f,ax= plt.subplots(figsize=(12,4))
    ax.plot(gc.index, gc.Value, 'g', label='12.9.3')
    ax.plot(gc.index, T91.Value, 'lightgreen', linestyle='--',label='12.9.1')
    ax.plot(gc.index, nohal.Value, 'purple',linestyle='--', label='12.9.1_main.NOHAL')
    plt.ylabel('%s'  %variable ) #' % (d[variable]['abbr'], d[variable]['unit']) )

    plt.gca().xaxis.set_major_locator(mdates.MonthLocator(interval=1))

    plt.legend()
    plt.savefig('./plots/timeseries_%s.png' %variable )


