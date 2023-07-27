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

a, lat, lon, lev, atime = GC.get_gc_var(rundir='tropchem_merra_4x5', variable=variable)#year='2015')
b, lat,lon,lev,btime = GC.get_gc_var(rundir='fullchem_4x5_LVOCfalse', variable=variable, version='GEOS-Chem')#,year='2015')


delta,interval=GC.find_timestep(atime)
y = rp.find_nearest(lat, 16.9)
x = rp.find_nearest(lon, -24.9)
print(lon[x], lat[y])

a_time=[] ; b_time=[]
for t in range(len(atime)):
    v=a[t,0,y,x]
    a_time.append(v)

for t in range(len(btime)):
    v=b[t,0,y,x]
    b_time.append(v)

a=pd.DataFrame({'Value':a_time}, index=atime)
#a=a.tz_localize("Atlantic/Cape_Verde")a=a.tz_localize("Atlantic/Cape_Verde")

b=pd.DataFrame({'Value':b_time}, index=btime)
#b=b.tz_localize("Atlantic/Cape_Verde")
b=b[:4]
print(a)
a=a[b.index[0]:b.index[-1]]

## Get Merge observations
try:
    df=CV.get_from_merge(d[variable])
    df=df[b.index[0]:b.index[-1]]
    df=df.resample('D').mean()

    f,ax= plt.subplots(figsize=(12,4))
    ax.plot(df.index, df.Value, 'k', label='CVAO')
    ax.plot(a.index, a.Value, 'g', label='v12.9.3')
    ax.plot(b.index, b.Value, 'red',linestyle='-', label='v13.0.0')
    plt.ylabel('%s (%s)' % (d[variable]['abbr'], d[variable]['unit']) )

    plt.gca().xaxis.set_major_locator(mdates.MonthLocator(interval=6))
    plt.legend()
    plt.savefig('./plots/timeseries_v13_%s.png' %variable )

except:
    f,ax= plt.subplots(figsize=(12,4))
    ax.plot(gc.index, gc.Value, 'g', label='12.9.3')
    ax.plot(T91.index, T91.Value, 'lightgreen', linestyle='--',label='12.9.1')
    ax.plot(nohal.index, nohal.Value, 'purple',linestyle='--', label='12.9.1_main.NOHAL')
    plt.ylabel('%s'  %variable ) #' % (d[variable]['abbr'], d[variable]['unit']) )

    plt.gca().xaxis.set_major_locator(mdates.MonthLocator(interval=1))

    plt.legend()
    plt.savefig('./plots/timeseries_lvoc_nowetdep_%s.png' %variable )
