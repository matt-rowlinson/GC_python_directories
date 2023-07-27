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
rundir=inputs.rundir
variable=inputs.var
version=inputs.version

variable='O3'
version='12.9.3'
rundir='tropchem_merra_4x5'
var, lat, lon, lev, time = GC.get_gc_var(rundir, variable, version)
print(var.shape)
delta,interval=GC.find_timestep(time)
y = rp.find_nearest(lat, 16.9)
x = rp.find_nearest(lon, -24.9)
var_time=[]
for t in range(len(time)):
    v=var[t,-1,y,x]
    var_time.append(v)
gc=pd.DataFrame({'Value':var_time}, index=time)
print(gc)
gc=gc.resample('D').mean()

## Get Merge observations
df=CV.get_from_merge(d[variable])
df=df[gc.index[0]:gc.index[-1]]
df=df.resample(delta).mean()
print(df)
f,ax= plt.subplots(figsize=(16,4))
ax.plot(df.index, df.Value, 'k', label='CVAO')
ax.plot(gc.index, gc.Value, 'g', label='GEOS-Chem v12.9.3')
plt.ylabel('%s (%s)' % (d[variable]['abbr'], d[variable]['unit']) )

plt.gca().xaxis.set_major_locator(mdates.MonthLocator(interval=interval))

plt.ylim((7.5,52.5))
plt.legend()
plt.savefig('/users/mjr583/GC/postprocess/plots/NEW_whole_%s_timeseries_%s.png' %(rundir,variable) )
plt.close()

gcmean=gc.resample('M').mean()
gc75=gc.resample('M').quantile(.75)
gc25=gc.resample('M').quantile(.25)

df=CV.get_from_merge(d[variable])
df=df[gc.index[0]:gc.index[-1]]
dfmean=df.resample('M').mean()
df75=df.resample('M').quantile(.75)
df25=df.resample('M').quantile(.25)

f,ax= plt.subplots(figsize=(6,4))
ax.plot(dfmean.index, dfmean.Value, 'k', label='CVAO')
ax.fill_between(dfmean.index, df25.Value,df75.Value, color='grey', alpha=.3)

ax.plot(gcmean.index, gcmean.Value, 'g', label='GEOS-Chem')
ax.fill_between(gcmean.index, gc25.Value, gc75.Value, color='g', alpha=.3)

plt.ylabel('%s (%s)' % (d[variable]['abbr'], d[variable]['unit']) )
plt.gca().xaxis.set_major_locator(mdates.MonthLocator(interval=interval))

plt.ylim((7.5,52.5))
plt.legend()
plt.savefig('/users/mjr583/GC/postprocess/plots/NEWWHOLE_comparison_%s.png' %variable )
plt.close()
