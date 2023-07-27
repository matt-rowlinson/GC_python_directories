#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#SBATCH --job-name=timeseries
#SBATCH --ntasks=1
#SBATCH --mem=10gb
#SBATCH --partition=nodes
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

version='12.9.3'

rundir='tropchem_merra_4x5'
var, lat, lon, lev, time = GC.get_gc_var(rundir, variable, version, year='2016')
delta,interval=GC.find_timestep(time)
y = rp.find_nearest(lat, 16.9)
x = rp.find_nearest(lon, -24.9)
var_time=[]
for t in range(len(time)):
    v=var[t,0,y,x]
    var_time.append(v)
merra=pd.DataFrame({'Value':var_time}, index=time)
if delta=='M':
    merra=merra.resample(delta).mean()

rundir='tropchem_geosfp_4x5'
var, lat, lon, lev, time = GC.get_gc_var(rundir, variable, version, year='2016')
delta,interval=GC.find_timestep(time)
y = rp.find_nearest(lat, 16.9)
x = rp.find_nearest(lon, -24.9)
var_time=[]
for t in range(len(time)):
    v=var[t,0,y,x]
    var_time.append(v)
geosfp=pd.DataFrame({'Value':var_time}, index=time)
if delta=='M':
    geosfp=geosfp.resample(delta).mean()

rundir='benchmark'
var, lat, lon, lev, time = GC.get_gc_var(rundir, variable, version, year='2016')
delta,interval=GC.find_timestep(time)
y = rp.find_nearest(lat, 16.9)
x = rp.find_nearest(lon, -24.9)
var_time=[]
for t in range(len(time)):
    v=var[t,0,y,x]
    var_time.append(v)
bm=pd.DataFrame({'Value':var_time}, index=time)
if delta=='M':
    bm=bm.resample(delta).mean()


## Get Merge observations
df=CV.get_from_merge(d[variable])
df=df[merra.index[0]:merra.index[-1]]
df=df.resample(delta).mean()

new_index=[]
for dt in df.index:
    dtt=dt.replace(day=1)
    new_index.append(dtt)
df.index=new_index

f,ax= plt.subplots(figsize=(16,4))
ax.plot(df.index, df.Value, 'k', label='CVAO')
ax.plot(df.index, merra.Value, 'g', label='GC Merra')
ax.plot(df.index, geosfp.Value, 'lightgreen', label='GC Geos-FP')
ax.plot(df.index, bm.Value, 'r', label='GC Benchmark')

plt.ylabel('%s (%s)' % (d[variable]['abbr'], d[variable]['unit']) )

plt.gca().xaxis.set_major_locator(mdates.MonthLocator(interval=1))

#plt.ylim((7.5,52.5))
plt.legend()
plt.savefig('/users/mjr583/GC/eval_12.9.3/plots/2016_%s_timeseries_%s.png' %(rundir,variable) )
plt.close()

sys.exit()
## Add on if more regular than monthly mean output
merramean=merra.resample('M').mean()
merra75=merra.resample('M').quantile(.75)
merra25=merra.resample('M').quantile(.25)

geosfpmean=geosfp.resample('M').mean()
geosfp75=geosfp.resample('M').quantile(.75)
geosfp25=geosfp.resample('M').quantile(.25)

df=CV.get_from_merge(d[variable])
df=df[merra.index[0]:merra.index[-1]]
dfmean=df.resample('M').mean()
df75=df.resample('M').quantile(.75)
df25=df.resample('M').quantile(.25)

f,ax= plt.subplots(figsize=(6,4))
ax.plot(dfmean.index, dfmean.Value, 'k', label='CVAO')
ax.fill_between(dfmean.index, df25.Value,df75.Value, color='grey', alpha=.3)

ax.plot(merramean.index, merramean.Value, 'g', label='GC Merra')
ax.fill_between(merramean.index, merra25.Value, merra75.Value, color='g', alpha=.3)

ax.plot(geosfpmean.index, geosfpmean.Value, 'g', label='GEOS-Chem')
ax.fill_between(geosfpmean.index, geosfp25.Value, geosfp75.Value, color='g', alpha=.3)

plt.ylabel('%s (%s)' % (d[variable]['abbr'], d[variable]['unit']) )
plt.gca().xaxis.set_major_locator(mdates.MonthLocator(interval=interval))

plt.ylim((7.5,52.5))
plt.legend()
plt.savefig('/users/mjr583/GC/eval_12.9.3/plots/2016_comparison_%s.png' %variable )
plt.close()
