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
from CVAO_dict import GC_dict as d
import CVAO_tools as CV
plt.style.use('seaborn-darkgrid')

inputs=GC.get_arguments()
rundir=inputs.rundir
variable=inputs.var
version=inputs.version

version='12.9.3'

rundir='tropchem_merra_4x5'
#import getvar as GC
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
bm129=bm.resample(delta).mean()

var, lat, lon, lev, time = GC.get_gc_var(rundir, variable, version='12.4.0', year='2016')
delta,interval=GC.find_timestep(time)
y = rp.find_nearest(lat, 16.9)
x = rp.find_nearest(lon, -24.9)
var_time=[]
for t in range(len(time)):
    v=var[t,0,y,x]
    var_time.append(v)
bm=pd.DataFrame({'Value':var_time}, index=time)
bm124=bm.resample(delta).mean()

var, lat, lon, lev, time = GC.get_gc_var(rundir, variable, version='12.6.0', year='2016')
delta,interval=GC.find_timestep(time)
y = rp.find_nearest(lat, 16.9)
x = rp.find_nearest(lon, -24.9)
var_time=[]
for t in range(len(time)):
    v=var[t,0,y,x]
    var_time.append(v)
bm=pd.DataFrame({'Value':var_time}, index=time)
bm126=bm.resample(delta).mean()

var, lat, lon, lev, time = GC.get_gc_var(rundir, variable, version='12.8.0', year='2016')
delta,interval=GC.find_timestep(time)
y = rp.find_nearest(lat, 16.9)
x = rp.find_nearest(lon, -24.9)
var_time=[]
for t in range(len(time)):
    v=var[t,0,y,x]
    var_time.append(v)
bm=pd.DataFrame({'Value':var_time}, index=time)
bm128=bm.resample(delta).mean()

## Get Merge observations
#df=CV.get_from_merge(d[variable])
#df=df[merra.index[0]:merra.index[-1]]
#df=df.resample(delta).mean()

new_index=[]
for dt in bm128.index:
    dtt=dt.replace(day=1)
    new_index.append(dtt)
bm128.index=new_index

f,ax= plt.subplots(figsize=(16,4))
#ax.plot(df.index, df.Value, 'k', label='CVAO')
#ax.plot(df.index, merra.Value , '#2ca25f', label='GC Merra'  )
#ax.plot(df.index, geosfp.Value, '#e5f5f9', label='GC Geos-FP')
ax.plot(bm128.index, bm124.Value, '#fef0d9', label='GC 12.4')
ax.plot(bm128.index, bm126.Value, '#fdcc8a', label='GC 12.6')
ax.plot(bm128.index, bm128.Value, '#fc8d59', label='GC 12.8')
ax.plot(bm128.index, bm129.Value, '#d7301f', label='GC 12.9')

plt.ylabel('%s (%s)' % (d[variable]['abbr'], d[variable]['unit']) )
plt.gca().xaxis.set_major_locator(mdates.MonthLocator(interval=1))
plt.legend()
plt.savefig('/users/mjr583/GC/eval_12.9.3/plots/benchmarks_%s.png' %(variable) )
plt.close()
