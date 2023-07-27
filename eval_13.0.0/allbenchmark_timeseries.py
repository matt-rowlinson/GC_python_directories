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
delta='MS'
if delta=='MS':
    merra=merra.resample(delta).mean()

var, lat, lon, lev, time = GC.get_gc_var(rundir, variable, version='12.9.1', year='2016')
delta,interval=GC.find_timestep(time)
y = rp.find_nearest(lat, 16.9)
x = rp.find_nearest(lon, -24.9)
var_time=[]
for t in range(len(time)):
    v=var[t,0,y,x]
    var_time.append(v)
bm=pd.DataFrame({'Value':var_time}, index=time)
bm1291=bm.resample('MS').mean()


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
bm1293=bm.resample(delta).mean()

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

var, lat, lon, lev, time = GC.get_gc_var(rundir='fullchem_hourly_default', variable=variable, version='GEOS-Chem', year='2016')
delta,interval=GC.find_timestep(time)
y = rp.find_nearest(lat, 16.9)
x = rp.find_nearest(lon, -24.9)
var_time=[]
for t in range(len(time)):
    v=var[t,0,y,x]
    var_time.append(v)
bm=pd.DataFrame({'Value':var_time}, index=time)
bm13=bm.resample('MS').mean()


## Get Merge observations
df=CV.get_from_merge(d[variable])
df=df[merra.index[0]:merra.index[-1]]
df=df.resample('MS').mean()

new_index=[]
for dt in df.index:
    dtt=dt.replace(day=1)
    new_index.append(dtt)
df.index=new_index

f,ax= plt.subplots(figsize=(16,4))
ax.plot(df.index, df.Value, 'k', label='CVAO')
ax.plot(df.index, bm124.Value, '#fef0d9', label='v12.4.0')
ax.plot(df.index, bm126.Value, '#fdcc8a', label='v12.6.0')
ax.plot(df.index, bm128.Value, '#fc8d59', label='v12.8.0')
ax.plot(df.index, bm1291.Value, 'green', label='v12.9.1')
ax.plot(df.index, bm1293.Value, '#d7301f', label='v12.9.3')
ax.plot(bm13.index, bm13.Value, 'purple', label='v13.0.0')


plt.ylabel('%s (%s)' % (d[variable]['abbr'], d[variable]['unit']) )
plt.gca().xaxis.set_major_locator(mdates.MonthLocator(interval=1))
plt.legend()
plt.savefig('/users/mjr583/GC/eval_13.0.0/plots/allbenchmarks_%s.png' %(variable) )
plt.close()

f,ax= plt.subplots(figsize=(16,4))
ax.plot(df.index, (df.Value - df.Value.mean() ) / df.Value.mean(), 'k', label='CVAO')
ax.plot(df.index, (bm124.Value - bm124.Value.mean() ) / bm124.Value.mean(), '#fef0d9', label='v12.4.0')
ax.plot(df.index, (bm126.Value - bm126.Value.mean() ) / bm126.Value.mean(), '#fdcc8a', label='v12.6.0')
ax.plot(df.index, (bm128.Value - bm128.Value.mean() ) / bm128.Value.mean(), '#fc8d59', label='v12.8.0')
ax.plot(df.index, (bm1291.Value - bm1291.Value.mean() ) / bm1291.Value.mean(), 'green', label='v12.9.1')
ax.plot(df.index, (bm1293.Value - bm1293.Value.mean() ) / bm1293.Value.mean(), '#d7301f', label='v12.9.3')
ax.plot(bm13.index, (bm13.Value - bm13.Value.mean() ) / bm13.Value.mean(), 'purple', label='v13.0.0')


plt.ylabel(' Fractional $O_3$ seasonality')
plt.gca().xaxis.set_major_locator(mdates.MonthLocator(interval=1))
plt.legend()
plt.savefig('/users/mjr583/GC/eval_13.0.0/plots/frac_annual_diff_allbenchmarks_%s.png' %(variable) )
plt.close()
