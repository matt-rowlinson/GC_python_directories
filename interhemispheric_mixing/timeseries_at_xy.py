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

try:
    variable=sys.argv[1]
    X=float(sys.argv[2])
    Y=float(sys.argv[3])
    out=sys.argv[4]
except:
    variable='CO'
    X=-24.9
    Y=16.9
    out='CV'

var, lat, lon, lev, time = GC.get_gc_var(rundir='irma_025x03125', variable=variable)
nobb, lat,lon,lev,bbtime = GC.get_gc_var(rundir='irma_025x03125_noBB', variable=variable)
noAf, lat,lon,lev,aftime = GC.get_gc_var(rundir='irma_025x03125_noAfBB', variable=variable)

print(nobb.min(), nobb.mean(), nobb.max())
print(noAf.min(), noAf.mean(), noAf.max())

delta,interval=GC.find_timestep(time)
y = rp.find_nearest(lat, Y)
x = rp.find_nearest(lon, X)
print(lon[x], lat[y])

var_time=[] ; nobb_time=[] ; noAf_time=[]
for t in range(len(time)):
    v=var[t,0,y,x]
    var_time.append(v)

for t in range(len(bbtime)):
    v=nobb[t,0,y,x]
    nobb_time.append(v)

for t in range(len(aftime)):
    v=noAf[t,0,y,x]
    noAf_time.append(v)

gc=pd.DataFrame({'Value':var_time}, index=time)
nobb=pd.DataFrame({'Value':nobb_time}, index=bbtime)
noAf=pd.DataFrame({'Value':noAf_time}, index=aftime)

print(gc)
print(nobb)
print(noAf)

## Get Merge observations
#df=CV.get_from_merge(d[variable])
#df=df[gc.index[0]:gc.index[-1]]
#df=df.resample(delta).mean()
#df=df['2013']
f,ax= plt.subplots(figsize=(12,4))
#ax.plot(df.index, df.Value, 'k', label='CVAO')
ax.plot(gc.index, gc.Value, 'g', label='GEOS-Chem')
ax.plot(nobb.index, nobb.Value, 'orange',linestyle='--', label='No biomass burning')
ax.plot(noAf.index, noAf.Value, 'purple',linestyle='--', label='No African biomass burning')
plt.ylabel('%s (%s)' % (d[variable]['abbr'], d[variable]['unit']) )

plt.gca().xaxis.set_major_locator(mdates.DayLocator(interval=interval))

plt.legend()
plt.savefig('/users/mjr583/GC/interhemispheric_mixing/plots/timeseries_%s_%s.png' % (out, variable) )
