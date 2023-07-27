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

variable='ethane'
version='12.9.3'
rundir='tropchem_merra_4x5'
var, lat, lon, lev, time = GC.get_gc_var(rundir, variable, version, year='2016')
bench, lat, lon, lev, benchtime = GC.get_gc_var(rundir='benchmark', variable=variable, version=version, year='2016')

delta,interval=GC.find_timestep(time)
y = rp.find_nearest(lat, 16.9)
x = rp.find_nearest(lon, -24.9)

var_time=[] ; bench_time=[]
for t in range(len(time)):
    v=var[t,0,y,x]
    var_time.append(v)

    b=bench[t,0,y,x]
    bench_time.append(b)
gc=pd.DataFrame({'Value':var_time}, index=time)
bc=pd.DataFrame({'Value':bench_time}, index=time)

print(gc)
print(bc)

## Get Merge observations
df=CV.get_from_merge(d[variable])
df=df[gc.index[0]:gc.index[-1]]
df=df.resample(delta).mean()
print(df)
f,ax= plt.subplots(figsize=(16,4))
ax.plot(df.index, df.Value, 'k', label='CVAO')
ax.plot(df.index, gc.Value, 'g', label='GEOS-Chem')
ax.plot(df.index, bc.Value, 'r', label='Benchmark')

plt.ylabel('%s (%s)' % (d[variable]['abbr'], d[variable]['unit']) )

plt.gca().xaxis.set_major_locator(mdates.MonthLocator(interval=interval))

plt.ylim((7.5,52.5))
plt.legend()
plt.savefig('/users/mjr583/GC/eval_12.9.3/plots/benchmark_%s_timeseries_%s.png' %(rundir,variable) )
plt.close()

gcmean=gc.resample('M').mean()
gc75=gc.resample('M').quantile(.75)
gc25=gc.resample('M').quantile(.25)

bcmean=bc.resample('M').mean()
bc75=bc.resample('M').quantile(.75)
bc25=bc.resample('M').quantile(.25)

df=CV.get_from_merge(d[variable])
df=df[gc.index[0]:gc.index[-1]]
dfmean=df.resample('M').mean()
df75=df.resample('M').quantile(.75)
df25=df.resample('M').quantile(.25)

f,ax= plt.subplots(figsize=(6,4))
ax.plot(dfmean.index, dfmean.Value, 'k', label='CVAO')
ax.fill_between(dfmean.index, df25.Value,df75.Value, color='grey', alpha=.3)

ax.plot(gcmean.index, gcmean.Value, 'g', label='12.9.3 (my GC)')
ax.fill_between(gcmean.index, gc25.Value, gc75.Value, color='g', alpha=.3)

ax.plot(bcmean.index, bcmean.Value, 'r', label='Benchmark')
ax.fill_between(bcmean.index, bc25.Value, bc75.Value, color='r', alpha=.3)

plt.ylabel('%s (%s)' % (d[variable]['abbr'], d[variable]['unit']) )
plt.gca().xaxis.set_major_locator(mdates.MonthLocator(interval=interval))

#plt.ylim((7.5,52.5))
plt.legend()
plt.savefig('/users/mjr583/GC/eval_12.9.3/plots/benchmark_comparison_%s.png' %variable )
plt.close()
