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
year='2014'
rundir='tropchem_merra_4x5'
var, lat, lon, lev, time = GC.get_gc_var(rundir, variable, version, year=year)
delta,interval=GC.find_timestep(time)
y = rp.find_nearest(lat, 13.)
x = rp.find_nearest(lon, -59.)
var_time=[]
for t in range(len(time)):
    v=var[t,0,y,x]
    var_time.append(v)
merra=pd.DataFrame({'Value':var_time}, index=time)
if delta=='M':
    merra=merra.resample(delta).mean()

## Get Merge observations
df=pd.read_csv('/users/mjr583/scratch/NCAS_CVAO/GAW_datasets/rag_ozone_2007_to_2016.csv', delimiter=',', index_col=0)
df.index = pd.to_datetime(df.index)
df.columns = ['Value']
df=df[merra.index[0]:merra.index[-1]]
df=df.resample(delta).mean()

new_index=[]
for dt in df.index:
    dtt=dt.replace(day=1)
    new_index.append(dtt)
df.index=new_index

f,ax= plt.subplots(figsize=(16,4))
ax.plot(df.index, df.Value, 'k', label='Ragged Point (Bar) Observations')
ax.plot(merra.index, merra.Value, 'g', label='GC Merra')
#ax.plot(geosfp.index, geosfp.Value, 'lightgreen', label='GC Geos-FP')
#ax.plot(bm.index, bm.Value, 'r', label='GC Benchmark')

plt.ylabel('%s (%s)' % (d[variable]['abbr'], d[variable]['unit']) )

plt.gca().xaxis.set_major_locator(mdates.MonthLocator(interval=1))

#plt.ylim((7.5,52.5))
plt.legend()
plt.savefig('/users/mjr583/GC/eval_12.9.3/plots/ragpoint_%s_%s.png' %(year, variable) )
plt.close()

