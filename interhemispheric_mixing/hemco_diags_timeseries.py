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
variable='CO'
version='12.9.3'

df,ems = GC.plot_all_gc_ems('irma_025x03125', 'CO')

print(df)

nobb, bbems = GC.plot_all_gc_ems('irma_025x03125_noAfrBB', 'CO')

print(nobb)

sys.exit()

delta,interval=GC.find_timestep(time)
y = rp.find_nearest(lat, 16.9)
x = rp.find_nearest(lon, -24.9)
print(lon[x], lat[y])

var_time=[] ; nobb_time=[]
for t in range(len(time)):
    v=var[t,0,y,x]
    var_time.append(v)

    if t<len(nobb):
        bb=nobb[t,0,y,x]
        nobb_time.append(bb)

gc=pd.DataFrame({'Value':var_time}, index=time)
nobb=pd.DataFrame({'Value':nobb_time}, index=bbtime)

plt.plot(gc.Value)
plt.plot(nobb.Value)
plt.savefig('/users/mjr583/scratch/GC/%s/rundirs/%s/plots/temp.png' % (version, rundir) )
plt.close()

## Get Merge observations
df=CV.get_from_merge(d[variable])
df=df[gc.index[0]:gc.index[-1]]
df=df.resample(delta).mean()
#df=df['2013']
f,ax= plt.subplots(figsize=(12,4))
ax.plot(df.index, df.Value, 'k', label='Congo')
ax.plot(gc.index, gc.Value, 'g', label='GEOS-Chem')
ax.plot(nobb.index, nobb.Value, 'orange',linestyle='--', label='No BB')
plt.ylabel('%s (%s)' % (d[variable]['abbr'], d[variable]['unit']) )

plt.gca().xaxis.set_major_locator(mdates.DayLocator(interval=interval))

plt.legend()
plt.savefig('/users/mjr583/scratch/GC/%s/rundirs/%s/plots/timeseries_%s.png' % (version, rundir, variable) )
