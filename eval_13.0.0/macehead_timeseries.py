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
var, lat, lon, lev, time = GC.get_gc_var(rundir, variable, version, year='2016')
delta,interval=GC.find_timestep(time)
y = rp.find_nearest(lat, 53.33)
x = rp.find_nearest(lon, -9.9)
var_time=[]
for t in range(len(time)):
    v=var[t,0,y,x]
    var_time.append(v)
merra=pd.DataFrame({'Value':var_time}, index=time)
merra=merra.resample('D').mean()

print('get v13')
rundir='fullchem_hourly_default'
var, lat, lon, lev, time = GC.get_gc_var(rundir, variable, version='GEOS-Chem', year='2016')
delta,interval=GC.find_timestep(time)
var_time=[]
for t in range(len(time)):
    v=var[t,0,y,x]
    var_time.append(v)
bm=pd.DataFrame({'Value':var_time}, index=time)
bm=bm.resample('D').mean()


## Get Merge observations
df=pd.read_csv('/users/mjr583/scratch/NCAS_CVAO/GAW_datasets/mc_ozone_2006_to_2018.csv', delimiter=',', index_col=0)
df.index = pd.to_datetime(df.index)
df.columns = ['Value']
df=df[merra.index[0]:merra.index[-1]]
df=df.resample('D').mean()

#new_index=[]
#for dt in df.index:
#    dtt=dt.replace(day=1)
#    new_index.append(dtt)
#df.index=new_index

f,ax= plt.subplots(figsize=(16,4))
ax.plot(df.index, df.Value, 'k', label='Mace Head Observations')
ax.plot(merra.index, merra.Value, 'g', label='v12.9.3')
ax.plot(bm.index, bm.Value, 'r', label='v13.0.0')

plt.ylabel('%s (%s)' % (d[variable]['abbr'], d[variable]['unit']) )

plt.gca().xaxis.set_major_locator(mdates.MonthLocator(interval=1))

#plt.ylim((7.5,52.5))
plt.legend()
plt.savefig('/users/mjr583/GC/eval_13.0.0/plots/macehead_2016_%s.png' %(variable) )
plt.close()

