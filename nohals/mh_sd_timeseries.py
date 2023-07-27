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

var, lat, lon, lev, time = GC.get_gc_var(rundir='tropchem_merra_4x5', variable=variable,year='2015')

nohal, lat,lon,lev,haltime = GC.get_gc_var(rundir='tropchem_merra_4x5', variable=variable, version='geos-chem',year='2015')


delta,interval=GC.find_timestep(time)
y = rp.find_nearest(lat, 53)
x = rp.find_nearest(lon, -9)
print(lon[x], lat[y])

var_time=[] ; nohal_time=[] ; var1=[]
for t in range(len(time)):
    v=var[t,0,y,x]
    var_time.append(v)

for t in range(len(haltime)):
    v=nohal[t,0,y,x]
    nohal_time.append(v)

gc=pd.DataFrame({'Value':var_time}, index=time)
nohal=pd.DataFrame({'Value':nohal_time}, index=haltime)

## Get Merge observations
gcmean=gc.resample('M').mean()
gc75=gc.resample('M').quantile(.75)
gc25=gc.resample('M').quantile(.25)

nomean=nohal.resample('M').mean()
nohal75=nohal.resample('M').quantile(.75)
nohal25=nohal.resample('M').quantile(.25)

df=pd.read_csv('/users/mjr583/scratch/NCAS_CVAO/GAW_datasets/mc_ozone_2006_to_2018.csv', delimiter=',', index_col=0)
df.index = pd.to_datetime(df.index)
df.columns = ['Value']

df=df[gc.index[0]:gc.index[-1]]
dfmean=df.resample('M').mean()
df75=df.resample('M').quantile(.75)
df25=df.resample('M').quantile(.25)

f,ax= plt.subplots(figsize=(11,4))
ax.plot(dfmean.index, dfmean.Value, 'k', label='Mace Head (Ire)')
ax.fill_between(dfmean.index, df25.Value,df75.Value, color='grey', alpha=.3)

ax.plot(gcmean.index, gcmean.Value, 'g', label='12.9.3')
ax.fill_between(gcmean.index, gc25.Value, gc75.Value, color='g', alpha=.3)

ax.plot(nomean.index, nomean.Value, 'purple', label='12.9.1_main.NOHALS')
ax.fill_between(nomean.index, nohal25.Value, nohal75.Value, color='purple', alpha=.3)

plt.ylabel('%s (%s)' % (d[variable]['abbr'], d[variable]['unit']) )
plt.gca().xaxis.set_major_locator(mdates.MonthLocator(interval=1))

plt.legend()
plt.savefig('plots/MH_timeseries_stddev_%s.png' %variable )
plt.close()
