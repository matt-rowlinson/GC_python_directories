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

var, lat, lon, lev, time = GC.get_gc_var(rundir='tropchem_merra_4x5', variable=variable,year='2016')

nohal, lat,lon,lev,haltime = GC.get_gc_var(rundir='tropchem_merra_4x5', variable=variable, version='geos-chem',year='2016')

thir, lat,lon,lev,ttime = GC.get_gc_var(rundir='fullchem_4x5_LVOCfalse', variable=variable, version='GEOS-Chem',year='2016')



delta,interval=GC.find_timestep(time)
y = rp.find_nearest(lat, 16.9)
x = rp.find_nearest(lon, -24.9)
print(lon[x], lat[y])

var_time=[] ; nohal_time=[] ; thirteen=[]
for t in range(len(time)):
    v=var[t,0,y,x]
    var_time.append(v)

for t in range(len(haltime)):
    v=nohal[t,0,y,x]
    nohal_time.append(v)

for t in range(len(ttime)):
    v=thir[t,0,y,x]
    thirteen.append(v)

gc=pd.DataFrame({'Value':var_time}, index=time)
nohal=pd.DataFrame({'Value':nohal_time}, index=haltime)
thirt=pd.DataFrame({'Value':thirteen}, index=ttime)


## Get Merge observations
gcmean=gc.resample('M').mean()
gc75=gc.resample('M').quantile(.75)
gc25=gc.resample('M').quantile(.25)

nomean=nohal.resample('M').mean()
nohal75=nohal.resample('M').quantile(.75)
nohal25=nohal.resample('M').quantile(.25)

Thmean=thirt.resample('M').mean()
Th75=thirt.resample('M').quantile(.75)
Th25=thirt.resample('M').quantile(.25)


df=CV.get_from_merge(d[variable])
df=df[gc.index[0]:gc.index[-1]]
dfmean=df.resample('M').mean()
df75=df.resample('M').quantile(.75)
df25=df.resample('M').quantile(.25)

f,ax= plt.subplots(figsize=(11,4))
ax.plot(dfmean.index, dfmean.Value, 'k', label='CVAO')
ax.fill_between(dfmean.index, df25.Value,df75.Value, color='grey', alpha=.3)

ax.plot(gcmean.index, gcmean.Value, 'g', label='v12.9.3')
ax.fill_between(gcmean.index, gc25.Value, gc75.Value, color='g', alpha=.3)

ax.plot(nomean.index, nomean.Value, 'purple', label='v12.9.1.NOHALS')
ax.fill_between(nomean.index, nohal25.Value, nohal75.Value, color='purple', alpha=.3)

#ax.plot(Thmean.index, Thmean.Value, 'r', label='v13.0.0')
#ax.fill_between(Thmean.index, Th25.Value, Th75.Value, color='r', alpha=.3)

plt.ylabel('%s (%s)' % (d[variable]['abbr'], d[variable]['unit']) )
plt.gca().xaxis.set_major_locator(mdates.MonthLocator(interval=1))

plt.legend()
plt.savefig('plots/timeseries_stddev_%s.png' %variable )
plt.close()
