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
import numpy as np
matplotlib.use('agg')
sys.path.append('/users/mjr583/python_lib')
import GC_tools as GC
import RowPy as rp
from CVAO_dict import CVAO_dict as d
import CVAO_tools as CV
plt.style.use('seaborn-darkgrid')

inputs=GC.get_arguments()
rundir='irma_025x03125_noBB'
variable='EmisCO_BioBurn'
version='12.9.3'

var, time, lat, lon, lev, area = GC.HEMCO_Diag_read(rundir='irma_025x03125', variable=variable)

nobb, bbtime, lat,lon,lev, area = GC.HEMCO_Diag_read(rundir='irma_025x03125_noBB', variable=variable)

noAf, aftime, lat,lon,lev,area = GC.HEMCO_Diag_read(rundir='irma_025x03125_noAfBB', variable=variable)

var=np.sum(np.sum(np.sum(var,-1),-1),-1)
nobb=np.sum(np.sum(np.sum(nobb,-1),-1),-1)

noAf=np.sum(np.sum(np.sum(noAf,-1),-1),-1)
print(var.shape)


delta,interval=GC.find_timestep(time)
y = rp.find_nearest(lat, 16.9)
x = rp.find_nearest(lon, -24.9)
print(lon[x], lat[y])


gc=pd.DataFrame({'Value':var}, index=time)
nobb=pd.DataFrame({'Value':nobb}, index=bbtime)
noAf=pd.DataFrame({'Value':noAf}, index=aftime)

print(gc)
print(nobb)
print(noAf)

#gc=gc[50:]
#nobb=nobb[50:]

plt.plot(gc.Value)
plt.plot(nobb.Value)
plt.savefig('/users/mjr583/scratch/GC/%s/rundirs/%s/plots/temp.png' % (version, rundir) )
plt.close()

variable=variable[4:6]
#df=df['2013']
f,ax= plt.subplots(figsize=(12,4))
ax.plot(gc.index, gc.Value, 'g', label='GEOS-Chem')
ax.plot(nobb.index, nobb.Value, 'orange',linestyle='--', label='No BB')
ax.plot(noAf.index, noAf.Value, 'purple',linestyle='--', label='No Af BB')
plt.ylabel('%s (%s)' % (d[variable]['abbr'], d[variable]['unit']) )

plt.gca().xaxis.set_major_locator(mdates.DayLocator(interval=interval))

plt.legend()
plt.savefig('/users/mjr583/GC/interhemispheric_mixing/plots/hemco_timeseries_%s.png' %variable )
