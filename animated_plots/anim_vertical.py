#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#SBATCH --job-name=pcolor_plot
#SBATCH --ntasks=1
#SBATCH --mem=15gb
#SBATCH --partition=nodes
#SBATCH --time=00:11:30
#SBATCH --output=LOGS/pcolormesh_%a.log
#SBATCH --array=1-2
import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('agg')
sys.path.append('/users/mjr583/python_lib')
import GC_tools as GC
import RowPy as rp
from CVAO_dict import CVAO_dict as d

jobid=int(os.getenv("SLURM_ARRAY_TASK_ID"))
#jobid=3
inputs=GC.get_arguments()
rundir=inputs.rundir
variable=inputs.var
version=inputs.version

var, lat, lon, lev, time = GC.get_gc_var(rundir, variable, version)
lev=pd.read_csv('./GC_vertical_levels.csv')['Pressure (hPa)'][:32]
#lev_72=pd.read_csv('./GC_72_vertical_levels.csv')['Pressure (hPa)'][:]

lat=lat[3:-3]
lodx=rp.find_nearest(lon,-24.9)
var = np.mean(var[:,:32,3:-3,lodx-10:lodx+10],3)
ms=GC.closest_met(time[jobid])

t=GC.hour_rounder(time[jobid])
diff=GC.dt_difference(t, ms)

if diff >= 1.:
    print('Need to interpolate winds')
    OMEGA, mlat,mlon=GC.get_interp_gc_input(ms, t, var='OMEGA', filetype='A3dyn')
else:
    OMEGA, mlat,mlon=GC.get_gc_input(ms, var='OMEGA', filetype='A3dyn')
mlat[0]=-90.
mlat[-1]=90.
aa = rp.find_nearest(mlat, lat[0])
bb = rp.find_nearest(mlat, lat[-1])
mlat=mlat[aa:bb+1]

aidx=rp.find_nearest(mlon, lon[lodx-10])
bidx=rp.find_nearest(mlon, lon[lodx+10])
OMEGA=np.mean(OMEGA[:32,aa:bb+1,aidx:bidx],2)*1e3

f,ax= plt.subplots(figsize=(8,8))
X,Y=np.meshgrid(lat,lev)
mX, mY = np.meshgrid(mlat, lev)

Min=0.
if variable=='O3' or variable =='NO2' or variable=='propane':
    Max=60.
elif variable=='CO':
    Max=150.
elif variable=='ethane':
    Max=800.
    Min=300.
elif variable=='meoh':
    Max=600.
elif variable=='CH4':
    Min=1820.
    Max=1900.
elif variable=='acetone':
    Min=200.
    Max=800.

im=ax.pcolormesh(X,Y, var[jobid, :,:], vmin=Min, vmax=Max, cmap='viridis')
ladx=rp.find_nearest(lat,16.9)

plt.scatter(lat[ladx],lev[0]-6, marker='x', color='r',zorder=11)

V=ax.contour(mX, mY, OMEGA,colors='w')#, 'w')
ax.clabel(V, inline=1, fmt='%1d', fontsize=10)

plt.gca().invert_yaxis()
plt.ylabel('Pressure (hPa)')
plt.xlabel('Latitude')
cbar = f.colorbar(im,orientation='horizontal')
cbar.ax.set_xlabel('%s (%s)' % (d[variable]['abbr'], d[variable]['unit']))
plt.title(time[jobid], fontsize=14)
plt.savefig('/users/mjr583/scratch/GC/%s/rundirs/%s/plots/vertical_%s_%s.png' % (version, rundir, variable, time[jobid]) )
