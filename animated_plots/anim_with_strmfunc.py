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
ms=GC.closest_met(time[jobid])

t=GC.hour_rounder(time[jobid])
diff=GC.dt_difference(t, ms)

if diff >= 1.:
    print('Need to interpolate winds')
    U , mlat,mlon=GC.get_interp_gc_input(ms, t, var='U', filetype='A3dyn')
    V , mlat,mlon=GC.get_interp_gc_input(ms, t, var='V', filetype='A3dyn')

else:
    U , mlat,mlon=GC.get_gc_input(ms, var='U', filetype='A3dyn')
    V , mlat,mlon=GC.get_gc_input(ms, var='V', filetype='A3dyn')
mlat[0]=-90.
mlat[-1]=90.
U=U[0]
V=V[0]
speed = np.sqrt(U**2 + V**2)
lw = 8*speed / speed.max()
print(speed.max(), speed.min())
print(lw.max(), speed.min())

f,ax= plt.subplots(figsize=(8,8))
X,Y=np.meshgrid(lon,lat)
m=rp.get_basemap(resolution='i', lines=True, lllat=lat.min(), lllon=lon.min(), urlat=lat.max(), urlon=lon.max(), ax=ax)

if variable=='O3':
    vmax=50.
elif variable=='CO':
    vmax=200.
elif variable=='ethane':
    vmax=5000.
elif variable=='propane':
    vmax=200.

im=ax.pcolormesh(X,Y, var[jobid, 0, :,:], vmin=0., vmax=var[:,0, 60:-6, 60:-60].max(), cmap='viridis')
strm = ax.streamplot(mlon,mlat,U,V,density=4, zorder=11, color='w', linewidth=lw)#, levels=levels,zorder=11)
#GC.makeStreamLegend(strm, ax, GC.LWToSpeed, nlines=5, color='k', fmt='{:6.1f}')

cbar = f.colorbar(im,orientation='horizontal')
cbar.ax.set_xlabel('%s (%s)' % (d[variable]['abbr'], d[variable]['unit']))
plt.title(time[jobid], fontsize=14)
plt.savefig('/users/mjr583/scratch/GC/%s/rundirs/%s/plots/pcolorm_strmfunc_%s_%s.png' % (version, rundir, variable, time[jobid]) )
