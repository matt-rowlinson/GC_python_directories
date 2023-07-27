#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#SBATCH --job-name=pcolor_plot
#SBATCH --ntasks=1
#SBATCH --mem=25Gb
#SBATCH --partition=nodes
#SBATCH --time=00:11:30
#SBATCH --output=LOGS/pcolormesh_%a.log
#SBATCH --array=1-2
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('agg')
sys.path.append('/users/mjr583/python_lib')
import GC_tools as GC
import RowPy as rp
from CVAO_dict import CVAO_dict as d

jobid=int(os.getenv("SLURM_ARRAY_TASK_ID"))
inputs=GC.get_arguments()
rundir=inputs.rundir
variable=inputs.var
version=inputs.version

var, lat, lon, lev, time = GC.get_gc_var(rundir, variable, version)
if variable=='CO':
    Max=180.
elif variable=='O3':
    Max=40
f,ax= plt.subplots(figsize=(8,8))
X,Y=np.meshgrid(lon,lat)
m=rp.get_basemap(resolution='i', lines=True, lllat=lat.min(), lllon=lon.min(), urlat=lat.max(), urlon=lon.max(), ax=ax)
im=ax.pcolormesh(X,Y, var[jobid, 0, :,:], vmin=0., vmax=Max ,cmap='viridis')#var[:,0, 60:-6, 60:-60].max(), cmap='viridis')
cbar = f.colorbar(im,orientation='horizontal')
cbar.ax.set_xlabel('%s (%s)' % (d[variable]['abbr'], d[variable]['unit']))
plt.title(time[jobid], fontsize=14)
plt.savefig('/users/mjr583/scratch/GC/%s/rundirs/%s/plots/pcolorm_%s_%s.png' % (version, rundir, variable, time[jobid]) )
