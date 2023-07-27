#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#SBATCH --job-name=pcolor_plot
#SBATCH --ntasks=1
#SBATCH --mem=2gb
#SBATCH --partition=nodes
#SBATCH --time=00:01:00
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
from datetime import datetime

rundir=sys.argv[1]
species=sys.argv[2]
dt=sys.argv[3]
dt=datetime.strptime(datetime.strptime(dt,'%Y-%m-%d-%H').strftime('%Y-%m-%d %H:%M'), '%Y-%m-%d %H:%M')
version='12.9.3'

var, lat, lon, lev, time = GC.get_gc_var(rundir, species, version)
jobid=rp.find_nearest(time, dt)

f,ax= plt.subplots(figsize=(8,8))
X,Y=np.meshgrid(lon,lat)
m=rp.get_basemap(resolution='i', lines=True, lllat=lat.min(), lllon=lon.min(), urlat=lat.max(), urlon=lon.max(), ax=ax)
im=ax.pcolormesh(X,Y, var[jobid, 0, :,:], vmin=0., vmax=var[:,0, 3:-3, 3:-3].max(), cmap='viridis')
cbar = f.colorbar(im,orientation='horizontal')
cbar.ax.set_xlabel('%s (%s)' % (d[species]['abbr'], d[species]['unit']))
plt.title(time[jobid], fontsize=14)
plt.savefig('/users/mjr583/scratch/GC/%s/%s/plots/%s_%s.png' % (version, rundir, species, str(time[jobid]).replace(" ","_")) )
