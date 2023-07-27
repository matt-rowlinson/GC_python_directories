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
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import matplotlib
matplotlib.use('agg')
sys.path.append('/users/mjr583/python_lib')
import GC_tools as GC
import RowPy as rp
import matplotlib.colors as mcolors
from CVAO_dict import CVAO_dict as d

current_dir = os.path.dirname(__file__)
rgb_WhGrYlRd = np.genfromtxt('/users/mjr583/python_lib/colormaps/WhGrYlRd.txt',
                                     delimiter=' ')
WhGrYlRd = mcolors.ListedColormap(rgb_WhGrYlRd/255.0)
cmap=WhGrYlRd

variable='CO'
out='newcontrol'

control_bc, lat, lon, lev, time = GC.get_gc_bc('irma_025x03125', var=variable)
fltr=np.where(control_bc < 0.)
control_bc[fltr] = np.nan

nobb_bc, lat, lon, lev, notime = GC.get_gc_bc('irma_025x03125_noBB', variable)
fltr=np.where(nobb_bc < 0.)
nobb_bc[fltr] = np.nan

noaf_bc, lat, lon, lev, aftime = GC.get_gc_bc('irma_025x03125_noAfBB', variable)
fltr=np.where(noaf_bc < 0.)
noaf_bc[fltr] = np.nan

control=np.mean(control_bc[124:126,0,:,:],0)
nobb=np.mean(nobb_bc[124:126,0,:,:],0)
noaf=np.mean(noaf_bc[124:126,0,:,:],0)

print(control.shape, control.mean())
print(nobb.shape, nobb.mean())
print(noaf.shape, noaf.mean())

fig,(ax1, ax2, ax3)= plt.subplots(1,3,figsize=(12,3))
X,Y=np.meshgrid(lon,lat)
m=rp.get_basemap(ax=ax1)
m=rp.get_basemap(ax=ax2)
m=rp.get_basemap(ax=ax3)

Max=np.nanmax(control)
im=ax1.pcolormesh(X,Y, control,  cmap=cmap, vmax=Max, norm=LogNorm())#, vmin=vmin,norm=norm)
im=ax2.pcolormesh(X,Y, noaf,  cmap=cmap, vmax=Max, norm=LogNorm())#, vmin=-4000., vmax=4000)#vmin,norm=norm)
im=ax3.pcolormesh(X,Y, nobb,  cmap=cmap, vmax=Max, norm=LogNorm())#vmin=-100., vmax=100.)#, vmin=vmin,norm=norm)

ax1.set_title('Control GC_BC')
ax2.set_title('No African biomass burning')
ax3.set_title('No biomass burning')

#3cbar_ax = fig.add_axes([0.15, 0.1, 0.75, 0.07])
#fig.colorbar(im, cax=cbar_ax, orientation='horizontal')

plt.suptitle('GC Boundary Conditions')

plt.savefig('/users/mjr583/GC/interhemispheric_mixing/plots/mapGC_BC_%s_%s.png' %(variable, out)) 
plt.close()
