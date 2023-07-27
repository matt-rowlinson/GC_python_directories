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
out='COconc'

control_bc, lat, lon, lev, time = GC.get_gc_var('irma_025x03125', 'CO')
fltr=np.where(control_bc < 0.)
control_bc[fltr] = np.nan

nobb_bc, lat, lon, lev, time = GC.get_gc_var('irma_025x03125_noBB', 'CO')
fltr=np.where(nobb_bc < 0.)
nobb_bc[fltr] = np.nan

noaf_bc, lat, lon, lev, time = GC.get_gc_var('irma_025x03125_noAfBB', 'CO')
fltr=np.where(noaf_bc < 0.)
noaf_bc[fltr] = np.nan

date='24-08-17_0000z'
n=168
control=control_bc[n,0,:,:]
nobb=nobb_bc[n,0,:,:]
noaf=noaf_bc[n,0,:,:]

print(control.shape)
print(nobb.shape)
print(noaf.shape)

fig,(ax1, ax2, ax3)= plt.subplots(1,3,figsize=(12,3))
X,Y=np.meshgrid(lon,lat)
m1=rp.get_basemap(lllon=-70, urlon=0., lllat=-15., urlat=30.,ax=ax1)
m2=rp.get_basemap(lllon=-70, urlon=0., lllat=-15., urlat=30.,ax=ax2)
m3=rp.get_basemap(lllon=-70, urlon=0., lllat=-15., urlat=30.,ax=ax3)

Max=np.nanmax(control)
im=ax1.pcolormesh(X,Y, control,  cmap=cmap, vmax=Max, norm=LogNorm())#, vmin=vmin,norm=norm)
im=ax2.pcolormesh(X,Y, noaf,  cmap=cmap, vmax=Max, norm=LogNorm())#, vmin=-4000., vmax=4000)#vmin,norm=norm)
im=ax3.pcolormesh(X,Y, nobb,  cmap=cmap, vmax=Max, norm=LogNorm())#vmin=-100., vmax=100.)#, vmin=vmin,norm=norm)

ax1.set_title('Control Surface %s' %variable)
ax2.set_title('No African biomass burning')
ax3.set_title('No biomass burning')

#cbar_ax = fig.add_axes([0.1, 0.1, 0.65, 0.07])
#fig.colorbar(im, cax=cbar_ax, orientation='horizontal')

plt.suptitle('SpeciesConc_%s  %s' %(variable, date))

plt.savefig('/users/mjr583/GC/interhemispheric_mixing/plots/SpeciesConc_map_%s_%s.png' %(variable, date)) 
plt.close()
