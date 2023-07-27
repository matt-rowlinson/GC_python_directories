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

variable='ethane'

for n in range(12):
    control, lat, lon, lev, time = GC.get_gc_var('tropchem_merra_4x5', variable=variable, version='12.9.3', year=2016)
    fltr=np.where(control < 0.)
    control[fltr] = np.nan

    nohals, lat, lon, lev, time = GC.get_gc_var('tropchem_merra_4x5', variable=variable,version='geos-chem', year='2016')
    fltr=np.where(nohals < 0.)
    nohals[fltr] = np.nan

    date=time[n].strftime(format='%m-%Y')
    control=control[n,0,:,:]
    nohals=nohals[n,0,:,:]

    Min, Max = rp.get_abs_max_min(control-nohals)

    fig,(ax1, ax2, ax3)= plt.subplots(1,3,figsize=(12,3))
    X,Y=np.meshgrid(lon,lat)
    m1=rp.get_basemap(ax=ax1, lines=False, freq=60.)#lllon=-70, urlon=0., lllat=-15., urlat=30.,ax=ax1)
    m2=rp.get_basemap(ax=ax2, lines=False, freq=60.)#lllon=-70, urlon=0., lllat=-15., urlat=30.,ax=ax2)
    m3=rp.get_basemap(ax=ax3, lines=False, freq=60.)#lllon=-70, urlon=0., lllat=-15., urlat=30.,ax=ax3)

    MAX=np.nanmax(control)
    diff = - (100 - (control/nohals * 100))
    Min, Max = rp.get_abs_max_min(diff)

    if variable == 'NO' or variable == 'NO2':

        im=ax1.pcolormesh(X,Y, control,  cmap=cmap, vmax=MAX, norm=LogNorm() )
        im2=ax2.pcolormesh(X,Y, nohals,  cmap=cmap, vmax=MAX, norm=LogNorm())
        im3=ax3.pcolormesh(X,Y, diff,  cmap='bwr', vmin=Min, vmax=Max)

    else:
        im=ax1.pcolormesh(X,Y, control,  cmap=cmap, vmax=MAX)
        im2=ax2.pcolormesh(X,Y, nohals,  cmap=cmap, vmax=MAX)
        im3=ax3.pcolormesh(X,Y, diff,  cmap='bwr', vmin=Min, vmax=Max)

    ax1.set_title('Control 12.9.3')
    ax2.set_title('12.9.1_main.NOHALS')
    ax3.set_title('Control - NOHALS')

    plt.subplots_adjust(bottom=0.24)
    cbar_ax = fig.add_axes([0.126, 0.15, 0.502, 0.07])
    fig.colorbar(im, cax=cbar_ax, orientation='horizontal')
    cbar_ax.set_xlabel('%s (%s)' %(variable, d[variable]['unit']))

    cbar_ax = fig.add_axes([0.673,0.15,0.225, 0.07])
    fig.colorbar(im3, cax=cbar_ax, orientation='horizontal')

    cbar_ax.set_xlabel(r'$\Delta$$O_3$ (%)')

    plt.suptitle('SpeciesConc_%s  %s' %(variable, date))

    plt.savefig('plots/SpeciesConc_map_pc_%s_%s.png' %(variable, date)) 
    plt.close()
