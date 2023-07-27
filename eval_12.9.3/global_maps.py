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

inputs=GC.get_arguments()
rundir=inputs.rundir
variable=inputs.var
version=inputs.version

variables=['O3','CO','NO2','NO','ethane','propane']
logscale=False
europeswitch=False
suff=''
if europeswitch:
    latmin=30. 
    latmax=75.
    lonmin=-10.
    lonmax=50.
    suff+='_europe'
else:
    var_array, lat, lon, lev, time = GC.get_gc_var(rundir, variable, version, year='2016')
    latmin=lat.min() 
    latmax=lat.max() 
    lonmin=lon.min() 
    lonmax=lon.max() 

if logscale:
    norm=matplotlib.colors.LogNorm()
    vmin=None
    suff+='_log'
else:
    norm=None
    vmin=0.
for variable in variables:
    print(variable)
    var_array, lat, lon, lev, time = GC.get_gc_var(rundir, variable, version, year='2016')
    fltr=np.where(var_array < 0.)
    var_array[fltr] = np.nan

    var=np.mean(var_array[:,0,:,:],0)

    f,ax= plt.subplots(figsize=(8,8))
    X,Y=np.meshgrid(lon,lat)
    m=rp.get_basemap(resolution='i', lines=False, lllat=latmin, lllon=lonmin, urlat=latmax, urlon=lonmax, ax=ax)
    im=ax.pcolormesh(X,Y, var,  cmap=cmap, vmin=vmin,norm=norm)
    cbar = f.colorbar(im,orientation='horizontal')
    cbar.ax.set_xlabel('%s (%s)' % (d[variable]['abbr'], d[variable]['unit']))
    #plt.title(time[jobid], fontsize=14)
    plt.savefig('/users/mjr583/GC/eval_12.9.3/plots/global_map_annualmean_%s%s.png' % (variable, suff ))
    plt.close()
    '''
    var=np.mean( np.concatenate([var_array[-31:,0,:,:],var_array[:60,0,:,:]]),0)
    f,ax= plt.subplots(figsize=(8,8))
    m=rp.get_basemap(resolution='i', lines=False, lllat=latmin, lllon=lonmin, urlat=latmax, urlon=lonmax, ax=ax)
    im=ax.pcolormesh(X,Y, var, vmin=vmin, cmap=cmap, norm=norm)
    cbar = f.colorbar(im,orientation='horizontal')
    cbar.ax.set_xlabel('%s (%s)' % (d[variable]['abbr'], d[variable]['unit']))
    #plt.title(time[jobid], fontsize=14)
    plt.savefig('/users/mjr583/scratch/GC/%s/%s/plots/global_map_DJFmean_%s%s.png' % (version, rundir, variable, suff))
    plt.close()

    var=np.mean( var_array[152:244,0,:,:],0 )
    f,ax= plt.subplots(figsize=(8,8))
    m=rp.get_basemap(resolution='i', lines=False, lllat=latmin, lllon=lonmin, urlat=latmax, urlon=lonmax, ax=ax)
    im=ax.pcolormesh(X,Y, var, vmin=vmin, cmap=cmap, norm=norm)
    cbar = f.colorbar(im,orientation='horizontal')
    cbar.ax.set_xlabel('%s (%s)' % (d[variable]['abbr'], d[variable]['unit']))
    #plt.title(time[jobid], fontsize=14)
    plt.savefig('/users/mjr583/scratch/GC/%s/%s/plots/global_map_JJAmean_%s%s.png' % (version, rundir, variable, suff ))
    plt.close()

    var=np.mean( var_array[91:122,0,:,:],0 )
    f,ax= plt.subplots(figsize=(8,8))
    m=rp.get_basemap(resolution='i', lines=False, lllat=latmin, lllon=lonmin, urlat=latmax, urlon=lonmax, ax=ax)
    im=ax.pcolormesh(X,Y, var, vmin=vmin, cmap=cmap, norm=norm)
    cbar = f.colorbar(im,orientation='horizontal')
    cbar.ax.set_xlabel('%s (%s)' % (d[variable]['abbr'], d[variable]['unit']))
    #plt.title(time[jobid], fontsize=14)
    plt.savefig('/users/mjr583/scratch/GC/%s/%s/plots/global_map_MAMmean_%s%s.png' % (version, rundir, variable, suff ))
    plt.close()

    var=np.mean( var_array[244:335,0,:,:],0 )
    f,ax= plt.subplots(figsize=(8,8))
    m=rp.get_basemap(resolution='i', lines=False, lllat=latmin, lllon=lonmin, urlat=latmax, urlon=lonmax, ax=ax)
    im=ax.pcolormesh(X,Y, var, vmin=vmin, cmap=cmap, norm=norm)
    cbar = f.colorbar(im,orientation='horizontal')
    cbar.ax.set_xlabel('%s (%s)' % (d[variable]['abbr'], d[variable]['unit']))
    #plt.title(time[jobid], fontsize=14)
    plt.savefig('/users/mjr583/scratch/GC/%s/%s/plots/global_map_SONmean_%s%s.png' % (version, rundir, variable, suff ))
    plt.close()
    '''
