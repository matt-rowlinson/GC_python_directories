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
from matplotlib.colors import BoundaryNorm
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
cmap=matplotlib.cm.bwr

inputs=GC.get_arguments()
rundir=inputs.rundir
variable=inputs.var
version='12.9.3'

variables=['O3','CO']#,'NO2','NO','ethane','propane']
logscale=False
europeswitch=False
suff=''
t, lat, lon, lev, time = GC.get_gc_var(rundir, 'O3', version, year='2016')
if europeswitch:
    latmin=30. 
    latmax=75.
    lonmin=-10.
    lonmax=50.
    suff+='_europe'
else:
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
    vmin=None


for variable in variables:
    print(variable)
    Max=10.
    benchmark, lat, lon, lev, time = GC.get_gc_var(rundir='benchmark', variable=variable, version='12.9.3')
    bench=np.mean(benchmark[:,0,:,:],0)

    var_array, lat, lon, lev, time = GC.get_gc_var(rundir, variable, version, year='2016')
    fltr=np.where(var_array < 0.)
    var_array[fltr] = np.nan
    var=np.mean(var_array[:,0,:,:],0)
    diff=100 - (bench/var*100)
    bounds = np.arange(-10,10.5,.5)
    norm = BoundaryNorm(bounds, cmap.N)

    f,ax= plt.subplots(figsize=(8,8))
    X,Y=np.meshgrid(lon,lat)
    m=rp.get_basemap(resolution='i', lines=False, lllat=latmin, lllon=lonmin, urlat=latmax, urlon=lonmax, ax=ax)
    im=ax.pcolormesh(X,Y, diff,  cmap=cmap, vmin=vmin,norm=norm)
    cbar = f.colorbar(im,orientation='horizontal')
    cbar.ax.set_xlabel('%s (%s)' % (d[variable]['abbr'], d[variable]['unit']))
    plt.title('Annual mean: Benchmark - 12.9.3', fontsize=14)
    plt.savefig('/users/mjr583/GC/eval_12.9.3/plots/benchmark-diff_annual_%s%s.png' % (variable, suff ))
    plt.close()

    var=np.mean( np.concatenate([var_array[-1:,0,:,:],var_array[:2,0,:,:]]),0)
    bench=np.mean( np.concatenate( [ benchmark[-1:,0,:,:], benchmark[:2,0,:,:] ] ),0)
    diff=100 - (bench/var*100)   

    f,ax= plt.subplots(figsize=(8,8))
    m=rp.get_basemap(resolution='i', lines=False, lllat=latmin, lllon=lonmin, urlat=latmax, urlon=lonmax, ax=ax)
    im=ax.pcolormesh(X,Y, diff, vmin=vmin, cmap=cmap, norm=norm)
    cbar = f.colorbar(im,orientation='horizontal')
    cbar.ax.set_xlabel('%s (%s)' % (d[variable]['abbr'], d[variable]['unit']))
    plt.title('DJF Benchmark - 12.9.3', fontsize=14)
    plt.savefig('/users/mjr583/GC/eval_12.9.3/plots/benchmark-diff_DJF_%s%s.png' % (variable, suff))
    plt.close()

    var=np.mean( var_array[5:8,0,:,:],0 )
    bench=np.mean(benchmark[5:8,0,:,:],0)
    diff=100 - (bench/var *100)

    f,ax= plt.subplots(figsize=(8,8))
    m=rp.get_basemap(resolution='i', lines=False, lllat=latmin, lllon=lonmin, urlat=latmax, urlon=lonmax, ax=ax)
    im=ax.pcolormesh(X,Y, diff, vmin=vmin, cmap=cmap, norm=norm)
    cbar = f.colorbar(im,orientation='horizontal')
    cbar.ax.set_xlabel('%s (%s)' % (d[variable]['abbr'], d[variable]['unit']))
    plt.title('JJA Benchmark - 12.9.3', fontsize=14)
    plt.savefig('/users/mjr583/GC/eval_12.9.3/plots/benchmark-diff_JJA_%s%s.png' % (variable, suff ))
    plt.close()

    var=np.mean( var_array[2:5,0,:,:],0 )
    bench=np.mean(benchmark[2:5,0,:,:],0)
    diff=100 - (bench/var * 100)

    f,ax= plt.subplots(figsize=(8,8))
    m=rp.get_basemap(resolution='i', lines=False, lllat=latmin, lllon=lonmin, urlat=latmax, urlon=lonmax, ax=ax)
    im=ax.pcolormesh(X,Y, diff, vmin=vmin, cmap=cmap, norm=norm)
    cbar = f.colorbar(im,orientation='horizontal')
    cbar.ax.set_xlabel('%s (%s)' % (d[variable]['abbr'], d[variable]['unit']))
    plt.title('MAM Benchmark - 12.9.3', fontsize=14)
    plt.savefig('/users/mjr583/GC/eval_12.9.3/plots/benchmark-diff_MAM_%s%s.png' % (variable, suff ))
    plt.close()

    var=np.mean( var_array[8:11,0,:,:],0 )
    bench=np.mean(benchmark[8:11,0,:,:],0)
    diff=100 - (bench/var * 100)

    f,ax= plt.subplots(figsize=(8,8))
    m=rp.get_basemap(resolution='i', lines=False, lllat=latmin, lllon=lonmin, urlat=latmax, urlon=lonmax, ax=ax)
    im=ax.pcolormesh(X,Y, diff, vmin=vmin, cmap=cmap, norm=norm)
    cbar = f.colorbar(im,orientation='horizontal')
    cbar.ax.set_xlabel('%s (%s)' % (d[variable]['abbr'], d[variable]['unit']))
    plt.title('SON Benchmark - 12.9.3', fontsize=14)
    plt.savefig('/users/mjr583/GC/eval_12.9.3/plots/benchmark-diff_SON_%s%s.png' % (variable, suff ))
    plt.close()

