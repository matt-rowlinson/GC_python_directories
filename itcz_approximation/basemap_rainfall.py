#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#SBATCH --job-name=strmfunction
#SBATCH --ntasks=1
#SBATCH --mem=10gb
#SBATCH --partition=interactive
#SBATCH --time=00:30:00
#SBATCH --output=LOGS/strmfunc.log
import matplotlib as mpl
import matplotlib.pyplot as plt
from netCDF4 import Dataset

mpl.rcParams['mathtext.default'] = 'regular'
import sys
sys.path.append('/users/mjr583/python_lib')
import cartopy_GC_tools as GC
import RowPy as rp
sys.path.append('/users/mjr583/scratch/python_archive')
from get_surface_area import surface_area_earth
import pandas as pd
import numpy as np
from matplotlib.colors import LogNorm

rundir='2018_4x5' ; version='12.9.3' ; jobid=0
month_lengths=[0,31,59,90,120,151,181,212,243,273,304,334,365]
months=['01','02','03','04','05','06','07','08','09','10','11','12']
days_in_month = [31,28,21,30,31,30,31,31,30,31,30,31]
mons=['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
variable='O3'
var, lat, lon, lev, time = GC.get_gc_var(rundir, variable, version)

hold=[]
for x in range(12):
    print(x)
    in_days=pd.date_range(time[month_lengths[x]], time[month_lengths[x+1]-1], freq='D')
    tp, mlat,mlon=GC.get_all_gc_input(in_days, var='PRECTOT', filetype='A1')
    print(np.sum(tp))
    mlat[0]=-90.
    mlat[-1]=90.
    tp=np.sum(tp, 0)

    tp = tp * 86400 * days_in_month[x]
    fltr=np.where(tp <= 5e4)
    tp[fltr]=np.nan

    X,Y=np.meshgrid(mlon, mlat)
    
    # Plot total precip
    #m=rp.get_basemap(lllon=-70, urlon=50, lllat=-10, urlat=75)
    m=rp.get_basemap()
    Max=np.round(np.nanmax(tp),2)
    levels=[1e-1, 1e1, 1e3,5e3,1e4,5e4, 1e5, 5e5]
    colors=['#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99']
    colors=['#ffffd9','#edf8b1','#c7e9b4','#7fcdbb','#41b6c4','#1d91c0','#225ea8','#253494','#081d58']
    tot = m.contourf(X,Y, tp, cmap=plt.cm.Blues)#,levels=levels, extend='min')
    #tot = m.contourf(X,Y, tp, cmap=plt.cm.tab10,levels=levels, norm=LogNorm(), extend='min')
    tot = m.contourf(X,Y, tp, colors=colors,levels=levels,extend='min')#, norm=LogNorm(), extend='min')
    #tot = m.pcolormesh(X,Y,tp, cmap=plt.cm.Blues, vmin=.1)
    tot.cmap.set_under('w')
    cbar=plt.colorbar(tot,orientation='horizontal')
    cbar.set_label('Total Precipitation (mm)')
    plt.title('%s 2018' %mons[x], fontsize=16)
    plt.savefig('plots/basemap_total_precipitation_%s.png' %months[x])
    plt.close()
