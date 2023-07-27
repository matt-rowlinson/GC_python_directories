#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#SBATCH --job-name=NOxtimeseries
#SBATCH --ntasks=1
#SBATCH --mem=18Gb
#SBATCH --partition=nodes
#SBATCH --time=00:45:00
#SBATCH --output=Logs/timeseries_%A.log
import sys
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib
matplotlib.use('agg')
sys.path.append('/users/mjr583/python_lib')
import GC_tools as GC
import RowPy as rp
import numpy as np
from CVAO_dict import CVAO_dict as d
import CVAO_tools as CV
#plt.style.use('seaborn-darkgrid')

inputs=GC.get_arguments()
variable=inputs.var
if variable==None:
    variable='ethane'

rundirs=['fullchem_4x5_LVOCfalse']#,'no_shipping_global','no_shipping_local']#shp_x0.99','shp_x1.01']
string='ACTUALNOSHIP_'

dfs=[]
for a in rundirs:
    print(a)
    var, lat, lon, lev, atime = GC.get_gc_var(rundir=a, variable=variable, version='GEOS-Chem', year='201701')
    
    a_time=[] 
print(var.shape)
var=var[0,0,:,:]

X,Y = np.meshgrid(lon, lat)
print(var.shape)
print(X.shape)

f,ax=plt.subplots()
m = rp.get_basemap(lllat=5., urlat=30., lllon=-40, urlon=-10, ax=ax, lines=False,resolution='i')
#m = rp.get_basemap(ax=ax)
im = m.pcolormesh(X,Y,var, vmin=25, vmax=35)
m.drawcoastlines()

cbar_ax = f.add_axes([0.9, 0.15, 0.05, 0.6])
f.colorbar(im, cax=cbar_ax)#, orientation='horizontal')

plt.savefig('plots/boxes.png')
plt.close()

sys.exit()
