#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#SBATCH --job-name=timeseries
#SBATCH --ntasks=1
#SBATCH --mem=10gb
#SBATCH --partition=interactive
#SBATCH --time=00:10:00
#SBATCH --output=LOGS/timeseries.log
import sys
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib
import glob
import numpy as np
import re
from netCDF4 import Dataset
matplotlib.use('agg')
sys.path.append('/users/mjr583/python_lib')
import GC_tools as GC
import RowPy as rp
from CVAO_dict import CVAO_dict as d
import CVAO_tools as CV
#plt.style.use('seaborn-darkgrid')

inpath = '/mnt/lustre/users/mjr583/GC/emissions/'
variable='NMVOC'
for infile in sorted(glob.glob('%s/%s*CMIP_CEDS*' %(inpath, variable))):
    print(infile)

    fh = Dataset(infile)
    lat=fh.variables['lat'][:]
    lon=fh.variables['lon'][:]
    sectors=fh.variables['sector'][:]
    keys=list(fh.variables.keys())

    em = fh.variables[keys[4]]
    
    time=fh.variables['time']
    t = re.findall(r'\b\d+\b', time.units)
    t = [ int(x) for x in t ] 
    dt = rp.days_to_datetime(time, start_year=t[0], start_month=t[1], start_day=t[2], start_hour=t[3])

em = np.sum(em, 1)
sae = np.swapaxes(rp.surface_area_earth(720, 360, verbose=False),0,1)
v=[]
for i in range(len(em)):
    Em=em[i,:,:]
    Em = Em * sae * (60*60*24*365/12) * 1e-9
    v.append(np.round(Em.sum(),2))

v = np.sum( np.array(v).reshape(-1, 12), axis=1 )
index = pd.to_datetime( dt )
x = index.year.unique()[1:]
x = pd.to_datetime( x, format="%Y" )
print(v[-1])

variable='VOC02'
for infile in sorted(glob.glob('%s/%s*CMIP_CEDS*' %(inpath, variable))):
    print(infile)

    fh = Dataset(infile)
    lat=fh.variables['lat'][:]
    lon=fh.variables['lon'][:]
    sectors=fh.variables['sector'][:]
    keys=list(fh.variables.keys())

    em = fh.variables[keys[4]]
    
    time=fh.variables['time']
    t = re.findall(r'\b\d+\b', time.units)
    t = [ int(x) for x in t ] 
    dt = rp.days_to_datetime(time, start_year=t[0], start_month=t[1], start_day=t[2], start_hour=t[3])

em = np.sum(em, 1)
v2=[]
for i in range(len(em)):
    Em=em[i,:,:]
    Em = Em * sae * (60*60*24*365/12) * 1e-9
    v2.append(np.round(Em.sum(),2))

v2 = np.sum( np.array(v2).reshape(-1, 12), axis=1 )
index = pd.to_datetime( dt )
x = index.year.unique()[1:]
x = pd.to_datetime( x, format="%Y" )

if variable=='VOC02':
    variable='C2H6'

print('Funky plot')
fig,ax = plt.subplots(figsize=(9,5))
fig.subplots_adjust(right=0.75)

twin1 = ax.twinx()
twin2 = ax.twinx()
twin2.spines['right'].set_position(("axes", 1.2))

p1, = ax.plot(  x, v , label='Total NMVOC' )
p2, = twin1.plot( x, v2, label='C2H6', color='r' )
p3, = twin2.plot( x, v2 / v, label='Ratio', color='g' )

ax.yaxis.label.set_color(p1.get_color())
twin1.yaxis.label.set_color(p2.get_color())
twin2.yaxis.label.set_color(p3.get_color())

ax.set_ylabel( 'NMVOCs (Tg $yr^{-1}$)'  )
twin1.set_ylabel( 'C2H6 (Tg $yr^{-1}$)'   )
twin2.set_ylabel( 'C2H6/NMVOCs Ratio'  )

tkw = dict(size=4, width=1.5)
ax.tick_params(axis='y', colors=p1.get_color(), **tkw)
twin1.tick_params(axis='y', colors=p2.get_color(), **tkw)
twin2.tick_params(axis='y', colors=p3.get_color(), **tkw)
ax.tick_params(axis='x', **tkw)

ax.legend(handles=[p1,p2,p3])
plt.savefig( 'plots/GLOBUS_CEDS_2019_both.png' )
plt.close()

sys.exit()
