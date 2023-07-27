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
sae = np.swapaxes(rp.surface_area_earth(720, 360),0,1)
v=[]
for i in range(len(em)):
    Em=em[i,:,:]
    Em = Em * sae * (60*60*24*365/12) * 1e-9
    v.append(np.round(Em.sum(),2))

v = np.sum( np.array(v).reshape(-1, 12), axis=1 )
index = pd.to_datetime( dt )
x = index.year.unique()[1:]
x = pd.to_datetime( x, format="%Y" )

if variable=='VOC02':
    variable='C2H6'
print(v[0], v[-1])
f,ax = plt.subplots()
ax.plot(x, v)
ax.set_ylabel( '%s (Tg $yr^{-1}$)' %variable )
plt.savefig( 'plots/GLOBUS_CEDS_2019_%s.png' %variable )
plt.close()

sys.exit()
