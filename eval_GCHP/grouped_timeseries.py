#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#SBATCH --job-name=timeseries
#SBATCH --ntasks=1
#SBATCH --mem=100Mb
#SBATCH --partition=nodes
#SBATCH --time=00:10:00
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
from CVAO_dict import CVAO_dict as d
import CVAO_tools as CV
plt.style.use('seaborn-darkgrid')

inputs=GC.get_arguments()
variable=inputs.var
if variable==None:
    variable='O3'

from GC_tools import Iy as D
variable='Iy'
a, lat, lon, lev, atime = GC.get_var_group(D, rundir='tropchem_merra_4x5', version='12.9.3', year='2016')
b, lat,lon,lev,btime = GC.get_var_group(D, rundir='fullchem_4x5_LVOCfalse', version='GEOS-Chem',year='2016')


delta,interval=GC.find_timestep(atime)
y = rp.find_nearest(lat, 16.9)
x = rp.find_nearest(lon, -24.9)

a_time=[] ; b_time=[]
for t in range(len(atime)):
    v=a[t,0,y,x]
    a_time.append(v)

for t in range(len(btime)):
    v=b[t,0,y,x]
    b_time.append(v)

a=pd.DataFrame({'v12.9.3':a_time}, index=atime)
a=a.tz_localize("Atlantic/Cape_Verde")
b=pd.DataFrame({'v13.0.0':b_time}, index=btime)
b=b.tz_localize("Atlantic/Cape_Verde")
a=a[b.index[0]:b.index[-1]]

rp.timeseries_plot([a,b], variable, cvao_data=False)
rp.timeseries_plot([a,b], variable, cvao_data=False, std=True)
