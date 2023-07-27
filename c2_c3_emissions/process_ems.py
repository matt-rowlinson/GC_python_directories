#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#SBATCH --job-name=timeseries
#SBATCH --ntasks=1
#SBATCH --mem=10gb
#SBATCH --partition=interactive
#SBATCH --time=00:10:00
#SBATCH --output=LOGS/timeseries.log
import sys
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib
from netCDF4 import Dataset
import numpy as np
matplotlib.use('agg')
sys.path.append('/users/mjr583/python_lib')
import RowPy as rp
import read as R
plt.style.use('seaborn-darkgrid')

fh = Dataset('/mnt/lustre/groups/chem-acm-2018/earth0_data/GEOS/HEMCO/CEDS/v2021-06/2017/C2H6-em-anthro_CMIP_CEDS_2017.nc')
lat = fh.variables['lat'][:]
lon = fh.variables['lon'][:]

v21_voc = np.load('temp_files/v21_voc.npy')
v20_voc = np.load('temp_files/v20_voc.npy')
v18_voc = np.load('temp_files/v18_voc.npy')

v21_ethane = np.load('temp_files/v21_ethane.npy')
v20_ethane = np.load('temp_files/v20_ethane.npy')
v18_ethane = np.load('temp_files/v18_ethane.npy')

## Sum over first dimension (years)


years = range(2000, 2020, 1)
years = [ str(xx) for xx in years ]
x = pd.to_datetime( years, format="%Y" )

fig,ax = plt.subplots(figsize=(9,5))

p1, = ax.plot(  x[:len(v18_ethane)], v18_ethane / v18_voc, ls='-', label='CEDS v2018 (to 2014)' )
p2, = ax.plot(  x[:len(v20_ethane)], v20_ethane / v20_voc, ls='--', label='CEDS v2020 (to 2017)' )
p3, = ax.plot(  x, v21_ethane / v21_voc, ls=':', label='CEDS v2021 (to 2019)' )
#ax.legend(handles=[L1a,L1b,L1c])

ax.legend()
ax.set_ylabel( 'C2H6 / NMVOCs Ratio'  )
plt.savefig( 'plots/global_ratio.png' )
plt.close()


