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
import glob
import re
import datetime
import numpy as np
from netCDF4 import Dataset
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

# Define a function to process time from "minutes since T0 to datetime format 
def GC_time_to_datetime(time):
    t0=fh.variables['time'].units
    t0=(int, re.findall(r'\d+', t0))[1]
    t0=datetime.datetime(int(t0[0]), int(t0[1]), int(t0[2]), int(t0[3]), int(t0[4]), int(t0[5]) )
    times=[]
    for dt in time[:]:
        times.append( t0 + datetime.timedelta(minutes=dt) )
    times=np.array(times)
    return times

variables=['CO','O3','C2H6','NO2','NO']
labels=['Control','No BB','No African BB']
bc = '/users/mjr583/scratch/GC/12.9.3/rundirs/irma_025x03125/GC_BC/'
bc_nobb = '/users/mjr583/scratch/GC/12.9.3/rundirs/irma_025x03125_noBB/GC_BC/'
bc_noaf = '/users/mjr583/scratch/GC/12.9.3/rundirs/irma_025x03125_noAfBB/GC_BC/'

x = -24.9
y = 16.9
CO=[] ; times=[]

for var in variables:
    print(var)
    CO=[] ; times=[]
    for infile in sorted(glob.glob('%s*Boundary*201708*control.nc4' %bc)):
        print(infile)
        fh=Dataset(infile)
        co = fh.variables['SpeciesBC_%s' %var][:]*1e9
        CO.append(co)

        time=fh.variables['time']
        times.append( GC_time_to_datetime(time) ) 

    lats=fh.variables['lat'][:]
    lons=fh.variables['lon'][:]
    lat_idx=rp.find_nearest(lats,y)
    lon_idx=rp.find_nearest(lons,x)
    print(lats[lat_idx])
    print(lons[lon_idx])

    CO=np.concatenate(CO)
    CO=CO[:,0,lat_idx, lon_idx]

    time=np.concatenate(times)

    df=pd.DataFrame({'BC_CO':CO}, index=time)

    CO=[] ; times=[]
    for infile in sorted(glob.glob('%s*Boundary*201708*NOBB.nc4' %bc_nobb)):
        print(infile)
        fh=Dataset(infile)
        co = fh.variables['SpeciesBC_%s' %var][:]*1e9
        CO.append(co)

        time=fh.variables['time']
        times.append( GC_time_to_datetime(time) ) 
    CO=np.concatenate(CO)
    print(CO.shape)
    CO=CO[:,0,lat_idx, lon_idx]
    time=np.concatenate(times)
    nobb=pd.DataFrame({'BC_CO':CO}, index=time)  

    CO=[] ; times=[]
    for infile in sorted(glob.glob('%s*Boundary*201708*noAfBB.nc4' %bc_noaf)):
        print(infile)
        fh=Dataset(infile)
        co = fh.variables['SpeciesBC_%s' %var][:]*1e9
        CO.append(co)

        time=fh.variables['time']
        times.append( GC_time_to_datetime(time) ) 
    CO=np.concatenate(CO)
    print(CO.shape)
    CO=CO[:,0,lat_idx, lon_idx]
    time=np.concatenate(times)
    noaf=pd.DataFrame({'BC_CO':CO}, index=time)
    
    
    nobb=nobb.resample('D').mean()
    noaf=noaf.resample('D').mean()
    df=df.resample('D').mean()
    df=df['2017-08':'2017-09']


    plt.plot(df.index,df.BC_CO,label='Default boundary conditions')
    plt.plot(nobb.index,nobb.BC_CO, label='With no biomass burning')
    plt.plot(noaf.index,noaf.BC_CO, label='With no African biomass burning')

    plt.ylim(bottom=0)
    plt.legend()
    plt.savefig('plots/GC_BC_%s_at%sN%sW.png' %(var, str(lats[lat_idx]), str(lons[lon_idx]) ) )
    plt.close()
