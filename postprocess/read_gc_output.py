#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('agg')
from netCDF4 import Dataset
import re
import datetime     
import glob
from mpl_toolkits.basemap import Basemap

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

#------------------MAIN SCRIPT-----------------------------------------------
path='/users/mjr583/scratch/GC/12.9.3/rundirs/global_4x5/OutputDir/'
variable=[] ; datetimes=[]

for infile in sorted(glob.glob('%s/GEOSChem.SpeciesConc.*.nc4' %path)):
    fh=Dataset(infile)
    
    var=fh.variables['SpeciesConc_CO'][:] 
    variable.append(var)

    time=fh.variables['time'] 
    datetimes.append( GC_time_to_datetime(time) ) 

lat=fh.variables['lat'][:]
lon=fh.variables['lon'][:]
lev=fh.variables['hyam'][:]

variable=np.concatenate(variable)  
variable=variable * 1e9 
time=np.concatenate(datetimes)

### Example 1: Quick timeseries plot ###
f = plt.figure()
plt.plot(time, variable[:,0,27,31], color='r')  
plt.ylabel('CO (ppb)')
plt.savefig('./timeseries.png')
plt.close()

### Example 2: Quick plot using basemap ###
f=plt.figure
m = Basemap(projection='cyl', llcrnrlon=-180,urcrnrlon=180,llcrnrlat=-90,urcrnrlat=90)
m.drawcoastlines()
X,Y = np.meshgrid(lon,lat)
mean_variable = np.mean(variable,0)            
m.pcolormesh(X,Y, mean_variable[0,:,:], cmap='viridis')
plt.savefig('./map.png')
plt.close()
