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
    # Get the reference date (T0) from unit information
    t0=fh.variables['time'].units
    t0=(int, re.findall(r'\d+', t0))[1]
    # Makes reference date a datetime object
    t0=datetime.datetime(int(t0[0]), int(t0[1]), int(t0[2]), int(t0[3]), int(t0[4]), int(t0[5]) )
    # Create empty list then loop over every GC time to get a datetime
    times=[]
    for dt in time[:]:
        times.append( t0 + datetime.timedelta(minutes=dt) )
    # Convert list back to an array and return the array
    times=np.array(times)
    return times

#------------------MAIN SCRIPT-----------------------------------------------
# Set path to GEOS-Chem output and create empty lists to store arrays
path='/users/mjr583/scratch/GC/12.9.3/rundirs/global_4x5/OutputDir/'
variable=[] ; datetimes=[]

# glob.glob finds all files which match the string (important to use sorted() as otherwise order will be random)
for infile in sorted(glob.glob('%s/GEOSChem.SpeciesConc.*.nc4' %path)):
    # Read netcdf using Dataset
    print(infile)
    fh=Dataset(infile
    
    # Read desired species and append to list
    var=fh.variables['SpeciesConc_CO'][:] 
    variable.append(var)

    # Read time variable and pass to function defined above to get datetime array
    time=fh.variables['time'] 
    datetimes.append( GC_time_to_datetime(time) 

# Get these variables out of the loop as will be the same in every file
lat=fh.variables['lat'][:]
lon=fh.variables['lon'][:]
lev=fh.variables['hyam'][:]

# Convert lists to arrays and convert from mol mol-1 to ppb
variable=np.concatenate(variable)  
variable=variable * 1e9 
time=np.concatenate(datetimes)

# Print shape to check dimensions match time array (variable.shape should be (time,levs,lats,lon)
print(variable.shape)
print(time.shape)


### Example 1: Quick timeseries plot ###
# This plots variable for all timesteps, at lev[0], lat[27] and lon[31], i.e. ~Cape Verde
f = plt.figure()
plt.plot(time, variable[:,0,27,31], color='r')  
plt.ylabel('CO (ppb)')
plt.savefig('./timeseries.png')
plt.close()

### Example 2: Quick plot using basemap ###
# Build the map
f=plt.figure
m = Basemap(projection='cyl', llcrnrlon=-180,urcrnrlon=180,llcrnrlat=-90,urcrnrlat=90)
m.drawcoastlines()
X,Y = np.meshgrid(lon,lat)
# Average variable across timesteps then plot on map for lowest model level
mean_variable = np.mean(variable,0)            
m.pcolormesh(X,Y, mean_variable[0,:,:], cmap='viridis'
plt.savefig('./map.png')
plt.close()
