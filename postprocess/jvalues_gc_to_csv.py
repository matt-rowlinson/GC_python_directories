
import glob
from netCDF4 import Dataset
import numpy as np
import pandas as pd
import sys
import datetime
import re
import matplotlib.pyplot as plt

sys.path.append('/users/mjr583/python_lib')
import GC_tools as GC
import RowPy as rp

jno2=[] ; jhono=[] ; jhno3=[] ; times=[] 
for i,infile in enumerate(sorted(glob.glob('/users/mjr583/scratch/GC/12.9.3/simone_Jvalues/GEOSChem.JV*2019*.nc4'))):
        print(infile) 
        fh=Dataset(infile)
        jno2.append( fh.variables['Jval_NO2'][:,0,:,:] )
        jhono.append(fh.variables['Jval_HNO2'][:,0,:,:])
        jhno3.append(fh.variables['Jval_HNO3'][:,0,:,:])  

        time=fh.variables['time'][:]
        t0=fh.variables['time'].units
        t0=(int, re.findall(r'\d+', t0))[1]
        t0=datetime.datetime(int(t0[0]), int(t0[1]), int(t0[2]), int(t0[3]), int(t0[4]), int(t0[5]) )
        for dt in time:
            times.append( t0 + datetime.timedelta(minutes=dt) )
lat= fh.variables['lat'][:]
lon= fh.variables['lon'][:]
time=np.array(times)

jno2=np.concatenate(jno2)
jhono=np.concatenate(jhono)
jhno3=np.concatenate(jhno3)

lat_idx = rp.find_nearest(lat, 16.9)
lon_idx = rp.find_nearest(lon, -24.9)

jno2=jno2[:,lat_idx, lon_idx]
jhono=jhono[:,lat_idx, lon_idx]
jhno3=jhno3[:,lat_idx, lon_idx]

df=pd.DataFrame( {'jNO2':jno2, 'jHONO':jhono, 'jHNO3':jhno3}, index=time )
print(df)
df.to_csv('/users/mjr583/scratch/GC/12.9.3/simone_Jvalues/2019_JValues.csv')
