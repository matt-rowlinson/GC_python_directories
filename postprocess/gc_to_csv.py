
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

spec='ALD2'

var=[] ; times=[] 
for i,infile in enumerate(sorted(glob.glob('/users/mjr583/scratch/GC/GEOS-Chem/rundirs/fullchem_4x5_LVOCfalse/OutputDir/GEOSChem.SpeciesConc*.nc4'))):
        print(infile) 
        
        fh=Dataset(infile)
        var.append( fh.variables['SpeciesConc_%s' %spec ][:,0,:,:] )
        unit=fh.variables['SpeciesConc_ALD2'].units
        longname=fh.variables['SpeciesConc_ALD2'].long_name

        time=fh.variables['time'][:]
        t0=fh.variables['time'].units
        t0=(int, re.findall(r'\d+', t0))[1]
        t0=datetime.datetime(int(t0[0]), int(t0[1]), int(t0[2]), int(t0[3]), int(t0[4]), int(t0[5]) )
        for dt in time:
            times.append( t0 + datetime.timedelta(minutes=dt) )
lat= fh.variables['lat'][:]
lon= fh.variables['lon'][:]
time=np.array(times)

var=np.concatenate(var)

lat_idx = rp.find_nearest(lat, 16.9)
lon_idx = rp.find_nearest(lon, -24.9)

var=var[:,lat_idx, lon_idx]

df=pd.DataFrame( {'%s (%s)' %(longname, unit): var }, index=time )
df.index.name='datetime'
print(df)
df.to_csv('/users/mjr583/GC/postprocess/csv_files/%s_at_cvao.csv' %spec)

df.plot()
plt.savefig('temp.png')
