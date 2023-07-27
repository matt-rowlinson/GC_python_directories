import pandas as pd
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import numpy as np
import sys
sys.path.append('/users/mjr583/python_lib/')
import RowPy as rp
import GC_tools as GC
import CVAO_tools as CV
from CVAO_dict import CVAO_dict as d


## Scaling file for GC
f='monthly.gridded.surface.methane.1979-2020.1x1.nc'
fh=Dataset(f)
lat=fh.variables['lat'][:]
lon=fh.variables['lon'][:]
time=fh.variables['time']
ch4=fh.variables['SFC_CH4']


ch4=ch4[-1]
f,ax=plt.subplots(1,1)
m=rp.get_basemap(ax=ax, lllat=-30, urlat=60.)
X,Y=np.meshgrid(lon,lat)

levels=np.arange(1930, 2000, 5)
im=m.contourf(X,Y,ch4)#,levels=levels)
cbar = f.colorbar(im,orientation='horizontal')
cbar.ax.set_xlabel('$CH_4$ (ppb)' )

plt.savefig('./plots/ch4_latbands.png')
plt.close()
