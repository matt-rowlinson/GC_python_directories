
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('agg')
from matplotlib.colors import LogNorm
import sys
sys.path.append('/users/mjr583/python_lib')
from CVAO_dict import CVAO_dict as d
import RowPy as rp
import GC_tools as GC

rundir='temp'
variable='EmisC2H6_Anthro'
date='201603010000'

ems, times, lat, lon, lev, area, var_list, unit = GC.HEMCO_Diag_read(rundir, version='GEOS-Chem',variable=variable, date=date)
ems=ems[0]

X,Y = np.meshgrid(lon,lat)
f, ax = plt.subplots()
m=rp.get_basemap(ax=ax)

print(X.shape,Y.shape)
print(ems.shape)
if len(ems.shape) == 3:
    ems=ems[0]
print(ems.shape)

im=ax.pcolormesh(X,Y,ems, norm=LogNorm())
cbar = f.colorbar(im,orientation='horizontal')
cbar.ax.set_xlabel('variable')

plt.savefig('./v12test.png')
plt.close()
