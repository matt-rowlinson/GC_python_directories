
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('agg')
from scipy.optimize import curve_fit
from sklearn.metrics import mean_squared_error, r2_score
from matplotlib.offsetbox import AnchoredText
import netCDF4
from netCDF4 import Dataset
import datetime     
import os
import glob
import re
from matplotlib.colors import LogNorm
import calendar
import argparse
import sys
sys.path.append('/users/mjr583/python_lib')
from CVAO_dict import CVAO_dict as d
import RowPy as rp
import GC_tools as GC

rundir='boundary_conds_run'
version='12.9.3'
em='EmisCO_BioBurn'

control, time, lat, lon, lev, area = GC.HEMCO_Diag_read(rundir=rundir, version=version,variable=em)
print(control.shape)
control = control[0,:,:]
print(control.shape)

X,Y=np.meshgrid(lon,lat)
f,ax1 = plt.subplots(1,1,figsize=(7,6))
m1=rp.get_basemap(ax=ax1)

im1=m1.pcolormesh(X,Y,control, norm=LogNorm())

ax1.set_title('%s - %s' %(rundir, em))

plt.savefig('./plots/HEMCODiag_%s_%s.png' %(rundir, em))
plt.close()
