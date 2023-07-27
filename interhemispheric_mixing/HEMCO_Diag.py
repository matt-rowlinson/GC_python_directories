
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

rundir='old_ceds_tropchem_merra_4x5'
version='12.9.3'
variables=['EmisCO_BioBurn', 'EmisCO_Total']

if rundir == None:
    sys.exit('Please give the run directory')
times=[]

em='EmisCO_BioBurn'

#control, time, lat, lon, lev, area = GC.HEMCO_Diag_read(rundir='irma_025x03125', version=version,variable=em)
control, time, lat, lon, lev, area = GC.HEMCO_Diag_read(rundir='boundary_conds_run', version=version,variable=em)
noAf, time, lat, lon, lev, area = GC.HEMCO_Diag_read(rundir='irma_025x03125_noAfBB', version=version,variable=em)
noBB, time, lat, lon, lev, area = GC.HEMCO_Diag_read(rundir='irma_025x03125_noBB', version=version,variable=em)
x=control.shape
n=0
if len(x)==4:
    control=control[n,0]*1e9
    noBB=noBB[n,0]*1e9
    noAf=noAf[n,0]*1e9
elif len(x)==3:
    control=control[n]*1e9
    noBB=noBB[n]*1e9
    noAf=noAf[n]*1e9

X,Y=np.meshgrid(lon,lat)

f,(ax1,ax2,ax3) = plt.subplots(1,3,figsize=(12,3))
m1=rp.get_basemap(lllon=-70, urlon=0., lllat=-15., urlat=30.,ax=ax1)
m2=rp.get_basemap(lllon=-70, urlon=0., lllat=-15., urlat=30.,ax=ax2)
m3=rp.get_basemap(lllon=-70, urlon=0., lllat=-15., urlat=30.,ax=ax3)

im1=m1.pcolormesh(X,Y,control, vmax=25., norm=LogNorm())
im2=m2.pcolormesh(X,Y,noAf, vmax=25.,norm=LogNorm())
im3=m3.pcolormesh(X,Y,noBB, vmax=25.,norm=LogNorm())

ax1.set_title('Control')
ax2.set_title('No African biomass burning')
ax3.set_title('No biomass burning')

plt.savefig('./plots/HEMCO_Diag_map_%s.png' %em)
plt.close()
