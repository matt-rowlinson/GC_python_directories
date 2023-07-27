#!/usr/bin/env python -W ignore::DeprecationWarning
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

a, time, lat, lon, lev, area, var_list,u = GC.HEMCO_Diag_read(rundir='tropchem_merra_4x5', variable='EmisCO_BioBurn', version='12.9.3', year='201601')
for em in var_list[::-1]:
    if "sC2H6" not in em:
        continue
    print(em)

    a, time, lat, lon, lev, area, var_list,unit = GC.HEMCO_Diag_read(rundir='tropchem_merra_4x5', variable=em, version='12.9.3', year='2016')
    try:
        b, time, lat, lon, lev, area, v,u = GC.HEMCO_Diag_read(rundir='fullchem_hourly_default', variable=em, version='GEOS-Chem',year='2016')
    except:
        continue
    
    secs= (60 * 60) * 24 * 365
    x=a.shape
    print(x)
    if len(x)==4:
        a_2d=np.sum(np.sum(a,0),0) * area * secs * 1e-9
        b_2d=np.sum(np.sum(b,0),0) * area * secs * 1e-9

        a_1d=np.sum(np.sum(np.sum(a,1),1),1) * secs/12 * np.sum(area) * 1e-9 * 1e-3 
        b_1d=np.sum(np.sum(np.sum(b,1),1),1) * secs/12 * np.sum(area) * 1e-9 * 1e-3

    elif len(x)==3:
        a_2d=np.sum(a,0) * area * secs * 1e-9
        b_2d=np.sum(b,0) * area * secs * 1e-9

        a_1d=np.sum(np.sum(a,1),1) * secs/12 * np.sum(area) * 1e-9 * 1e-3
        b_1d=np.sum(np.sum(b,1),1) * secs/12 * np.sum(area) * 1e-9 * 1e-3
    print(a_1d.sum())
    print(b_1d.sum())
    rp.difference_map(a_2d,b_2d, em, lon, lat, labels=['v12.9.3','v13.0.0','v12 - v13'], HEMCO=True, unit='Tg')
    rp.difference_map(a_2d,b_2d, em, lon, lat, labels=['v12.9.3','v13.0.0','v12 - v13c'], pc=True,HEMCO=True,unit='Tg')
    
    a_1d = pd.DataFrame({'v12.9.3 : %s Tg' %int(np.sum(a_1d)) :a_1d}, index=time)
    b_1d = pd.DataFrame({'v13.0.0 : %s Tg' %int(np.sum(b_1d)) :b_1d}, index=time)

    rp.timeseries_plot([a_1d, b_1d], em, cvao_data=False, HEMCO=True,unit='Annual global emission (Tg)')
