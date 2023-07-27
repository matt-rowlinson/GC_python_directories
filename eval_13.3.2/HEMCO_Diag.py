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


def emission_to_tg(ems, area, unit="kg/m2/s"):
    hold=[] 
    for em in ems:
        x = em * area           ## kg/m2/s to kg/s per gridbox
        x = x * 86400 * 30      ## kg/s/gridbox to kg/gridbox/month
        em =  x * 1e-9          ## kg to Tg
        hold.append( em )
    em = np.sum( np.sum( np.array(hold) , 1 ) , 1 )
    return em

a, time, lat, lon, lev, area, var_list,unit = GC.HEMCO_Diag_read(rundir='control', variable="EmisC2H6_Anthro", version='13.1.2', year='2016')

em = np.nansum( a , 0)
print(em.shape, area.shape)
em = em * area * (3600*24)*365/12 * 1e-9
print(em.shape)
print(np.nansum(em))

sys.exit()

em = emission_to_tg( ems=a, area=area )
print(em.shape)
    
