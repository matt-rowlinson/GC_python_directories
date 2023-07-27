#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#SBATCH --job-name=timeseries
#SBATCH --ntasks=1
#SBATCH --mem=10gb
#SBATCH --partition=interactive
#SBATCH --time=00:10:00
#SBATCH --output=LOGS/timeseries.log
import sys
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib
import glob
import numpy as np
import re
from netCDF4 import Dataset
matplotlib.use('agg')
sys.path.append('/users/mjr583/python_lib')
import GC_tools as GC
import RowPy as rp
from CVAO_dict import CVAO_dict as d
import CVAO_tools as CV
sae = np.swapaxes(rp.surface_area_earth(720, 360, verbose=False),0,1)

years = range(2000, 2020, 1)
years = [ str(xx) for xx in years ]
x = pd.to_datetime( years, format="%Y" )

ann_sum=[] ; ann_em=[]
print('Get v2021')
for y in years:
    inpath = '/mnt/lustre/groups/chem-acm-2018/earth0_data/GEOS/HEMCO/CEDS/v2021-06/%s/' %y
    em=0
    for infile in sorted(glob.glob('%s/*anthro*' %inpath)):
        fh = Dataset(infile)
        try:
            voc_check = fh.VOC_name
            print(voc_check)
        except:
            continue
        keys=list(fh.variables.keys())
        for key in keys[3:]:
            em += fh.variables[key][:]
    em = np.sum(em, 0)
    em = em * sae * (3600*24)*365/12 * 1e-9
    ann_em.append(em)
    ann_sum.append(np.nansum(em))
    sys.exit()
v21=np.array(ann_em)
v21_sums=ann_sum


ann_sum=[] ; ann_em=[]
print('Get v2020')
for y in years[:-2]:
    inpath = '/mnt/lustre/groups/chem-acm-2018/earth0_data/GEOS/HEMCO/CEDS/v2020-08/%s/' %y
    em=0
    for infile in sorted(glob.glob('%s/*anthro*' %inpath)):
        fh = Dataset(infile)
        try:
            voc_check = fh.VOC_name
        except:
            continue
        keys=list(fh.variables.keys())
        for key in keys[3:]:
            em += fh.variables[key][:]
    em = np.sum(em, 0)
    em = em * sae * (3600*24)*365/12 * 1e-9
    ann_em.append(em)
    ann_sum.append(np.nansum(em))
v20=np.array(ann_em)
v20_sums=ann_sum

ann_sum=[] ; ann_em=[]
print('Get v2018')
for y in years[:-5]:
    inpath = '/mnt/lustre/groups/chem-acm-2018/earth0_data/GEOS/HEMCO/CEDS/v2018-08/%s/' %y
    em=0
    for infile in sorted(glob.glob('%s/*anthro*' %inpath)):
        fh = Dataset(infile)
        try:
            voc_check = fh.VOC_name
        except:
            continue
        keys=list(fh.variables.keys())
        for key in keys[3:]:
            em += fh.variables[key][:]
    em = np.sum(em, 0)
    em = em * sae * (3600*24)*365/12 * 1e-9
    ann_em.append(em)
    ann_sum.append(np.nansum(em))
v18=np.array(ann_em)
v18_sums=ann_sum

