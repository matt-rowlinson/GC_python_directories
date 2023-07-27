#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#SBATCH --job-name=NMVOCs
#SBATCH --ntasks=1
#SBATCH --mem=8gb
#SBATCH --partition=interactive
#SBATCH --time=00:10:00
#SBATCH --output=Logs/NMVOCs.log
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

inpath = '/mnt/lustre/users/mjr583/GC/emissions/'
variable='CO'

for infile in sorted(glob.glob('%s/%s-em-anthro*CMIP_CEDS*' %(inpath, variable))):
    print(infile)

    fh = Dataset(infile)
    lat=fh.variables['lat'][:]
    lon=fh.variables['lon'][:]
    sectors=fh.variables['sector'][:]
    keys=list(fh.variables.keys())

    em = fh.variables[keys[4]]
    
    time=fh.variables['time']
    t = re.findall(r'\b\d+\b', time.units)
    t = [ int(x) for x in t ] 
    dt = rp.days_to_datetime(time, start_year=t[0], start_month=t[1], start_day=t[2], start_hour=t[3])
NMVOC = np.sum(em, 1)

v=[]
for i in range(len(em)):
    Em=NMVOC[i,:,:]
    Em = Em * sae * (60*60*24*365/12) * 1e-9
    v.append(np.round(Em.sum(),2))
NMVOC=v

for infile in sorted(glob.glob('%s/%s*BIO*CMIP_CEDS*' %(inpath, variable))):
    print(infile)

    fh = Dataset(infile)
    lat=fh.variables['lat'][:]
    lon=fh.variables['lon'][:]
    keys=list(fh.variables.keys())
    bf = fh.variables[keys[4]]
bf = np.sum(bf, 1)
v=[]
for i in range(len(bf)):
    Em=bf[i,:,:]
    Em = Em * sae * (60*60*24*365/12) * 1e-9
    v.append(np.round(Em.sum(),2))

for infile in sorted(glob.glob('%s/%s*AIR*CMIP_CEDS*' %(inpath, variable))):
    print(infile)

    fh = Dataset(infile)
    lat=fh.variables['lat'][:]
    lon=fh.variables['lon'][:]
    keys=list(fh.variables.keys())
    air = fh.variables[keys[4]]
air = np.sum(air, 1)
Air=[]
for i in range(len(air)):
    Em=air[i,:,:]
    Em = Em * sae * (60*60*24*365/12) * 1e-9
    Air.append(np.round(Em.sum(),2))
NMVOC_bf= np.array(v) + np.array(NMVOC)

NMVOC = np.sum( np.array(NMVOC).reshape(-1, 12), axis=1 )
NMVOC_bf = np.sum( np.array(NMVOC_bf).reshape(-1, 12), axis=1 )

years = range(2000, 2020, 1)
years = [ str(xx) for xx in years ]
x = pd.to_datetime( years, format="%Y" )

older=[]
for y in years:
    print(y)
    inpath = '/mnt/lustre/groups/chem-acm-2018/earth0_data/GEOS/HEMCO/CEDS/v2021-06/%s/' %y
    em=0
    for infile in sorted(glob.glob('%s/CO-*' %inpath)):
        fh = Dataset(infile)
        keys=list(fh.variables.keys())
        for key in keys[3:]:
            em += fh.variables[key][:]
    em = np.sum(em, 0)
    em = em * sae * (3600*24)*365/12 * 1e-9
    older.append(np.nansum(em))
v21=older #+ v

older=[]
for y in years[:-2]:
    inpath = '/mnt/lustre/groups/chem-acm-2018/earth0_data/GEOS/HEMCO/CEDS/v2020-08/%s/' %y
    em=0
    for infile in sorted(glob.glob('%s/CO-*' %inpath)):
        print(infile)
        fh = Dataset(infile)
        keys=list(fh.variables.keys())
        for key in keys[3:]:
            em += fh.variables[key][:]
    em = np.sum(em, 0)
    em = em * sae * (3600*24)*365/12 * 1e-9
    older.append(np.nansum(em))
v20=older#[:-2]


older=[]
for y in years[:-5]:
    inpath = '/mnt/lustre/groups/chem-acm-2018/earth0_data/GEOS/HEMCO/CEDS/v2018-08/%s/' %y
    em=0
    for infile in sorted(glob.glob('%s/CO-*' %inpath)):
        print(infile)
        fh = Dataset(infile)
        keys=list(fh.variables.keys())
        for key in keys[3:]:
            em += fh.variables[key][:]
    em = np.sum(em, 0)
    em = em * sae * (3600*24)*365/12 * 1e-9
    older.append(np.nansum(em))
v18=older

fig,ax = plt.subplots(figsize=(9,5))

p3, = ax.plot(  x[:len(v18)], v18 , c='#4daf4a', label='v2017-05-18' )
p2, = ax.plot(  x[:len(v20)], v20 , c='#377eb8', label='v2020_v1.0')
p5, = ax.plot(  x,            NMVOC   , c='#984ea3', label='v2021-04-21 CO-em-anthro' )
p4, = ax.plot(  x,            NMVOC_bf, c='#e41a1c',label='v2021-04-21 CO-em-anthro + CO-em-SOLID-BIOFUEL')#+ NMVOC-em-AIR' )
p1, = ax.plot(  x, v21 , ls='-.', c='yellow', label='v2021-04-21 HEMCO CO-em-anthro' )

ax.set_ylabel( 'Anthropogenic NMVOCs (Tg $yr^{-1}$)'  )
ax.set_ylim(bottom=0.)
ax.legend()
plt.savefig( 'plots/full_CEDS_versions_CO.png' )
plt.close()

