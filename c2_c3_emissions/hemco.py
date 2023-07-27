
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
import calendar
import argparse
import sys
sys.path.append('/users/mjr583/python_lib')
from CVAO_dict import CVAO_dict as d
import RowPy as rp

rundir='geo_ethane_NEI'
version='12.9.3'
variable=None

if rundir == None:
    sys.exit('Please give the run directory')
times=[]
for infile in sorted(glob.glob('/users/mjr583/scratch/GC/GEOS-Chem/rundirs/%s/OutputDir/HEMCO_diag*15*.nc' % (rundir))):
    fh=Dataset(infile)
    print(infile)
    
    time=fh.variables['time'][:]
    AREA=fh.variables['AREA'][:]

    t0=fh.variables['time'].units
    t0=(int, re.findall(r'\d+', t0))[1]
    t0=datetime.datetime(int(t0[0]), int(t0[1]), int(t0[2]), int(t0[3]), int(t0[4]), int(t0[5]) )
    for dt in time:
        times.append( t0 + datetime.timedelta(minutes=dt) )

times=np.array(times)
lat=fh.variables['lat'][:]
lon=fh.variables['lon'][:]
lev=fh.variables['lev'][:]

NMVOCs = ["BENZ","TOLU","ALK4","C2H6","MONX","MTPO","C2H4","C3H6","C3H8", "PRPE","ALK2","MEOH","ACET", "XYLE", "LVOC", "LVOCOA","MEK", "ISOP", "ROH"]
names=list(fh.variables.keys())[8:]
names = [s for s in names if "Anthro" in s]
names = [s for s in names if s.replace("Emis","").replace("_Anthro","") in NMVOCs]

total_nmvocs=0
if variable != None:
    names= [s for s in names if variable+'_' in s]
for var in names:  
    gc_em=[]
    for i,infile in enumerate(sorted(glob.glob('/users/mjr583/scratch/GC/GEOS-Chem/rundirs/%s/OutputDir/HEMCO_diag*15*.nc' %rundir))):
        fh=Dataset(infile)
        em=fh.variables[var]
        unit=em.units
        long_name=em.long_name

        em=np.array(np.squeeze(em[:]))
        days_in_month=calendar.monthrange(times[i].year, times[i].month)[1]
        
        unit='Tg'
        if len(em.shape) == 3:  ## Convert from kg/m2/s to Tg month-1
            hold_em=[]
            for i in range(47):
                x = em[i] * AREA             ## kg/m2/s to kg/s per gridbox
                x = x * 86400 * days_in_month           ## kg/s/gridbox to kg/gridbox/month
                hold_em.append( x * 1e-9)    ## kg to Tg
            em=np.array(hold_em)
        else:
            x = em * AREA           ## kg/m2/s to kg/s per gridbox
            x = x * 86400 * 30      ## kg/s/gridbox to kg/gridbox/month
            em =  x * 1e-9          ## kg to Tg
        gc_em.append(em)
    ems=np.array(gc_em)
    ems=np.sum(np.sum(ems, 0),0)

    import read as R

    tot = R.regional_totals(ems, lat, lon, region="North America")
    if "C2H6" in var:
        c2h6=tot.sum()
    #print(tot.sum())
    total_nmvocs += tot.sum()
print("TOTAL = ", total_nmvocs)
print("Ratio = ", c2h6 / total_nmvocs )

print( c2h6 / 21.76 )

'''
    if ems.ndim == 4:
        ems=np.sum(np.sum(np.sum(ems[:,:,:,:],1),1),1)
    elif ems.ndim == 3:
        ems=np.sum(np.sum(ems,1),1)
    else: # This shouldn't happen
        print('What now?')
        sys.exit()
    df=pd.DataFrame({'Value':ems}, index=times)
    print('Annual total =', df.groupby(df.index.year).sum())
'''
