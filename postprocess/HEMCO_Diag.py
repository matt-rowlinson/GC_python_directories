
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
import GC_tools as GC

rundir='boundary_conds_run'
version='12.9.3'
variables=['EmisCO_BioBurn', 'EmisCO_Total']

if rundir == None:
    sys.exit('Please give the run directory')
times=[]


#emissions, time, lat, lon, lev, area = GC.HEMCO_Diag_read(rundir='tropchem_merra_4x5', version=version,variable='EmisCO_Total')
#old_emissions, time, lat, lon, lev, area = GC.HEMCO_Diag_read(rundir=rundir, version=version,variable='EmisCO_Total')
#print(emissions.sum())
#print(old_emissions.sum())

#sys.exit()

for infile in sorted(glob.glob('/users/mjr583/scratch/GC/%s/rundirs/%s/OutputDir/HEMCO_diagnostics.201701010000.nc' % (version, rundir))):
    fh=Dataset(infile)
    time=fh.variables['time'][:]
    lat=fh.variables['lat'][:]
    lon=fh.variables['lon'][:]
    lev=fh.variables['lev'][:]
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

names=list(fh.variables.keys())[8:]
names = [s for s in names if "CO_Anthro" in s]
print(names)
if variable != None:
    names= [s for s in names if variable+'_' in s]
for var in names:  
    print(var)
    mw=False
    gc_em=[]
    for i,infile in enumerate(sorted(glob.glob('/users/mjr583/scratch/GC/%s/%s/output/HEMCO_diag*.nc' % (version, rundir)))):
        fh=Dataset(infile)
        em=fh.variables[var]
        unit=em.units
        long_name=em.long_name

        em=np.array(np.squeeze(em[:]))
        days_in_month=calendar.monthrange(times[i].year, times[i].month)[1]
        #if days_in_month==29:
        #    days_in_month=28
        
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
    #sys.exit()
    ems=np.array(gc_em)
    if ems.ndim == 4:
        ems=np.sum(np.sum(np.sum(ems[:,:,:,:],1),1),1)
    elif ems.ndim == 3:
        ems=np.sum(np.sum(ems,1),1)
    else: # This shouldn't happen
        print('What now?')
        sys.exit()
    df=pd.DataFrame({'Value':ems}, index=times)
    #print('Annual total =', df.groupby(df.index.year).sum())
    #print('All = ', df)
    df.to_csv('./%s.csv' %variable)
    #sys.exit()

    years=df.index.year.drop_duplicates().tolist()[1:]
    for year in years:
        if year==2013:
            continue
        plt.scatter(df[str(year)].index.month, df.Value[str(year)], label=str(year))
    plt.legend()
    #plt.savefig('./test.png')
    #plt.close()
    #sys.exit()


    f,ax= plt.subplots(figsize=(12,4))
    ax.plot(times, ems, 'orchid', label=long_name)
    plt.ylabel(unit)
    plt.legend()
    plt.savefig('./plots/%s_%s.png' % (rundir, long_name) )
    plt.close()

