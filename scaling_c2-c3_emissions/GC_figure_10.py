#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#SBATCH --job-name=timeseries
#SBATCH --ntasks=1
#SBATCH --mem=2Gb
#SBATCH --partition=interactive
#SBATCH --time=00:05:00
#SBATCH --output=Logs/timeseries_%A.log
import sys
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib
matplotlib.use('agg')
sys.path.append('/users/mjr583/python_lib')
import GC_tools as GC
import RowPy as rp
from CVAO_dict import CVAO_dict as d
import CVAO_tools as CV
plt.style.use('seaborn-darkgrid')


inputs=GC.get_arguments()
variable=inputs.var
if variable==None:
    variable='ethane'

rundirs=['ceds_only','scale_all_vocs']
df0=[] ; cvs=[]
variables=['ethane','propane','toluene','benzene','PRPE', 'ALK4']
scales=[1e12,1e12,1e9,1e12,1e12,1e12]
for variable in variables:
    print( variable )
    dfs=[]
    for a in rundirs:
        var, lat, lon, lev, atime = GC.get_gc_var(rundir=a, variable=variable, version='13.1.2', year='2017')
        
        a_time=[] 
        y = rp.find_nearest(lat, 16.9)
        x = rp.find_nearest(lon, -24.9)
        for t in range(len(atime)):
            v=var[t,0,y,x] * 1e12
            a_time.append(v)
        a=pd.DataFrame({a:a_time}, index=atime)
        dfs.append(a)
    
    cv_df=pd.concat(dfs, axis=1, ignore_index=False)
    dfs[0]=dfs[0].resample('D').mean()
    labels=['Base','Scaled VOC emissions']
    colors=['#1b9e77','darkorange','#7570b3']
    
    ## Turn into single plot with subplots
    df0.append( cv_df )


f, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2,figsize=(12,6),sharex=False)
import matplotlib.dates as mdates
myFmt = mdates.DateFormatter('%b')
ax1.xaxis.set_major_formatter(myFmt)
plt.gca().xaxis.set_major_locator(mdates.MonthLocator(interval=2))

df = CV.get_from_merge(d['toluene'])
df = df.loc['2017'].resample('D').mean()
ax3.plot( df0[2].ceds_only, c=colors[0], label='Base' )
ax3.plot( df0[2].scale_all_vocs, color=colors[1], label='Scaled CEDS' )
ax3.plot( df, c='k', label='Observed' )
ax3.legend()
ax3.set_ylabel( 'Toluene (pptv)' )
ax3.xaxis.set_major_formatter(myFmt)
plt.gca().xaxis.set_major_locator(mdates.MonthLocator(interval=2))

df = CV.get_from_merge(d['benzene'])
df = df.loc['2017'].resample('D').mean()
ax1.plot( df0[3].ceds_only, c=colors[0], label='Base' )
ax1.plot( df0[3].scale_all_vocs, color=colors[1], label='Scaled CEDS' )
ax1.plot( df, c='k', label='Observed' )
ax1.legend()
ax1.set_ylabel( 'Benzene (pptv)' )
ax1.xaxis.set_major_formatter(myFmt)
plt.gca().xaxis.set_major_locator(mdates.MonthLocator(interval=2))

iB = CV.get_from_merge(d['iso_butane']).fillna(0)
nB = CV.get_from_merge(d['n_butane']).fillna(0)
iP = CV.get_from_merge(d['iso_pentane']).fillna(0)
nP = CV.get_from_merge(d['n_pentane']).fillna(0)
df = iB + nB + iP + nP
df = df.loc['2017'].resample('D').mean()
ax2.plot( df0[5].ceds_only, c=colors[0], label='Base' )
ax2.plot( df0[5].scale_all_vocs, color=colors[1], label='Scaled CEDS' )
ax2.plot( df, c='k', label='Observed' )
ax2.legend()
ax2.set_ylabel( 'ALK4 (pptv)' )
ax2.xaxis.set_major_formatter(myFmt)
plt.gca().xaxis.set_major_locator(mdates.MonthLocator(interval=2))

plt.delaxes( ax4 )
plt.tight_layout()
plt.savefig( 'IGC10/plots/figure_10.png')   
plt.close()

