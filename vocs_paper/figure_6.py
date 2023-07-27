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

rundirs=['ceds_only','geo_only','all_scaled']
df0=[] ; cvs=[]
variables=['ethane','propane','toluene','benzene','PRPE','ALK4']
for variable in variables:
    print( variable )
    dfs=[]
    for a in rundirs:
        var, lat, lon, lev, atime = GC.get_gc_var(rundir=a, variable=variable, version='13.4.0', year='2017')
        
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
    labels=['Base','No VOC emissions']
    colors= ['y','#0000ff', '#ff0000']
    #['#1b9e77','darkorange','#7570b3']
    
    ## Turn into single plot with subplots
    df0.append( cv_df )

f, ((ax1,ax2,ax3),(ax4,ax5,ax6)) = plt.subplots(2,3,figsize=(15,5),sharex=True)

df = CV.get_from_merge(d['ethane'])
df = df.loc['2017'].resample( 'D' ).mean()
#ax1.plot( df0[0].ceds_only, color=colors[0], label='Base' )
ax1.plot( df0[0]["geo_only"], color=colors[1], label='Base' )
ax1.plot( df0[0]["all_scaled"], color=colors[2], label='Scaled CEDS' )
ax1.plot( df, c='k', label='Observations' )
ax1.legend()
ax1.set_ylabel( 'Ethane (pptv)' )

df = CV.get_from_merge(d['propane'])
df = df.loc['2017'].resample('D').mean()
#ax2.plot( df0[1].ceds_only, c=colors[0], label='Base' )
ax2.plot( df0[1]["geo_only"], color=colors[1], label='Base' )
ax2.plot( df0[1]["all_scaled"], color=colors[2], label='Scaled CEDS' )
ax2.plot( df, c='k', label='Observed' )
ax2.legend()
ax2.set_ylabel( 'Propane (pptv)' )

df = CV.get_from_merge(d['toluene'])
df = df.loc['2017'].resample('D').mean()
#ax5.plot( df0[2].ceds_only, c=colors[0], label='Base' )
ax5.plot( df0[2]["geo_only"], color=colors[1], label='Base' )
ax5.plot( df0[2]["all_scaled"], color=colors[2], label='Scaled CEDS' )
ax5.plot( df, c='k', label='Observed' )
ax5.legend()
ax5.set_ylabel( 'Toluene (pptv)' )

df = CV.get_from_merge(d['benzene'])
df = df.loc['2017'].resample('D').mean()
#ax4.plot( df0[3].ceds_only, c=colors[0], label='Base' )
ax4.plot( df0[3]["geo_only"], color=colors[1], label='Base' )
ax4.plot( df0[3]["all_scaled"], color=colors[2], label='Scaled CEDS' )
ax4.plot( df, c='k', label='Observed' )
ax4.legend()
ax4.set_ylabel( 'Benzene (pptv)' )

df = CV.get_from_merge(d['PRPE'])
df = df.loc['2017'].resample('D').mean()
#ax6.plot( df0[4].ceds_only, c=colors[0], label='Base' )
ax6.plot( df0[4]["geo_only"], color=colors[1], label='Base' )
ax6.plot( df0[4]["all_scaled"], color=colors[2], label='Scaled CEDS' )
ax6.plot( df, c='k', label='Observed' )
ax6.legend()
ax6.set_ylabel( 'Propene (pptv)' )

iB = CV.get_from_merge(d['iso_butane']).fillna(0)
nB = CV.get_from_merge(d['n_butane']).fillna(0)
iP = CV.get_from_merge(d['iso_pentane']).fillna(0)
nP = CV.get_from_merge(d['n_pentane']).fillna(0)
df = iB + nB + iP + nP
df = df.loc['2017'].resample('D').mean()
#ax3.plot( df0[5].ceds_only, c=colors[0], label='Base' )
ax3.plot( df0[5]["geo_only"], color=colors[1], label='Base' )
ax3.plot( df0[5]["all_scaled"], color=colors[2], label='Scaled CEDS' )
ax3.plot( df, c='k', label='Observed' )
ax3.legend()
ax3.set_ylabel( 'ALK4 (pptv)' )

# xticks
ax4.xaxis.set_major_locator(mdates.MonthLocator(interval=3))
ax5.xaxis.set_major_locator(mdates.MonthLocator(interval=3))
ax6.xaxis.set_major_locator(mdates.MonthLocator(interval=3))

plt.tight_layout()
plt.savefig( 'plots/figure_6.png')   
plt.close()

