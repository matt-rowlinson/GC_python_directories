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

rundirs=['geo_2x25','all_2x25']
df0=[] ; cvs=[]
variables=['ethane','propane','TOLU','BENZ','PRPE', 'ALK4']
scales=[1e12,1e12,1e9,1e12,1e12,1e12]
for variable in variables:
    print( variable )
    dfs=[]
    for a in rundirs:
        var, lat, lon, lev, atime = GC.get_gc_var(rundir=a, variable=variable, version='14.0.1', year='2017')
        
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
    colors=['forestgreen','darkorange']
    
    ## Turn into single plot with subplots
    df0.append( cv_df )


f, ((ax1,ax2),(ax3,ax4),(ax5,ax6)) = plt.subplots(3,2,figsize=(14,10),sharex=True)

df = CV.get_from_merge(d['ethane'])
df = df.loc['2017'].resample( 'D' ).mean()
ax1.plot( df0[0][rundirs[0]], color=colors[0], label='Base' )
ax1.plot( df0[0][rundirs[1]], color=colors[1], label='Scaled CEDS emissions' )
ax1.plot( df, c='k', label='Observations' )
ax1.legend()
ax1.set_ylabel( 'C2H6 (pptv)' )

df = CV.get_from_merge(d['propane'])
df = df.loc['2017'].resample('D').mean()
ax2.plot( df0[1][rundirs[0]], c=colors[0], label='Base' )
ax2.plot( df0[1][rundirs[1]], color=colors[1], label='Scaled CEDS emissions' )
ax2.plot( df, c='k', label='Observed' )
ax2.legend()
ax2.set_ylabel( 'C3H8 (pptv)' )

df = CV.get_from_merge(d['TOLU'])
df = df.loc['2017'].resample('D').mean()
ax3.plot( df0[2][rundirs[0]], c=colors[0], label='Base' )
ax3.plot( df0[2][rundirs[1]], color=colors[1], label='Scaled CEDS emissions' )
ax3.plot( df, c='k', label='Observed' )
ax3.legend()
ax3.set_ylabel( 'Toluene (pptv)' )

df = CV.get_from_merge(d['BENZ'])
df = df.loc['2017'].resample('D').mean()
ax4.plot( df0[3][rundirs[0]], c=colors[0], label='Base' )
ax4.plot( df0[3][rundirs[1]], color=colors[1], label='Scaled CEDS emissions' )
ax4.plot( df, c='k', label='Observed' )
ax4.legend()
ax4.set_ylabel( 'Benzene (pptv)' )

df = CV.get_from_merge(d['PRPE'])
df = df.loc['2017'].resample('D').mean()
ax6.plot( df0[4][rundirs[0]], c=colors[0], label='Base' )
ax6.plot( df0[4][rundirs[1]], color=colors[1], label='Scaled CEDS emissions' )
ax6.plot( df, c='k', label='Observed' )
ax6.legend()
ax6.set_ylabel( 'Propene (pptv)' )

iB = CV.get_from_merge(d['iso_butane']).fillna(0)
nB = CV.get_from_merge(d['n_butane']).fillna(0)
iP = CV.get_from_merge(d['iso_pentane']).fillna(0)
nP = CV.get_from_merge(d['n_pentane']).fillna(0)
df = iB + nB + iP + nP
df = df.loc['2017'].resample('D').mean()
ax5.plot( df0[5][rundirs[0]], c=colors[0], label='Base' )
ax5.plot( df0[5][rundirs[1]], color=colors[1], label='Scaled CEDS emissions' )
ax5.plot( df, c='k', label='Observed' )
ax5.legend()
ax5.set_ylabel( 'ALK4 (pptv)' )

plt.delaxes(ax6)
plt.tight_layout()
plt.savefig( 'plots/figure_5.png')   
plt.close()

