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

#rundirs=['control','ceds_only', 'nei_scaling','emep_scaling','nei_emep_scaling', 'inc_asia_scaling']
rundirs=['ceds_only','new_scale_all_vocs']#ceds_only', 'asia_scale_2.5','global_scale_2.22']
df0=[] ; cvs=[]
variables=['ethane','propane']
for variable in variables:
    dfs=[]
    for a in rundirs:
        print(a)
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
    labels=['Base','Scaled_CEDS','No VOC emissions']#NEI/EMEP/Asia scale','Global scale']#,'CEDS + EMEP scaling','CEDS + NEI scaling','CEDS + NEI/EMEP scaling', 'Inc. Asia scale']
    colors=['#1b9e77','darkorange','#7570b3']
    
    ## Turn into single plot with subplots
    df0.append( cv_df )

df = CV.get_from_merge(d['ethane'])
df = df['2017'].resample( 'D' ).mean()

f, (ax1,ax2) = plt.subplots(2,1,figsize=(7,6),sharex=True)

ax1.plot( df0[0].ceds_only, color=colors[0], label='Base (CEDS)' )
ax1.plot( df0[0].new_scale_all_vocs, color=colors[1], label='Scaled CEDS' )
ax1.plot( df, c='k', label='Observations' )
ax1.legend()
ax1.set_ylabel( 'C$_2$H$_6$ (pptv)' )

df = CV.get_from_merge(d['propane'])
df = df['2017'].resample('D').mean()
ax2.plot( df0[1].ceds_only, c=colors[0], label='Base (CEDS)' )
ax2.plot( df0[1].new_scale_all_vocs, c=colors[1], label='Scaled CEDS' )
ax2.plot( df, c='k', label='Observed' )
ax2.legend()
ax2.set_ylabel( 'C$_3$H$_8$ (pptv)' )

import matplotlib.dates as mdates
myFmt = mdates.DateFormatter('%b')
ax2.xaxis.set_major_formatter(myFmt)
plt.gca().xaxis.set_major_locator(mdates.MonthLocator(interval=2))


plt.tight_layout()
plt.savefig( 'new_plots/figure_7.png')    
plt.close()
