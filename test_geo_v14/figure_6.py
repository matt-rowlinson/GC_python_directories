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

rundirs=['geo','scale_all','geo_2x25','all_2x25']
df0=[] ; cvs=[]
variables=['ethane','propane','TOLU','BENZ','PRPE', 'ALK4']
scales=[1e12,1e12,1e12,1e12,1e12,1e12]
for variable in variables:
    print( variable )
    dfs=[]
    for a in rundirs:
        var, lat, lon, lev, atime = GC.get_gc_var(rundir=a, variable=variable, version='14.0.1', year='2017')
        
        a_time=[] 
        y = rp.find_nearest(lat, 16.9)
        x = rp.find_nearest(lon, -24.9)
        for t in range(len(atime)):
            v=var[t,0,y+1,x-1] #* 1e9
            a_time.append(v)
        a=pd.DataFrame({a:a_time}, index=atime)
        dfs.append(a)
    
    cv_df=pd.concat(dfs, axis=1, ignore_index=False)
    dfs[0]=dfs[0].resample('D').mean()
    
    ## Turn into single plot with subplots
    df0.append( cv_df )

labels=['Base 4x5','Scaled 4x5','Base 2x2.5','Scaled 2x2.5']
colors=['blue','red','forestgreen','darkorange']
style =['-','-','--','--']
ylabel=['C2H6','C3H8','Toluene','Benzene','Propene'] 
f, ((ax1,ax2),(ax3,ax4),(ax5,ax6)) = plt.subplots(3,2,figsize=(14,10),sharex=True)
axes=[ax1,ax2,ax3,ax4,ax5]
for n, v in enumerate(variables[:-1]):
    df = CV.get_from_merge(d[v])
    df = df.loc['2017'].resample( 'D' ).mean()
    for nn, run in enumerate(rundirs):
        axes[n].plot( df0[n][run]*scales[n], color=colors[nn], linestyle=style[nn], label=labels[nn] )
    axes[n].plot( df, c='k', label='Observations' )
    axes[n].legend()
    axes[n].set_ylabel( ylabel[n]+' (pptv)' )

iB = CV.get_from_merge(d['iso_butane']).fillna(0)
nB = CV.get_from_merge(d['n_butane']).fillna(0)
iP = CV.get_from_merge(d['iso_pentane']).fillna(0)
nP = CV.get_from_merge(d['n_pentane']).fillna(0)
df = iB + nB + iP + nP
df = df.loc['2017'].resample('D').mean()
#ax6.plot( df0[5][rundirs[0]], c=colors[0], label='Base' )

ax6.plot( df0[0][rundirs[0]]*1e12, color=colors[0], linestyle=style[0], label=labels[0] )
ax6.plot( df0[0][rundirs[1]]*1e12, color=colors[1], linestyle=style[1], label=labels[1] )
ax6.plot( df0[0][rundirs[2]]*1e12, color=colors[2], linestyle=style[2], label=labels[2] )
ax6.plot( df0[0][rundirs[3]]*1e12, color=colors[3], linestyle=style[3], label=labels[3] )

ax6.plot( df, c='k', label='Observed' )
ax6.legend()
ax6.set_ylabel( 'ALK4 (pptv)' )

plt.tight_layout()
plt.savefig( 'plots/2x25.figure_6.png')   
plt.close()

