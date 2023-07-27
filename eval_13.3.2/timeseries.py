#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#SBATCH --job-name=NOxtimeseries
#SBATCH --ntasks=1
#SBATCH --mem=2Gb
#SBATCH --partition=interactive
#SBATCH --time=00:15:00
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

rundirs=['control', 'scale_all_vocs','control','all_scale']
versions=['13.1.2','13.1.2','13.3.2','13.3.2']

string='v13_Scaling_'
variables=['OH']#,'OH','ethane','propane','benzene','toluene']
for n, variable in enumerate(variables):
    dfs=[]
    for n, a in enumerate(rundirs):
        print(a, versions[n], variable)
        var, lat, lon, lev, atime = GC.get_gc_var(rundir=a, variable=variable, version=versions[n], year='2016')
        
        a_time=[] 
        y = rp.find_nearest(lat, 16.9)
        x = rp.find_nearest(lon, -24.9)
        for t in range(len(atime)):
            v=var[t,0,y,x]
            a_time.append(v)
        a=pd.DataFrame({a:a_time}, index=atime)
        dfs.append(a)
    
    labels=['v13.1.3 Control', 'v13.1.2 Scaled','v13.3.2 Control', 'v13.3.2 Scaled']
    styles=['solid','dashed','solid','dashed','dashdot','dotted','dashed','dashdot']
    colors=['g','r','b','y','purple','orange']

    f,ax= plt.subplots(figsize=(12,4))
    for n,df in enumerate(dfs):
        ax.plot(df.index, df, color=colors[n],linestyle=styles[n],label=labels[n])

    #cv=CV.get_from_merge(d[variable])
    #cv=cv[dfs[0].index[0]:dfs[0].index[-1]]
    #cv=cv.resample('D').mean()
    #ax.plot(cv.index, cv.Value, 'k', label='CVAO')

    plt.ylabel( 'ppb' )
    plt.ylabel('%s (%s)' % (d[variable]['abbr'], d[variable]['unit']) )

    plt.gca().xaxis.set_major_locator(mdates.MonthLocator(interval=2))
    plt.legend(ncol=2)
    plt.savefig('./plots/%stimeseries_v13_%s.png' %(string, variable) )
    plt.close()
    print( 'done' )

    f,ax= plt.subplots(figsize=(12,4))
    ax.plot(dfs[1].index, dfs[1].values - dfs[0].values, label="v13.1.2")
    ax.plot(dfs[3].index, dfs[3].values - dfs[2].values, label="v13.3.2")

    #cv=CV.get_from_merge(d[variable])
    #cv=cv[dfs[0].index[0]:dfs[0].index[-1]]
    #cv=cv.resample('D').mean()
    #ax.plot(cv.index, cv.Value, 'k', label='CVAO')

    plt.ylabel(f'delta {d[variable]["abbr"]} {d[variable]["unit"]}' )

    plt.gca().xaxis.set_major_locator(mdates.MonthLocator(interval=2))
    plt.legend()#ncol=2)
    plt.savefig(f'./plots/delta{variable}_timeseries_v13.png')
    plt.close()
    print( 'done' )





