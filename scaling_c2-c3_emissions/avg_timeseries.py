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
rundirs=['ceds_only','scale_all_vocs']#ceds_only', 'asia_scale_2.5','global_scale_2.22']

string='NEWTALK_'
variables=['ALK4','toluene','PRPE','propane','ethane','benzene','O3','CO','OH']
#variables=['XYLE','MEK','CH2O']
#variables=['O3','CO','ethane','propane','benzene','toluene','eoh','propene', 'ALK4']
variables=['toluene']
for variable in variables:
    dfs=[]
    for a in rundirs:
        print(a)
        var, lat, lon, lev, atime = GC.get_gc_var(rundir=a, variable=variable, version='13.1.2', year='2017')
        
        a_time=[] 
        y = rp.find_nearest(lat, 16.9)
        x = rp.find_nearest(lon, -24.9)
        print( lat[y], lat[y+1] )
        print( lon[x], lon[x+1] )
        for t in range(len(atime)):
            v=var[t,0,y,x]
            a_time.append(v)
        a=pd.DataFrame({a:a_time}, index=atime)
        dfs.append(a)
    
    cv_df=pd.concat(dfs, axis=1, ignore_index=False)
    print( cv_df )
    dfs[0]=dfs[0].resample('D').mean()
    labels=['Base','Scaled_CEDS','No VOC emissions']#NEI/EMEP/Asia scale','Global scale']#,'CEDS + EMEP scaling','CEDS + NEI scaling','CEDS + NEI/EMEP scaling', 'Inc. Asia scale']
    styles=['solid','dashed','dashdot','dotted','dashed','dashdot']
    colors=['#1b9e77','#d95f02','#7570b3']

    df = CV.get_from_merge(d[variable])
    df = df.groupby( df.index.dayofyear ).mean()
    df = df.drop(df.index[29])
    df.index = dfs[0].index
    
    f,ax = plt.subplots(figsize=(12,4))
    df=df.resample('D').mean()
    ax.plot(df.index, df.Value, 'k', label='CVAO')

    for n,df in enumerate(dfs):
        #df=df.resample(resample).mean()
        ax.plot(df.index, df, color=colors[n],linestyle=styles[n],label=labels[n])
    #df=df.resample('D').mean()
    #ax.plot(df.index, df.Value, 'k', label='CVAO')
    plt.ylabel(f'{variable} (pptv)')
    plt.tight_layout()
    plt.savefig(f'plots/{string}timeseries_v13_{variable}.png')
    plt.close()
    
    #rp.timeseries_plot(dfs, variable, cvao_data=True,  string=string, labels=labels, styles=styles, colors=colors)
    #rp.timeseries_plot(dfs, variable, cvao_data=True, std=True, resample='D', string=string, labels=labels, styles=styles, colors=colors)
    #plt.close()
