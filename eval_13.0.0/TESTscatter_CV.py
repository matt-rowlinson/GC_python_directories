#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#SBATCH --job-name=timeseries
#SBATCH --ntasks=1
#SBATCH --mem=10gb
#SBATCH --partition=interactive
#SBATCH --time=00:10:00
#SBATCH --output=LOGS/timeseries.log
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.lines as mlines
import matplotlib.transforms as mtransforms
import matplotlib
matplotlib.use('agg')
sys.path.append('/users/mjr583/python_lib')
import GC_tools as GC
import RowPy as rp
from CVAO_dict import CVAO_dict as d
import CVAO_tools as CV
plt.style.use('seaborn-darkgrid')

inputs=GC.get_arguments()
rundir=inputs.rundir
variable=inputs.var
version=inputs.version

if variable=='all':
    for variable in d:
        print(variable)
        try:
            var, lat, lon, lev, time = GC.get_gc_var(rundir, variable, version)

            delta,interval=GC.find_timestep(time)
            y = rp.find_nearest(lat, 16.9)
            x = rp.find_nearest(lon, -24.9)
            var_time=[]
            for t in range(len(time)):
                #print(t)
                v=var[t,0,y,x]
                #print(v)
                var_time.append(v)
            gc=pd.DataFrame({'Value':var_time}, index=time)

            if delta=='M':
                gc=gc.resample(delta).mean()

            ## Get Merge observations
            df=CV.get_from_merge(d[variable])
            df=df[gc.index[0]:gc.index[-1]]
            df=df.resample(delta).mean()

            f,ax= plt.subplots(figsize=(8,8))
            ax.scatter(df.Value,gc.Value)
            plt.ylabel('GEOS-Chem %s (%s)' % (d[variable]['abbr'], d[variable]['unit']) )
            plt.xlabel('CVAO %s (%s)' % (d[variable]['abbr'], d[variable]['unit']) )

            Max=np.max([gc.Value.max(), df.Value.max()])
            Min=np.min([gc.Value.min(), df.Value.min()])
            plt.xlim(Min,Max)
            plt.ylim(Min,Max)

            line = mlines.Line2D([0, 1], [0, 1], color='red', alpha=0.2)
            transform = ax.transAxes
            line.set_transform(transform)
            ax.add_line(line)#, alpha=0.6)

            absError= gc.Value - df.Value
            SE = np.square(absError)
            MSE= np.mean(SE)
            RMSE=np.round( np.sqrt(MSE), 3)
            R2 = np.round( 1. - (np.var(absError) / np.var(df.Value)), 3)

            plt.text(.1,.9, 'RMSE: %s' %RMSE, fontsize=20,horizontalalignment='left',
                         verticalalignment='center',
                          transform = ax.transAxes)
            plt.text(.1,.8, '$R^2$: %s' %R2, fontsize=20,horizontalalignment='left',
                         verticalalignment='center',
                          transform = ax.transAxes)

            plt.savefig('/users/mjr583/scratch/GC/%s/%s/plots/scatter_GC-CV_%s.png' % (version, rundir, variable) )


        except:
            print(variable, 'failed to plot')
            continue


else:
    var, lat, lon, lev, time = GC.get_gc_var(rundir, variable, version)

    delta,interval=GC.find_timestep(time)
    y = rp.find_nearest(lat, 16.9)
    x = rp.find_nearest(lon, -24.9)
    var_time=[]
    for t in range(len(time)):
        v=var[t,0,y,x]
        #print(v)
        var_time.append(v)
    gc=pd.DataFrame({'Value':var_time}, index=time)

    if delta=='M':
        gc=gc.resample(delta).mean()

    ## Get Merge observations
    df=CV.get_from_merge(d[variable])
    df=df[gc.index[0]:gc.index[-1]]
    df=df.resample(delta).mean()

    ## Split seasonally here
    ggc=gc
    ddf=df
    djf=[1,2,12] ; mam=[3,4,5] ; jja=[6,7,8] ; son=[9,10,11]
    labs=['DJF','MAM','JJA','SON']
    seas=[djf, mam, jja, son]
    for n in range(4):
        months=seas[n]
        df=ddf[ddf.index.month.isin(months)]
        gc=ggc[ggc.index.month.isin(months)]
    
        f,ax= plt.subplots(figsize=(8,8))
        ax.scatter(df.Value,gc.Value)
        plt.ylabel('GEOS-Chem %s (%s)' % (d[variable]['abbr'], d[variable]['unit']) )
        plt.xlabel('CVAO %s (%s)' % (d[variable]['abbr'], d[variable]['unit']) )

        Max=np.max([gc.Value.max(), df.Value.max()])
        Min=np.min([gc.Value.min(), df.Value.min()])
        plt.xlim(Min,Max)
        plt.ylim(Min,Max)

        line = mlines.Line2D([0, 1], [0, 1], color='red', alpha=0.2)
        transform = ax.transAxes
        line.set_transform(transform)
        ax.add_line(line)#, alpha=0.6)
        
        df.index=gc.index
        absError= gc.Value - df.Value
        
        SE = np.square(absError)
        MSE= np.nanmean(SE)
        RMSE=np.round( np.sqrt(MSE), 3)
        R2 = np.round( 1. - (np.var(absError) / np.var(df.Value)), 3)

        plt.text(.1,.9, 'RMSE: %s' %RMSE, fontsize=20,horizontalalignment='left',
                     verticalalignment='center',
                          transform = ax.transAxes)
        plt.text(.1,.8, '$R^2$: %s' %R2, fontsize=20,horizontalalignment='left',
                     verticalalignment='center',
                          transform = ax.transAxes)
        plt.title(labs[n], fontsize=22)

        plt.savefig('./plots/scatter_%s_%s_%s.png' % (version, variable, labs[n]) )

ddf.index=ggc.index

x=len(ddf)
ten=int(x/10)

df=ddf.Value.nlargest(ten).sort_index()#.reset_index()
gc=ggc.reindex(df.index)
df=pd.DataFrame({'Value':df}, index=gc.index)
gc.columns=['Value']

f,ax= plt.subplots(figsize=(8,8))
ax.scatter(df.Value,gc.Value)
plt.ylabel('GEOS-Chem %s (%s)' % (d[variable]['abbr'], d[variable]['unit']) )
plt.xlabel('CVAO %s (%s)' % (d[variable]['abbr'], d[variable]['unit']) )

Max=np.max([gc.Value.max(), df.Value.max()])
Min=np.min([gc.Value.min(), df.Value.min()])
plt.xlim(Min,Max)
plt.ylim(Min,Max)

line = mlines.Line2D([0, 1], [0, 1], color='red', alpha=0.2)
transform = ax.transAxes
line.set_transform(transform)
ax.add_line(line)#, alpha=0.6)

df.index=gc.index
absError= gc.Value - df.Value

SE = np.square(absError)
MSE= np.nanmean(SE)
RMSE=np.round( np.sqrt(MSE), 3)
R2 = np.round( 1. - (np.var(absError) / np.var(df.Value)), 3)

plt.text(.1,.9, 'RMSE: %s' %RMSE, fontsize=20,horizontalalignment='left',
             verticalalignment='center',
                  transform = ax.transAxes)
plt.text(.1,.8, '$R^2$: %s' %R2, fontsize=20,horizontalalignment='left',
             verticalalignment='center',
                  transform = ax.transAxes)
plt.title('Highest 10%', fontsize=22)
plt.savefig('./plots/scatter_%s_%s_%s.png' % (version, variable, 'top10') )


## Smallest 10 %
df=ddf.Value.nsmallest(ten).sort_index()
gc=ggc.reindex(df.index)
df=pd.DataFrame({'Value':df}, index=gc.index)
gc.columns=['Value']

f,ax= plt.subplots(figsize=(8,8))
ax.scatter(df.Value,gc.Value)
plt.ylabel('GEOS-Chem %s (%s)' % (d[variable]['abbr'], d[variable]['unit']) )
plt.xlabel('CVAO %s (%s)' % (d[variable]['abbr'], d[variable]['unit']) )

Max=np.max([gc.Value.max(), df.Value.max()])
Min=np.min([gc.Value.min(), df.Value.min()])
plt.xlim(Min,Max)
plt.ylim(Min,Max)

line = mlines.Line2D([0, 1], [0, 1], color='red', alpha=0.2)
transform = ax.transAxes
line.set_transform(transform)
ax.add_line(line)#, alpha=0.6)

df.index=gc.index
absError= gc.Value - df.Value

SE = np.square(absError)
MSE= np.nanmean(SE)
RMSE=np.round( np.sqrt(MSE), 3)
R2 = np.round( 1. - (np.var(absError) / np.var(df.Value)), 3)

plt.text(.1,.9, 'RMSE: %s' %RMSE, fontsize=20,horizontalalignment='left',
             verticalalignment='center',
                  transform = ax.transAxes)
plt.text(.1,.8, '$R^2$: %s' %R2, fontsize=20,horizontalalignment='left',
             verticalalignment='center',
                  transform = ax.transAxes)
plt.title('Lowest 10%', fontsize=22)

plt.savefig('./plots/scatter_%s_%s_%s.png' % (version, variable, 'bottom10') )

