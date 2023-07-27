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

#rundir=sys.argv[1]
#species=sys.argv[2]
#version='12.9.3'
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

            plt.savefig('/users/mjr583/scratch/GC/%s/%s/plots/scatter_GC-CV_%s.png' % (version, rundir, variable) )


        except:
            print(variable, 'failed to plot')
            continue


else:
    var, lat, lon, lev, time = GC.get_gc_var(rundir, 'O3', version)
    #ms=GC.closest_met(time[jobid])

    #t=GC.hour_rounder(time[jobid])
    #diff=GC.dt_difference(t,ms)
    
    GC_T=[]
    for t in time:
        try:
            T,mlat,mlon = GC.get_gc_input(t, var='T2M',filetype='A1')
        except:
            continue
        ladx = rp.find_nearest(mlat, 16.9)
        lodx = rp.find_nearest(mlon, -24.9)
        T=T[ladx,lodx] - 273.15
        GC_T.append(T)
    
    met=pd.DataFrame({'T2M':GC_T}, index=time)
    print(met)

    stop

    stop 


    delta,interval=GC.find_timestep(time)
    y = rp.find_nearest(lat, 16.9)
    x = rp.find_nearest(lon, -24.9)
    var_time=[]
    for t in range(len(time)):
        print(t)
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

    plt.savefig('/users/mjr583/scratch/GC/%s/%s/plots/scatter_GC-CV_%s.png' % (version, rundir, variable) )

    '''
    f,ax= plt.subplots(figsize=(12,4))
    ax.plot(df.index, df.Value, 'k', label='CVAO')
    ax.plot(gc.index, gc.Value, 'g', label='GEOS-Chem')
    plt.ylabel('%s (%s)' % (d[variable]['abbr'], d[variable]['unit']) )

    if delta=='H' or delta=='D':
        plt.gca().xaxis.set_major_locator(mdates.DayLocator(interval=interval))
    if delta=='M' or delta=='Y':
        plt.gca().xaxis.set_major_locator(mdates.MonthLocator(interval=interval))

    plt.legend()
    plt.savefig('/users/mjr583/scratch/GC/%s/%s/plots/timeseries_%s.png' % (version, rundir, variable) )
    '''
