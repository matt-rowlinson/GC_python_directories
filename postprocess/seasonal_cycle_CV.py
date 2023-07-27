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
                v=var[t,0,y,x]
                var_time.append(v)
            gc=pd.DataFrame({'Value':var_time}, index=time)

            if delta=='M':
                gc=gc.resample(delta).mean()

            ## Get Merge observations
            df=CV.get_from_merge(d[variable])
            df=df[gc.index[0]:gc.index[-1]]
            df=df.resample(delta).mean()
            ## Get mean seasonal cycle
            cv_meanseas = df.groupby(df.index.month).mean()
            gc_meanseas = gc.groupby(gc.index.month).mean()

            ## Get cycle in each year and plot
            f,ax= plt.subplots(figsize=(12,4))
            x=['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
            years=np.unique(df.index.year)
            annmean=[]
            for year in years:
                seas=gc[str(year)]
                ax.plot(x,seas,color='g',alpha=.3,linestyle='--')
                seas=df[str(year)]
                ax.plot(x,seas,color='grey',alpha=.3,linestyle='--')

            
            ax.plot(x, cv_meanseas, 'k', label='CVAO')
            ax.plot(x, gc_meanseas, 'darkgreen', label='GEOS-Chem')
            plt.ylabel('%s (%s)' % (d[variable]['abbr'], d[variable]['unit']) )

            #if delta=='H' or delta=='D':
            #    plt.gca().xaxis.set_major_locator(mdates.DayLocator(interval=interval))
            #if delta=='M' or delta=='Y':
            #    plt.gca().xaxis.set_major_locator(mdates.MonthLocator(interval=interval))

            plt.legend()
            plt.savefig('/users/mjr583/scratch/GC/%s/%s/plots/all_years_seas_%s.png' % (version, rundir, variable) )
            plt.close()

            gc_per75 = gc.groupby(gc.index.month).quantile(.9)
            gc_per25 = gc.groupby(gc.index.month).quantile(.1)

            cv_per75 = df.groupby(df.index.month).quantile(.9)
            cv_per25 = df.groupby(df.index.month).quantile(.1)
            
            f,ax= plt.subplots(figsize=(12,4))
            ax.plot(x, cv_meanseas, 'k', label='CVAO')
            ax.fill_between(x, cv_per25.Value, cv_per75.Value, color='grey', alpha=.2)

            ax.plot(x, gc_meanseas, 'darkgreen', label='GEOS-Chem')
            ax.fill_between(x, gc_per25.Value, gc_per75.Value, color='g', alpha=.2)

            plt.ylabel('%s (%s)' % (d[variable]['abbr'], d[variable]['unit']) )

            plt.legend()
            plt.savefig('/users/mjr583/scratch/GC/%s/%s/plots/seas_with_perc_%s.png' % (version, rundir, variable) )


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
        print(t)
        v=var[t,0,y,x]
        var_time.append(v)
    gc=pd.DataFrame({'Value':var_time}, index=time)

    if delta=='M':
        gc=gc.resample(delta).mean()

    ## Get Merge observations
    df=CV.get_from_merge(d[variable])
    df=df[gc.index[0]:gc.index[-1]]
    df=df.resample(delta).mean()
    
    df=df['2015':]
    gc=gc['2015':]
    ## Get mean seasonal cycle
    cv_meanseas = df.groupby(df.index.month).mean()
    gc_meanseas = gc.groupby(gc.index.month).mean()

    ## Get cycle in each year and plot
    f,ax= plt.subplots(figsize=(12,4))
    x=['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
    years=np.unique(df.index.year)
    annmean=[]
    for year in years:
        seas=gc[str(year)]
        ax.plot(x,seas,color='g',alpha=.3,linestyle='--')
        seas=df[str(year)]
        ax.plot(x,seas,color='grey',alpha=.3,linestyle='--')


    ax.plot(x, cv_meanseas, 'k', label='CVAO')
    ax.plot(x, gc_meanseas, 'darkgreen', label='GEOS-Chem')
    plt.ylabel('%s (%s)' % (d[variable]['abbr'], d[variable]['unit']) )

    plt.legend()
    plt.savefig('/users/mjr583/scratch/GC/%s/%s/plots/all_years_seas_%s.png' % (version, rundir, variable) )
    plt.close()

    gc_per75 = gc.groupby(gc.index.month).quantile(.9)
    gc_per25 = gc.groupby(gc.index.month).quantile(.1)

    cv_per75 = df.groupby(df.index.month).quantile(.9)
    cv_per25 = df.groupby(df.index.month).quantile(.1)
    
    f,ax= plt.subplots(figsize=(12,4))
    ax.plot(x, cv_meanseas, 'k', label='CVAO')
    ax.fill_between(x, cv_per25.Value, cv_per75.Value, color='grey', alpha=.2)

    ax.plot(x, gc_meanseas, 'darkgreen', label='GEOS-Chem')
    ax.fill_between(x, gc_per25.Value, gc_per75.Value, color='g', alpha=.2)

    plt.ylabel('%s (%s)' % (d[variable]['abbr'], d[variable]['unit']) )

    plt.legend()
    plt.savefig('/users/mjr583/scratch/GC/%s/%s/plots/seas_with_perc_%s.png' % (version, rundir, variable) )
