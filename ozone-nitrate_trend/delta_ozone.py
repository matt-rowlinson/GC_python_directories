#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#SBATCH --job-name=NOxtimeseries
#SBATCH --ntasks=1
#SBATCH --mem=18Gb
#SBATCH --partition=nodes
#SBATCH --time=00:45:00
#SBATCH --output=Logs/timeseries_%A.log
import sys
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
matplotlib.use('agg')
sys.path.append('/users/mjr583/python_lib')
import GC_tools as GC
import RowPy as rp
plt.style.use('seaborn-darkgrid')

def get_model_as_dataframe(rundirs, variable, version="13.1.2", collection="SpeciesConc", year="2017", x=-24.9, y=16.9, z=0):  
    dfs=[]
    for a in rundirs:
        print(a, year)
        var, lat, lon, lev, atime = GC.get_gc_var(rundir=a, variable=variable, collection=collection, version='13.1.2', year=year)
        a_time=[] 
        idy = rp.find_nearest(lat, y)
        idx = rp.find_nearest(lon, x)
        for t in range(len(atime)):
            v=var[t,z,idy,idx]
            a_time.append(v)
        a=pd.DataFrame({a:a_time}, index=atime)
        dfs.append(a)
    df=pd.concat(dfs, axis=1, ignore_index=False)
    print(df)
    return df

def dataframe_by_year(year):
    rundirs=[f'nitrate_photol_{year}_control',f'nitrate_photol_{year}_scale-25',
             f'nitrate_photol_{year}_scale-50',f'nitrate_photol_{year}_scale-100']
    df = get_model_as_dataframe( rundirs, variable="O3", year=year )
    df.columns=['control','j25','j50','j100']
    df=df.resample('M').mean()
    return df

def plot_by_scale(ax, df, label, c='r', first=False):
    idx, z, p, X = trendline( df )
    if first:
        ax.plot( df.index, df, c,label=f'{label}  (y=%.2f ppb decade$^{-1}$)' %z)
    else:
        ax.plot( df.index, df, c,label=f'{label}  (y=%.2f)' %z)
    ax.plot( df.index[idx], p(X), '--', c=c)
    return

def main():
    #df_1980 = dataframe_by_year('1980')
    #df_2000 = dataframe_by_year('2000')
    #df_2005 = dataframe_by_year('2005')
    df_2010 = dataframe_by_year('2010')
    df_2017 = dataframe_by_year('2017')
    #print( df_1980.count)
    #print( df_2000.count)
    #print( df_2005.count)
    print( df_2010.count)
    print( df_2017.count)
    
    js=['j25','j50','j100']
    for j in js:
        f,ax=plt.subplots(figsize=(10.5,3.75))
        ax.plot( df_1980.index, df_1980[j] , label='1980' )
        ax.plot( df_2000.index, df_2000[j] , label='2000' )
        ax.plot( df_2005.index, df_2005[j] , label='2005' )
        ax.plot( df_2010.index, df_2010[j] , label='2010' )
        ax.plot( df_2017.index, df_2017[j] , label='2017' )

        ax.set_ylabel( '$O_3$ (ppbv)')  
        plt.suptitle('$O_3$ at J$_{%sscale}$' %j[1:])
        plt.legend()
        plt.savefig(f'plots/o3_by_year_{j}.png')
        plt.close()

    for j in js:
        f,ax=plt.subplots(figsize=(10.5,3.75))
        ax.plot( df_1980.index, df_1980[j] - df_1980['control'], label='1980' )
        ax.plot( df_2000.index, df_2000[j] - df_2000['control'], label='2000' )
        ax.plot( df_2005.index, df_2005[j] - df_2005['control'], label='2005' )
        ax.plot( df_2010.index, df_2010[j] - df_2010['control'], label='2010' )
        ax.plot( df_2017.index, df_2017[j] - df_2017['control'], label='2017' )

        ax.set_ylabel( '$O_3$ (ppbv)')  
        plt.suptitle('Delta $O_3$ at J$_{%sscale}$' %j[1:])
        plt.legend()
        plt.savefig(f'plots/delta_o3_over_time_{j}.png')
        plt.close()

    js=['j25','j50','j100']
    for j in js:
        f,ax=plt.subplots(figsize=(10.5,3.75))
        ax.plot( df_1980.index, df_1980[j] - df_1980['control'], label='1980' )
        ax.plot( df_1980.index, df_2000[j] - df_2000['control'], label='2000' )
        ax.plot( df_1980.index, df_2005[j] - df_2005['control'], label='2005' )
        ax.plot( df_1980.index, df_2010[j] - df_2010['control'], label='2010' )
        ax.plot( df_1980.index, df_2017[j] - df_2017['control'], label='2017' )

        import matplotlib.dates as mdates
        myFmt = mdates.DateFormatter('%m')
        ax.xaxis.set_major_formatter(myFmt)

        ax.set_ylabel( 'delta $O_3$ (ppbv)')  
        plt.suptitle('Delta $O_3$ at J$_{%sscale}$' %j[1:])
        plt.legend()
        plt.savefig(f'plots/delta_o3_by_year_{j}.png')
        plt.close()

    for j in js:
        f,ax=plt.subplots(figsize=(10.5,3.75))
        ax.plot( df_1980.index, df_1980[j] / df_1980['control'] * 100, label='1980' )
        ax.plot( df_1980.index, df_2000[j] / df_2000['control'] * 100, label='2000' )
        ax.plot( df_1980.index, df_2005[j] / df_2005['control'] * 100, label='2005' )
        ax.plot( df_1980.index, df_2010[j] / df_2010['control'] * 100, label='2010' )
        ax.plot( df_1980.index, df_2017[j] / df_2017['control'] * 100, label='2017' )

        import matplotlib.dates as mdates
        myFmt = mdates.DateFormatter('%m')
        ax.xaxis.set_major_formatter(myFmt)

        ax.set_ylabel( 'j / control * 100')  
        plt.suptitle('Delta $O_3$ at J$_{%sscale}$' %j[1:])
        plt.legend()
        plt.savefig(f'plots/ratio_o3_by_year_{j}.png')
        plt.close()

if __name__=="__main__":
    main()
