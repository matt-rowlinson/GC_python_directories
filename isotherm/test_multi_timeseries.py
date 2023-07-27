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
        print(a)
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
    return df

def trendline(df):
    XX=np.arange(len(df))
    idx=np.isfinite(df)
    Y=df[idx]
    X=XX[idx]
    time=df.index[idx]

    z = np.polyfit(X, Y, 1)
    p = np.poly1d(z)
    z = z[0] * 12. * 10.
    return idx, z, p, X

def dataframe_by_year(year):
    rundirs=[f'nitrate_photol_{year}_control',f'base_run_{year}',f'j100_run_{year}']
    try:
        df = get_model_as_dataframe( rundirs, variable="O3", year=year )
    except:
        df = get_model_as_dataframe( rundirs, variable="O3", year='' )
    df.columns=['No_dust','Base','j100']
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
    df_1980 = dataframe_by_year('1980')
    df_1990 = dataframe_by_year('1990')
    df_2000 = dataframe_by_year('2000')
    df_2005 = dataframe_by_year('2005')
    df_2010 = dataframe_by_year('2010')
    df_2017 = dataframe_by_year('2017')

    merged = pd.concat( [df_1980, df_2000, df_2005, df_2010, df_2017 ] )
    merged=merged.resample('M').mean()
    
    f,ax=plt.subplots(figsize=(14,5))
    colors=['#1b9e77','#d95f02','#7570b3','#e7298a']
    plot_by_scale(ax, merged['control'], label='Control', c=colors[0], first=True)
    plot_by_scale(ax, merged['j25'],    label=r'NITs J$_{25scale}$',  c=colors[1])
    plot_by_scale(ax, merged['j50'],     label=r'J$_{50scale}$',  c=colors[2])
    plot_by_scale(ax, merged['j100'],    label=r'J$_{100scale}$', c=colors[3])
    
    ax.set_ylabel( r'$O_3$ (ppbv)')    
    plt.legend()
    plt.savefig('plots/trends_by_Jscaling.png')
    plt.close()

if __name__=="__main__":
    main()
