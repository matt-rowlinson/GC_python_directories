#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#SBATCH --job-name=cv-nitrate
#SBATCH --ntasks=1
#SBATCH --mem=18Gb
#SBATCH --partition=interactive
#SBATCH --time=00:45:00
#SBATCH --output=Logs/cv-nitrate_%A.log
import pandas as pd
import sys
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('agg')
sys.path.append('/users/mjr583/python_lib')
import GC_tools as GC
import RowPy as rp
plt.style.use('seaborn-darkgrid')
sys.path.append('/users/mjr583/cvao')
import read_nitrate
import numpy as np

def get_model_as_dataframe(rundirs, variable, version="13.1.2", collection="SpeciesConc", year="2017", x=-24.9, y=16.9, z=0):  
    dfs=[]
    for a in rundirs:
        print(a)
        var, lat, lon, lev, atime = GC.get_gc_var(rundir=a, variable=variable, collection=collection, version='13.1.2', year=year)
        #var = var * 1e12
        a_time=[] 
        idy = rp.find_nearest(lat, y)
        idx = rp.find_nearest(lon, x)
        for t in range(len(atime)):
            v=var[t,z,idy,idx]
            a_time.append(v)
        a=pd.DataFrame({a:a_time}, index=atime)
        dfs.append(a)
    cv_df=pd.concat(dfs, axis=1, ignore_index=False)
    return dfs

def gc_nitrate(rundir, year='20'):
    NITs = get_model_as_dataframe( rundir, variable="NITs", year=year)[0] 
    NIT = get_model_as_dataframe( rundir, variable="NIT", year=year)[0]
    NO3 = (NITs + NIT) #* 1e9 * 2.54 
    NO3.columns = ["NO3 / ppb"]
    NO3=NO3['2007':'2011']
    return NO3

def main():
    collection='SpeciesConc'
    rundir=['nitrate_photol_control']
    year="20"
    
    ds0 = gc_nitrate(['nitrate_photol_control'], year='2010')
    #ds1 = gc_nitrate(['nitrate_photol_all_Jscale-25'], year='20')
    
    print( ds0 )
    sys.exit()
    #NITs = get_model_as_dataframe( rundir, variable="NITs", year=year)[0] 
    #NIT = get_model_as_dataframe( rundir, variable="NIT", year=year)[0]
    #NO3 = (NITs + NIT) * 1e9 
    #NO3.columns = ["NO3 / ppb"] 
    #NO3=NO3['2007':'2011']

    df = read_nitrate.read_nitrate_dataset(convert=False)
    
    f,ax=plt.subplots(1,1, figsize=(10,4))
    ax.plot( ds0.index, ds0, c='g', label='GEOS-Chem v13.1.2')
    ax.plot( df.index,  df, c='k', label='Obs')

    plt.legend()
    plt.ylabel("NO3 (ppb)")
    plt.savefig('plots/cv-gc_nitrate.png')
    plt.close()

    f,ax=plt.subplots(1,1, figsize=(10,4))
    ax.plot( ds0['2010'].index, ds0['2010'], c='g', label='GEOS-Chem v13.1.2')
    ax.plot( df['2010'].index,  df['2010'], c='k', label='Obs')

    plt.legend()
    plt.ylabel("NO3 (ug m-3)")
    plt.savefig('plots/cv-gc_nitrate_2010.png')
    plt.close()

    ds0=ds0.resample('M').mean()
    df=df.resample('M').mean()
    f,ax=plt.subplots(1,1, figsize=(10,4))
    ax.plot( ds0.index, ds0, c='g', label='GEOS-Chem v13.1.2')
    ax.plot( df.index,  df, c='k', label='Obs')

    plt.legend()
    plt.ylabel("NO3 (ppb)")
    plt.savefig('plots/cv-gc_nitrate_monthly.png')
    plt.close()

if __name__=="__main__":
    main()

