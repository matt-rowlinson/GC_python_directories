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
import matplotlib
matplotlib.use('agg')
sys.path.append('/users/mjr583/python_lib')
import RowPy as rp
from CVAO_dict import CVAO_dict as d
import CVAO_tools as CV
plt.style.use('seaborn-darkgrid')
import numpy as np 
import xarray as xr
import os

def find_file_list(path, substrs):
    file_list =[]
    for root, directory, files in os.walk(path):
        for f in files:
            for s in substrs:
                if s in f:
                    file_list.append(os.path.join(root, f))
    file_list.sort()
    return file_list

def get_data_as_xr(rundir, year='', lev=0, version='14.0.1'):
    path=f'/mnt/lustre/users/mjr583/GC/{version}/rundirs/{rundir}'
    ds = xr.open_mfdataset( find_file_list(path, [f'SpeciesConc.{year}']), combine='by_coords' )
    return ds

def output_stats(var, obs, model, rundirs):
    first=True
    for rundir, m in zip(rundirs,model):
        m = np.squeeze(m).values
        o = np.squeeze(obs).values

        absError= m - o
        SE = np.square(absError)
        MSE= np.mean(SE)
        RMSE=np.round( np.sqrt(MSE/len(o)), 3)
        R2 = np.round( 1. - (np.var(absError) / np.var(o)), 3)
        nmb = 100 * (np.sum( absError - o ) / np.sum( o ))
        bias = ((m-o) / o).mean() * 100
        if first:
            df = pd.DataFrame( {'Mean AbsError': absError.mean(), 'bias %': bias, 'NMB':nmb,'r$^2$':R2,'RMSE':RMSE},
                    index=[rundir]    )
            first=False
        else:
            d  = pd.DataFrame( {'Mean AbsError': absError.mean(), 'bias %': bias, 'NMB':nmb,'r$^2$':R2,'RMSE':RMSE},
                    index=[rundir]    )
            df = pd.concat([df,d], axis=0, ignore_index=True)
    df=df.transpose()
    df.columns = rundirs
    to_table( var, df.round(2) )
    return df 

def to_table(var, df):
    fig, ax =plt.subplots(figsize=(3.5,1.5))
    ax.axis('tight')
    ax.axis('off')
    rcolors = plt.cm.BuPu(np.full(len(df.index), 0.1))
    ccolors = plt.cm.BuPu(np.full(len(df.columns), 0.1))
    the_table = ax.table(cellText=df.values,
            rowLabels=df.index,
            rowColours=rcolors,
            rowLoc='center',
            colLoc='center',
            cellLoc='center',
            colColours=ccolors,  
            colLabels=df.columns,
            loc='center')
    plt.title(var+': Model performance stats')
    plt.tight_layout
    plt.savefig(f'Output/cvao.{var}.stats.png', dpi=300, bbox_inches='tight')
    plt.close()
    return

####-------------------------------------------------------

def main():
    rundirs=['geo_2x25','all_2x25']
    variables=['C2H6','C3H8','O3','TOLU','BENZ','PRPE', 'ALK4']

    df_=[]
    for r in rundirs:
        print( r )
        ds = get_data_as_xr(r, '2017')
        y = rp.find_nearest(ds.lat, 16.9)
        x = rp.find_nearest(ds.lon, -24.9)
        ds = ds.isel(lat=y).isel(lon=x).isel(lev=0)
        dfs=[]
        for variable in variables:
            df = (ds[f'SpeciesConc_{variable}'] * 1e12).to_dataframe()
            df.index = pd.to_datetime( df.index, format='%Y-%m-%d %H:%M:%S').strftime('%Y-%m-%d')
            df.index = pd.to_datetime( df.index, format='%Y-%m-%d')
            dfs.append( df )
        
        df=pd.concat(dfs, axis=1, ignore_index=False)
        df = df.loc[:,~df.columns.duplicated()].copy()
        df.columns = df.columns.str.replace("SpeciesConc_","")
        df_.append( df )

    for n, var in enumerate(variables): 
        print( var )
        df = CV.get_from_merge(d[var])
        df = df.loc['2017'].resample( 'D' ).mean()
        df = df.dropna()
        
        v1 =  df_[0][var]
        v2 =  df_[1][var]
        v1 = v1.reindex( df.index )
        v2 = v2.reindex( df.index )

        output_stats( var, df, [v1,v2], ['Base','Scaled Emissions'] )
        
if __name__=="__main__":
    main()
