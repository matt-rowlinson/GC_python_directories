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
import matplotlib.dates as mdates
import matplotlib
matplotlib.use('agg')
sys.path.append('/users/mjr583/python_lib')
import GC_tools as GC
import RowPy as rp
from CVAO_dict import CVAO_dict as d
import CVAO_tools as CV
plt.style.use('seaborn-darkgrid')

def get_inputs():
    inputs=GC.get_arguments()
    variable=inputs.var
    if variable==None:
        variable='O3'
    return variable

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
    cv_df=pd.concat(dfs, axis=1, ignore_index=False)
    return dfs

def plot_dataframes(dfs, labels, std=False,
                         styles=['solid','dotted','dashed','dashdot','dotted','dashed','dashdot'] ):
    if 'NO' in variable:
        LOG=True
    else:
        LOG=False
    labels=rundirs
    styles=['solid','dotted','dashed','dashdot','dotted','dashed','dashdot']
    rp.timeseries_plot(dfs, variable, cvao_data=False,  string=string, labels=labels, styles=styles, logscale=LOG)
    if std:
        rp.timeseries_plot(dfs, variable, cvao_data=False, std=True, resample='D', string=string, labels=labels, styles=styles, logscale=LOG)


if __name__=="__main__":
    variable=get_inputs()
    year='2009'
    collection='SpeciesConc'
    rundirs=[f'nitrate_photol_2010_control',f'dec_new_base','dev_both']
             #f'nitrate_photol_2010_scale-50',f'nitrate_photol_2010_scale-100']
    string=f'{year}_'

    dfs = get_model_as_dataframe( rundirs, variable, collection=collection, year=year )
    plot_dataframes(dfs, labels=rundirs, std=True)
