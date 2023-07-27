#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#SBATCH --job-name=hemco_timeseries
#SBATCH --ntasks=1
#SBATCH --mem=1gb
#SBATCH --partition=nodes
#SBATCH --time=00:00:10
#SBATCH --output=Logs/hemco_timeseries_%A.log
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('agg')
import os
import numpy as np
import sys
import pandas as pd
sys.path.append('/users/mjr583/python_lib/')
import RowPy as rp
plt.style.use('seaborn-darkgrid')

def find_file_list(path, substrs):
    file_list =[]
    for root, directory, files in os.walk(path):
        for f in files:
            for s in substrs:
                if s in f:
                    file_list.append(os.path.join(root, f))
    file_list.sort()
    return file_list

days=np.array([31,28,31,30,31,30,31,31,30,31,30,31])
regions={
    False : False,
    #'China' : { 'minlon' : 70, 'maxlon' : 165, 'minlat' : 0., 'maxlat' : 50. },
    'US'   : { 'minlon' : -139.95, 'maxlon' : -50.05, 'minlat' : 20.05, 'maxlat' : 59.95 },
    'EU'   : { 'minlon' : -30, 'maxlon' : 90, 'minlat' : 50., 'maxlat' : 90. },
    'EU2'   : { 'minlon' : -30, 'maxlon' : 73.05, 'minlat' : 29.95, 'maxlat' : 50.05 },
    #'NH'   : { 'minlon' : -180, 'maxlon' : 180, 'minlat' : 0., 'maxlat' : 90. },
    #'SH'   : { 'minlon' : -180, 'maxlon' : 180, 'minlat' : -90., 'maxlat' : 0. },
        }
rundirs = ['all_scaled']#, 'asia-meic_scale']
labels = ['?']
key='C3H8'
region = False

for region in regions:
    print(region)
    for r, rundir in enumerate(rundirs):
        path=f'/users/mjr583/scratch/GC/13.4.0/rundirs/{rundir}/OutputDir/'
        flist=find_file_list(path, ['HEMCO_diagnostics.2016'])
        a = xr.open_mfdataset(flist, combine='by_coords')

        if region:
            out=region
            minlon,maxlon,minlat,maxlat=(regions[region]['minlon'],regions[region]['maxlon'], 
                                         regions[region]['minlat'], regions[region]['maxlat'])
            a = a.\
                    where(a.lon>= minlon, drop=True).\
                    where(a.lon<= maxlon, drop=True).\
                    where(a.lat>= minlat, drop=True).\
                    where(a.lat<= maxlat, drop=True)
        else:
            out='Global'
        
        AREA=a.AREA
        HEMCO=a.sum(dim='lev', keep_attrs=True)
        

        h=[] ; tops=[]
        for k in HEMCO.variables.keys():
            if key in k:
                #print( k )
                pass
            else: 
                continue
            em = HEMCO[k]
            em=np.array(em)
            days=np.array([31,28,31,30,31,30,31,31,30,31,30,31])
            new_em=[]
            for n in range(len(em)):
                x = em[n] * AREA[n] * ( 3600 * 24 * days[n]) * 1e-9
                new_em.append(x)
            em = np.array(new_em)
            
            total = np.round(em.sum(),2)  
            h.append( total ) 
            tops.append( k )
        if r==0:
            df = pd.DataFrame( {rundir:h},index=tops) 
        else:
            df[rundir] = h
    df.columns = ['Base CEDS']#,'Scaled CEDS EOH','Scaled CEDS EOH Upper']
    df.index = df.index.str.replace('EmisEOH_','')
    print( df )
    ax = df.plot.bar(rot=0.)
    plt.ylabel( 'Tg yr$^{-1}$' )
    plt.savefig(f'plots/{out}_{key}_emissions.png')
    plt.close()
    sys.exit()
