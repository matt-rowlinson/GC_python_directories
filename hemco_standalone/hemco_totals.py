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
import pandas as pd
import matplotlib
matplotlib.use('agg')
import os
import numpy as np
import sys
sys.path.append('/users/mjr583/python_lib/')
import RowPy as rp

def find_file_list(path, substrs):
    file_list =[]
    for root, directory, files in os.walk(path):
        for f in files:
            for s in substrs:
                if s in f:
                    file_list.append(os.path.join(root, f))
    file_list.sort()
    return file_list

def annual_means( em , times=False):
    try:
        em = np.sum(em.reshape(-1, 12), axis=1)
    except:
        for i in range(1,12):
            try:
                em = np.sum(em[:-i].reshape(-1, 12), axis=1)
                break
            except:
                continue
    if type(times) != bool:
        times = pd.to_datetime(times)
        times = times.year.drop_duplicates().tolist()
        return em, times
    return em


regions={
    False : False,
    #'Asia' : { 'minlon' : 70, 'maxlon' : 165, 'minlat' : 0., 'maxlat' : 50. },
    #'US'   : { 'minlon' : -160, 'maxlon' : -60, 'minlat' : 15., 'maxlat' : 72. },
    #'EU'   : { 'minlon' : -30, 'maxlon' : 90, 'minlat' : 30., 'maxlat' : 90. },
    #'NH'   : { 'minlon' : -180, 'maxlon' : 180, 'minlat' : 0., 'maxlat' : 90. },
    #'SH'   : { 'minlon' : -180, 'maxlon' : 180, 'minlat' : -90., 'maxlat' : 0. },
        }


if __name__=="__main__":
    rundirs = ['hemco_standalone']
    labels = ['Check']
    key='NO'
    region = False
    days=np.array([31,28,31,30,31,30,31,31,30,31,30,31])
    
    light_nox=[] ; times=[]
    for region in regions:
        for r, rundir in enumerate(rundirs):
            path=f'/users/mjr583/scratch/GC/13.1.2/rundirs/{rundir}/OutputDir/'
            flist=find_file_list(path, ['HEMCO_sa_diagnostics.20'])
            #print(flist)

            for f in flist:
                print( f )
                fh = xr.open_dataset(f)
                AREA=fh.AREA
                time=fh.coords['time'].values

                em = np.squeeze( fh[f'Emis{key}_Lightning'] )
                em=em.sum(dim='lev', keep_attrs=True)
                 
                em= np.array(em)

                em = em * AREA * ( 3600 * 24 * 365/12 ) * 1e-9
                total = np.nansum( em ) * 0.46680 
                light_nox.append( total)
                times.append( time )

            
    light_nox = np.array(light_nox)
    times = np.array( times )

    f, ax = plt.subplots(1,1,figsize=(15,5))
    plt.plot( times, light_nox )
    
    plt.ylabel( "Lightning $NO_x$ (Tg N)")
    plt.savefig('plots/mon_mean_lightnox.png')
    plt.close()


    annual_light_nox = annual_means( light_nox )
    times = np.squeeze( times[0::12] )

    f, ax = plt.subplots(1,1,figsize=(15,5))
    plt.plot( times, annual_light_nox )
    
    plt.ylabel( "Lightning $NO_x$ (Tg N)")
    plt.savefig('plots/ann_mean_lightnox.png')
    plt.close()

            
            
            
'''
            a = xr.open_mfdataset(flist, combine='by_coords')
            print("Success")

            
            AREA=a.AREA
            em=a.sum(dim='lev', keep_attrs=True)
            em = np.array( em[f'Emis{key}_Lightning'] )
            days=np.array([31,28,31,30,31,30,31,31,30,31,30,31])
            new_em=[]
            print(em.shape)
            for n in range(len(em)):
                x = em[n] * AREA[n] * ( 3600 * 24 * days[n]) * 1e-9
                new_em.append(x)
            em = np.array(new_em)
            print(rundir, out+': ', np.round(em.sum(),2) * 0.46680, 'Tg ')

            a=a.sum(dim='lat', keep_attrs=True)
            a=a.sum(dim='lon', keep_attrs=True)
            a=a.sum(dim='lev', keep_attrs=True)
            a1d = a['Emis%s_Lightning' %key]
            
            
            a1d.plot.step("-", label=labels[r])
        plt.legend()
        plt.savefig(f'plots/{rundir}_{key}__emissions.png')
        plt.close()

'''
