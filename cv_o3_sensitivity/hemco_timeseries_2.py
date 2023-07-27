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

days=np.array([31,28,31,30,31,30,31,31,30,31,30,31])
regions={
    False : False,
    'Asia' : { 'minlon' : 70, 'maxlon' : 165, 'minlat' : 0., 'maxlat' : 50. },
    'US'   : { 'minlon' : -160, 'maxlon' : -60, 'minlat' : 15., 'maxlat' : 72. },
    'EU'   : { 'minlon' : -20, 'maxlon' : 40, 'minlat' : 36., 'maxlat' : 72. },
    'NH'   : { 'minlon' : -180, 'maxlon' : 180, 'minlat' : 0., 'maxlat' : 90. },
    'SH'   : { 'minlon' : -180, 'maxlon' : 180, 'minlat' : -90., 'maxlat' : 0. },
        }
rundirs = ['fullchem_hourly_default','fullchem_onlyCEDS','geo_ethane', 'geo_ethane_NEI']#, 
labels = ['hemco_default (Tzompa)', 'Only CEDS','CEDS + Geological','CEDS + Geo + NEI']
#            'fullchem_onlyNEI','fullchem_onlyDICE','fullchem_onlyMIX' ]
#rundirs=['fullchem_4x5_LVOCfalse']
region = False
#for rundir in rundirs:
#    print(rundir)
for region in regions:
    print(region)
    for r, rundir in enumerate(rundirs):
        path='/users/mjr583/scratch/GC/GEOS-Chem/rundirs/%s/OutputDir/' %rundir
        flist=find_file_list(path, ['HEMCO_diagnostics.201'])
        a = xr.open_mfdataset(flist)

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
            out='global'
        
        AREA=a.AREA
        #AREA = rp.surface_area_earth(46,72)
        #print(AREA.values.sum())
        em=a.sum(dim='lev', keep_attrs=True)
        em = em['EmisC2H6_Total']
        em=np.array(em)
        print(em.shape)
        days=np.array([31,28,31,30,31,30,31,31,30,31,30,31])
        new_em=[]
        for n in range(len(em)):
            x = em[n] * AREA[n] * ( 3600 * 24 * days.mean()) * 1e-9
            new_em.append(x)
        em = np.array(new_em)
        print(rundir, out+': ', np.round(em.sum(),2), 'Tg ')

        a=a.sum(dim='lat', keep_attrs=True)
        a=a.sum(dim='lon', keep_attrs=True)
        a=a.sum(dim='lev', keep_attrs=True)
        a1d = a['EmisC2H6_Total']
        
        if rundir=='fullchem_hourly_default':
            rundir='hemco_default'
        #else:
        #    label='only_ceds'
        
        a1d.plot.step("-", label=labels[r])
    plt.legend()
    plt.savefig('plots/allGEO_EMS_%s_anthro_ethane_emissions.png' %out)
    plt.close()
    #sys.exit()
