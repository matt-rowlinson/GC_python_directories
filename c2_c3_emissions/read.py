#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#SBATCH --job-name=read.py
#SBATCH --ntasks=1
#SBATCH --mem=10gb
#SBATCH --partition=interactive
#SBATCH --time=00:15:00
#SBATCH --output=Logs/read.py.log
import sys
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib
import glob
import numpy as np
from netCDF4 import Dataset
matplotlib.use('agg')
sys.path.append('/users/mjr583/python_lib')
import RowPy as rp
#plt.style.use('seaborn-darkgrid')


def read_em(source, variable, year='2010', month=False, verbose=False):
    HEMCO_path = '/mnt/lustre/groups/chem-acm-2018/earth0_data/GEOS/HEMCO/'

    if source=='CEDS':
        inpath = HEMCO_path + source + '/v2020-08/' + year
        for infile in sorted(glob.glob('%s/%s-*total-anthro*' %(inpath, variable) )):
            if verbose:
                print(infile)
            fh = Dataset(infile)
            lat=fh.variables['lat'][:]
            lon=fh.variables['lon'][:]
            keys=list(fh.variables.keys())
            em=0
            for key in keys[3:]:
                em += fh.variables[key][:]
        em = np.sum(em, 0)
        sae = np.swapaxes(rp.surface_area_earth(720, 360, verbose=False),0,1)
        em = em * sae * (3600*24)*365/12 * 1e-9
        if verbose:
            print(np.nansum(em), 'Tg yr-1')

    elif source=='Tzompa' or source.title()=='Tzompa' or source=='TZOMPA':
        inpath = HEMCO_path + '/C2H6_2010/v2019-06/'
        hold=np.zeros((12,91,144))
        for n, infile in enumerate(sorted(glob.glob('%s/C2H6*anth*2010*2x25*' %inpath))):
            if verbose:
                print(infile)
            fh = Dataset(infile)
            lat=fh.variables['lat'][:]
            lon=fh.variables['lon'][:]
            keys=list(fh.variables.keys())
            anthro = fh.variables['ANTHR_C2H6'][:]
            biof = fh.variables['BIOFUEL_C2H6'][:]
            total=anthro+biof
            hold[n]=total
        factor = 30.07 / (2*12.011)
        em = hold * factor
        sae = np.swapaxes(rp.surface_area_earth(len(lon), len(lat), verbose=False),0,1)
        new=[]
        for n in range(12):
            x = em[n] * sae * (3600 * 24 * 30.1) * 1e-9
            new.append(x)
        em = np.sum(np.array(new),0)
        if verbose:
            print(em.sum(), 'Tg yr-1')

    elif source=='NEI':
        em=0
        inpath = HEMCO_path + '/NEI2011/v2017-02-MM/'
        hold=np.zeros((12,400,900))
        for n, infile in enumerate(sorted(glob.glob('%s/08/NEI11_0.1x0.1_2011*_monmean_egupk.nc' %inpath))):
            
            print(infile)
            fh = Dataset(infile)
            lat=fh.variables['lat'][:]
            lon=fh.variables['lon'][:]
            keys=list(fh.variables.keys())
            
            if variable=='all':
                VOC_list = ['HONO', 'FORM', 'PAR', 'IOLE', 'OLE', 'ETH', 'ETHA', 'ETOH', 'MEOH', 'BENZENE', 'TOL',
                            'XYL', 'PEC', 'POC', 'PSO4', 'ACROLEIN', 'ALD2', 'ALDX', 'ISOP']
                for key in keys:
                    if key in VOC_list:
                        em += np.squeeze( fh.variables[key][:] )
                        print( key , np.sum( fh.variables[key][:] ) ) 
                total = np.sum(em, 0)
            else:
                total = np.squeeze( fh.variables[variable][:] )
            total = np.nansum(total)
            hold[n]=total
            print( np.round( total.sum() * 2. , 2) , "Tg C" )
            #sys.exit()
        em=hold
        sae = np.swapaxes(rp.surface_area_earth(len(lon), len(lat), verbose=False),0,1)
        new=[]
        for n in range(12):
            x = em[n] * sae * (3600 * 24 ) * 1e-9
            new.append(x)
        em = np.sum(np.array(new),0)
        #print(em.sum(), 'Tg yr-1')
        print( np.round( em.sum() * 2. , 2) , "Tg C" )


    elif source=='EDGAR' or source.title()=='edgar' or source=='Edgar':
        inpath = '/mnt/lustre/users/mjr583/GC/emissions/'
        for n, infile in enumerate(sorted(glob.glob('%s/v432*.nc' %inpath))):
            if verbose:
                print(infile)
            fh = Dataset(infile)
            lat=fh.variables['lat'][:]
            lon=fh.variables['lon'][:]
            
            edgar = fh.variables['emi_voc2'][:]
        sae = np.swapaxes(rp.surface_area_earth(len(lon), len(lat), verbose=False),0,1)
        
        em = edgar * sae * (3600 * 24 * 365) * 1e-9
        if verbose:
            print(em.sum(), 'Tg yr-1')

    elif source=='GFED':
        inpath = '/mnt/lustre/users/mjr583/GC/GEOS-Chem/rundirs/fullchem_hourly_default/OutputDir'
        hold=np.zeros((12,46,72))
        for n, infile in enumerate(sorted(glob.glob('%s/HEMCO_diagnostic*2018*' %inpath))):
            if verbose:
                print(infile)
            fh = Dataset(infile)
            lat=fh.variables['lat'][:]
            lon=fh.variables['lon'][:]
            keys=list(fh.variables.keys())
            total = fh.variables['EmisC2H6_BioBurn'][:]
            hold[n]=total
        em = hold
        sae = np.swapaxes(rp.surface_area_earth(len(lon), len(lat), verbose=False),0,1)
        new=[]
        for n in range(12):
            x = em[n] * sae * (3600 * 24 * 30.1) * 1e-9
            new.append(x)
        em = np.sum(np.array(new),0)
        if verbose:
            print(em.sum(), 'Tg yr-1')

    else:
        print('ERROR: Enter valid source for emission read')
        return -1
    return em, lon, lat


region_dict = {
        'Global' : {
            'lllat' : -90.,
            'urlat' : 90., 
            'lllon' : -180., 
            'urlon' : 180.
                },
        'Europe' : {
            'lllat' : 30.05,
            'urlat' : 89.95, 
            'lllon' : -29.95, 
            'urlon' : 89.95 
                },
        'North America' : {
            'lllat' : 12.7,
            'urlat' : 71.5, 
            'lllon' : -169., 
            'urlon' : -52.2
                },
        'South America' : {
            'lllat' : -56.3 ,
            'urlat' : 12.7 , 
            'lllon' : -82.9, 
            'urlon' : -33.01
                },
        'Asia' : {
            'lllat' : 0. ,
            'urlat' : 55.1 , 
            'lllon' : 50.8 , 
            'urlon' : 150.9 
                },
        'Africa' : {
            'lllat' : -35.2 ,
            'urlat' : 35.5 , 
            'lllon' : -17.3 , 
            'urlon' : 50.8 
                },
        'Oceania' : {
            'lllat' : -47.1 ,
            'urlat' : 0. , 
            'lllon' : 99.1 , 
            'urlon' : 179. 
                },
        'UK' : {
            'lllat' : 49.8 ,
            'urlat' : 59.3 , 
            'lllon' : -10.56 , 
            'urlon' : 1.9 
                },

        }
import copy
def regional_totals(em, lat ,lon, region):
    reg_em = copy.deepcopy(em)
    d = region_dict
    for ilon in range(len(lon)):
        for ilat in range(len(lat)):
            if d[region]['lllon'] < lon[ilon] <  d[region]['urlon']:
                if  d[region]['lllat'] < lat[ilat] <  d[region]['urlat']:
                    pass
                else:
                    reg_em[ilat, ilon] = np.nan
            else:
                reg_em[ilat, ilon] = np.nan
    total=np.nansum(reg_em)

    return total

def regional(em, lat ,lon, region):
    reg_em = copy.deepcopy(em)
    d = region_dict
    for ilon in range(len(lon)):
        for ilat in range(len(lat)):
            if d[region]['lllon'] < lon[ilon] <  d[region]['urlon']:
                if  d[region]['lllat'] < lat[ilat] <  d[region]['urlat']:
                    pass
                else:
                    reg_em[ilat, ilon] = np.nan
            else:
                reg_em[ilat, ilon] = np.nan
    return reg_em


