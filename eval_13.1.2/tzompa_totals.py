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
region_dict={
    False : False,
    'Asia' : { 'minlon' : 70, 'maxlon' : 165, 'minlat' : 0., 'maxlat' : 50. },
    'US'   : { 'minlon' : -160, 'maxlon' : -60, 'minlat' : 15., 'maxlat' : 72. },
    'EU'   : { 'minlon' : -30, 'maxlon' : 90, 'minlat' : 30., 'maxlat' : 90. },
    'NH'   : { 'minlon' : -180, 'maxlon' : 180, 'minlat' : 0., 'maxlat' : 90. },
    'SH'   : { 'minlon' : -180, 'maxlon' : 180, 'minlat' : -90., 'maxlat' : 0. },
        }

def get_region(a, region, regions=region_dict):
    minlon,maxlon,minlat,maxlat=(regions[region]['minlon'],regions[region]['maxlon'], 
            regions[region]['minlat'], regions[region]['maxlat'])
    a = a.\
            where(a.lon>= minlon, drop=True).\
            where(a.lon<= maxlon, drop=True).\
            where(a.lat>= minlat, drop=True).\
            where(a.lat<= maxlat, drop=True)
    return a


def convert_em_to_tg(a, key='C2H6'):
    AREA=a.AREA
    em=a.sum(dim='lev', keep_attrs=True)
    em = em['Emis%s_Anthro' %key]
    em=np.array(em)
    days=np.array([31,28,31,30,31,30,31,31,30,31,30,31])
    new_em=[]
    for n in range(len(em)):
        x = em[n] * AREA[n] * ( 3600 * 24 * days[n]) * 1e-9
        new_em.append(x)
    em = np.array(new_em)
    return em

def print_tot(em):
    print( rundir, region, ': ', np.round(em.sum(),2), 'Tg')
    return 

def open_diagnostics( path, rundir, year='2017'):
    flist = find_file_list( path, [f'HEMCO_diagnostics.{year}'])
    a = xr.open_mfdataset( flist, combine='by_coords')
    return a

if __name__=="__main__":
    ### Set these options
    rundir='control'
    key='C2H6'
    region='EU'
    path= f'/mnt/lustre/users/mjr583/GC/13.1.2/rundirs/{rundir}/OutputDir/'
    
    ### Get Tzompa files
    import read_tzompa_xiao as tz_xi
    lonlat_region = (region_dict[region]['minlon'], region_dict[region]['maxlon'],
                        region_dict[region]['minlat'], region_dict[region]['maxlat'])
    Tz, Tlon, Tlat = tz_xi.get('tzompa', lonlat=(lonlat_region), only_land=False )
    Tz = get_region( a , region=region, regions=region_dict)
    print( np.nansum( Tz * 1e-3) , 'Tg' )
    
    ### Get HEMCO_diag
    a = open_diagnostics( path, rundir, year='2017' )
    a = get_region(a, region=region, regions=region_dict)
    em = convert_em_to_tg(a, key=key)
    print( np.nansum( em ), 'Tg' )




