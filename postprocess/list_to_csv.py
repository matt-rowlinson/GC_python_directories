#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#SBATCH --job-name=process_to_csv
#SBATCH --ntasks=1
#SBATCH --mem=3gb
#SBATCH --partition=interactive
#SBATCH --time=01:30:00
#SBATCH --output=Logs/to_csv.log
import glob
from netCDF4 import Dataset
import numpy as np
import pandas as pd
import sys
import datetime
import re
import matplotlib.pyplot as plt

sys.path.append('/users/mjr583/python_lib')
import GC_tools as GC
import RowPy as rp

variables = { 'O3' : {
                        'Collection' : 'SpeciesConc',
                        'varname' : 'SpeciesConc_O3'},
                'NO' : {
                        'Collection' : 'SpeciesConc',
                        'varname' : 'SpeciesConc_NO'},
                'NO2' : {
                        'Collection' : 'SpeciesConc',
                        'varname' : 'SpeciesConc_NO2'},
                'H2O' : {
                        'Collection' : 'SpeciesConc',
                        'varname' : 'SpeciesConc_H2O'},
                'OH' : {
                        'Collection' : 'ConcAfterChem',
                        'varname' : 'OHconcAfterChem'},
                'HO2' : {
                        'Collection' : 'ConcAfterChem',
                        'varname' : 'HO2concAfterChem'},
                'JO1d' : {
                        'Collection' : 'JValues',
                        'varname' : 'JvalO3O1D'},
                'JNO2' : {
                        'Collection' : 'JValues',
                        'varname' : 'Jval_NO2'}
                        }

for n,v in enumerate(variables):
    #if n != 0:
    #    if collection ==  :
    collection=variables[v]['Collection']
    varname=variables[v]['varname']

    var=[] ; times=[] 
    for i,infile in enumerate(sorted(glob.glob('/users/mjr583/scratch/GC/GEOS-Chem/rundirs/fullchem_hourly_default/OutputDir/GEOSChem.*%s*.nc4' %collection))):
        print(infile) 
        
        fh=Dataset(infile)
        var.append( fh.variables['%s' %(varname) ][:,0,:,:] )
        unit=fh.variables['%s' %(varname)].units
        longname=fh.variables['%s' %(varname)].long_name
        
        if n==0:
            time=fh.variables['time'][:]
            t0=fh.variables['time'].units
            t0=(int, re.findall(r'\d+', t0))[1]
            t0=datetime.datetime(int(t0[0]), int(t0[1]), int(t0[2]), int(t0[3]), int(t0[4]), int(t0[5]) )
            for dt in time:
                times.append( t0 + datetime.timedelta(minutes=dt) )
            time=np.array(times)
    lat= fh.variables['lat'][:]
    lon= fh.variables['lon'][:]
    var=np.concatenate(var)

    lat_idx = rp.find_nearest(lat, 16.9)
    lon_idx = rp.find_nearest(lon, -24.9)

    var=var[:,lat_idx, lon_idx]

    if n == 0:
        df=pd.DataFrame( {'%s (%s)' %(longname, unit): var }, index=time )
    else:
        df['%s (%s)' %(longname, unit)] = var
df.index.name='datetime'
print(df)
df.to_csv('/users/mjr583/GC/postprocess/csv_files/GC_variables_hourly_2016-2018_with_H2O.csv')
