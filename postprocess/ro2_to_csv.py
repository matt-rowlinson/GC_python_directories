#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#SBATCH --job-name=process_to_csv
#SBATCH --ntasks=1
#SBATCH --mem=6gb
#SBATCH --partition=nodes
#SBATCH --time=08:30:00
#SBATCH --output=Logs/RO_to_csv_%a.log
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

peroxy_list1=['ATO2', 'A3O2','BRO2','C4HVP1','C3HVP2','DHPCARP','DIBOO','HC5OO','ICHOO','ICNOO','IDHNBOO','IDHNDOO1','IDHNDOO2']
peroxy_list2=['IEPOXAOO','IEPOXBOO','IHOO1','IHOO4','IHPNBOO','IHPNDOO','IHPOO1','IHPOO2','IHPOO3','INO2B','INO2D','ISNOOA','ISNOOA']
peroxy_list3=['ISNOHOO','ISOPNOO1','ISOPNOO2','KO2','LIMO2','MACR1OO','MAO3','MAOPO2','MAP','MCO3','MCRHP','MCROHOO','MO2','MVKOHOO']
peroxy_list4=['NPRNO3','NRO2','OTHRO2','PIO2','PO2','PRN1','R4N1','R4O2','RCO3', 'TRO2','XRO2']
var=[] ; times=[]
counter=0

infile='/users/mjr583/scratch/GC/GEOS-Chem/rundirs/fullchem_hourly_default/OutputDir/GEOSChem.SpeciesConc.20160101_0000z.nc4'
fh=Dataset(infile)
keys=list(fh.variables.keys())
for k in keys[:]:
    var=[]
    varname = k.replace('SpeciesConc_','')
    if varname in peroxy_list4:
        print(varname)
    else:
        continue
    for n,infile in enumerate(sorted(glob.glob('/users/mjr583/scratch/GC/GEOS-Chem/rundirs/fullchem_hourly_default/OutputDir/GEOSChem.SpeciesConc*.nc4'))):
        print(infile) 
        
        fh=Dataset(infile)
        keys=list(fh.variables.keys())
        
        var.append( fh.variables[k][:,0,:,:] )
        unit=fh.variables['%s' %(k)].units
        longname=fh.variables['%s' %(k)].long_name

        if counter==0:
            time=fh.variables['time'][:]
            t0=fh.variables['time'].units
            t0=(int, re.findall(r'\d+', t0))[1]
            t0=datetime.datetime(int(t0[0]), int(t0[1]), int(t0[2]), int(t0[3]), int(t0[4]), int(t0[5]) )
            for dt in time:
                times.append( t0 + datetime.timedelta(minutes=dt) )
            time=np.array(times)
    lat= fh.variables['lat'][:]
    lon= fh.variables['lon'][:]
    variable=np.concatenate(var)

    lat_idx = rp.find_nearest(lat, 16.9)
    lon_idx = rp.find_nearest(lon, -24.9)

    variable=variable[:,lat_idx, lon_idx]

    if counter == 0:
        df=pd.DataFrame( {'%s (%s)' %(longname, unit): variable }, index=time )
        counter+=1
    else:
        df['%s (%s)' %(longname, unit)] = variable
print(counter)
df.index.name='datetime'
print(df)
df.to_csv('/users/mjr583/GC/postprocess/csv_files/GC_peroxy_radicals_4_hourly_2016-2018.csv')
