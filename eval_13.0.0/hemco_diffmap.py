#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#SBATCH --job-name=gcpy_diffmap
#SBATCH --ntasks=1
#SBATCH --mem=1gb
#SBATCH --partition=nodes
#SBATCH --time=00:00:10
#SBATCH --output=Logs/gcpy_diffmap_%A.log
import xarray as xr
import gcpy.plot as gcplot
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('agg')
import os
import sys

def find_file_list(path, substrs):
    file_list = []
    for root, directory, files in os.walk(path):
        for f in files:
            for s in substrs:
                if s in f:
                    file_list.append(os.path.join(root, f))
    file_list.sort()
    return file_list

path='/users/mjr583/scratch/GC/GEOS-Chem/rundirs/fullchem_hourly_default/OutputDir/'
flist=find_file_list(path, ['HEMCO_diagnostics.20160301'])
a = xr.open_mfdataset(flist)

path='/users/mjr583/scratch/GC/GEOS-Chem/rundirs/fullchem_hourly_default/OutputDir/'
flist=find_file_list(path, ['HEMCO_diagnostics.20160301'])
b = xr.open_mfdataset(flist)

a=a.sum(dim='time', keep_attrs=True)
b=b.sum(dim='time', keep_attrs=True)

akeys=a.keys()
bkeys=b.keys()

keys=list(set(akeys).intersection(bkeys))
for key in keys:
    if 'EmisC2H6' in  key:
        print(key)
        gcplot.compare_single_level(a, 'v12.9.3', b, 'v13.0.0', varlist=[key])
        plt.savefig('plots/HEMCO_%s.png' %key)
        plt.close()
