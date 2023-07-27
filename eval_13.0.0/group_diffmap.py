#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#SBATCH --job-name=gcpy_diffmap
#SBATCH --ntasks=1
#SBATCH --mem=8gb
#SBATCH --partition=nodes
#SBATCH --time=00:10:00
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

path='/users/mjr583/scratch/GC/12.9.3/rundirs/tropchem_merra_4x5/'
flist=find_file_list(path, ['SpeciesConc.2016'])
a = xr.open_mfdataset(flist)

path='/users/mjr583/scratch/GC/GEOS-Chem/rundirs/fullchem_4x5_LVOCfalse/'
flist=find_file_list(path, ['SpeciesConc.2016'])
b = xr.open_mfdataset(flist)

a=a.mean(dim='time', keep_attrs=True)
b=b.mean(dim='time', keep_attrs=True)

a['SpeciesConc_NO'] += a['SpeciesConc_NO2']
b['SpeciesConc_NO'] += b['SpeciesConc_NO2']

gcplot.compare_single_level(a, 'v12.9.3', b, 'v13.0.0', varlist=['SpeciesConc_NO'], extra_title_txt='NOx')
plt.savefig('plots/surfacemap_annmean_NOx.png' )

