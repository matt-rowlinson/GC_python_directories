#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#SBATCH --job-name=gcpy_diffmap
#SBATCH --ntasks=1
#SBATCH --mem=8gb
#SBATCH --partition=nodes
#SBATCH --time=00:25:00
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

path='/users/mjr583/scratch/GC/GEOS-Chem/rundirs/fullchem_4x5_LVOCfalse/'
flist=find_file_list(path, ['SpeciesConc.201907'])
a = xr.open_mfdataset(flist)

path='/users/mjr583/scratch/GC/GCHP/rundirs/benny/GCHP_benchmark/' 
flist=find_file_list(path, ['SpeciesConc.201907'])
b = xr.open_mfdataset(flist)

a=a.mean(dim='time', keep_attrs=True)
b=b.mean(dim='time', keep_attrs=True)

species=['O3','CO', 'NO','NO2','TOLU','ISOP','HNO2','NO3','ACET','SO2', 'C2H6','CH4','C3H8','BENZ','ALD2']
titles=['O3','CO', 'NO','NO2','TOLU','ISOP','HNO2','NO3','ACET','SO2', 'C2H6','CH4','C3H8','BENZ','ALD2']

for n,s in enumerate(species):
    print(s)
    gcplot.compare_single_level(a, 'GCC_13.0.1', b, 'GCHP_13.0.1', varlist=['SpeciesConc_%s' %s],extra_title_txt=titles[n])
    plt.savefig('plots/%s_diffmap_%s.png' %('benchmark', s))

