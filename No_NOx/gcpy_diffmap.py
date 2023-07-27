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

path='/users/mjr583/scratch/GC/GEOS-Chem/rundirs/fullchem_hourly_default/'
flist=find_file_list(path, ['SpeciesConc.2017'])
a = xr.open_mfdataset(flist)

name='no_NOx'
path='/users/mjr583/scratch/GC/GEOS-Chem/rundirs/%s/' %name
flist=find_file_list(path, ['SpeciesConc.2017'])
b = xr.open_mfdataset(flist)

a=a.mean(dim='time', keep_attrs=True)
b=b.mean(dim='time', keep_attrs=True)
keys=list(a.keys())

for n,s in enumerate(keys):
    if 'N' in s:
        print(s)
        gcplot.compare_single_level(a, 'Ref', b, name, varlist=[s],extra_title_txt='N Species (%s)' %s)
        plt.savefig('Nspecies_plots/%s_diffmap_%s.png' %(name, s))

