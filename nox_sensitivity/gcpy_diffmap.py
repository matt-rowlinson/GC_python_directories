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
print('get control')
path='/users/mjr583/scratch/GC/GEOS-Chem/rundirs/fullchem_hourly_default/'
flist=find_file_list(path, ['SpeciesConc.2017'])
a = xr.open_mfdataset(flist)
print('get dev')
path='/users/mjr583/scratch/GC/GEOS-Chem/rundirs/fullchem_2xEU_ethane/'
flist=find_file_list(path, ['SpeciesConc.2017'])
b = xr.open_mfdataset(flist)

a=a.mean(dim='time', keep_attrs=True)
b=b.mean(dim='time', keep_attrs=True)

species=['C2H6','O3','CO']#,'CH4','NO','NO2','C2H6','C3H8', 'BENZ','TOLU','ACET','ALD2','PRPE','MOH']
titles=['C2H6','O3','CO']#,'CH4','NO','NO2','C2H6','C3H8','Benzene','Toluene','Acetone','Acetaldehyde','Propene','Methanol']
print('plot by species')
for n,s in enumerate(species):
    print(s)
    gcplot.compare_single_level(a, 'v12.9.3', b, 'v13.0.0', varlist=['SpeciesConc_%s' %s], extra_title_txt=titles[n])
    plt.savefig('plots/surfacemap_annmean_%s.png' %s)

