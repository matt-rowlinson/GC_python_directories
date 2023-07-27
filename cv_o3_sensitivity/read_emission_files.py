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
sys.path.append('/users/mjr583/python_lib/')
import RowPy as rp
from netCDF4 import Dataset
import numpy as np

def find_file_list(path, substrs):
    file_list = []
    for root, directory, files in os.walk(path):
        for f in files:
            for s in substrs:
                if s in f:
                    file_list.append(os.path.join(root, f))
    file_list.sort()
    return file_list

def tzompa():
    path='/mnt/lustre/groups/chem-acm-2018/earth0_data/GEOS/ExtData/HEMCO/C2H6_2010/v2019-06/'
    flist=find_file_list(path, ['.2x25'])

    hold=np.zeros((12,91,144))
    for n, infile in enumerate(flist):
        fh=Dataset(infile)
        anth = fh.variables['ANTHR_C2H6'][:]
        biof = fh.variables['BIOFUEL_C2H6'][:]
        
        total = anth + biof
        hold[n] = total

        alat=fh.variables['lat'][:]
        alon=fh.variables['lon'][:]

    tzompa=hold
    factor = 30.07 / (2*12.011)
    tzompa = tzompa * factor

    return tzompa, alon, alat

def ceds():
    path='/mnt/lustre/groups/chem-acm-2018/earth0_data/GEOS/ExtData/HEMCO/CEDS/v2020-08/2017/'
    flist=find_file_list(path, ['C2H6-em-total-anthro'])
    for infile in flist:
        fh=Dataset(infile)
        keys = fh.variables.keys()
        for key in keys:
            if key in ['lat','lon','time']:
                continue
            em = fh.variables[key][:]
            try:
                emissions+=em
            except:
                emissions=em
        blat=fh.variables['lat'][:]
        blon=fh.variables['lon'][:]
    ceds= emissions
    return ceds, blon, blat


def plot_emissions(em, lon, lat):
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    levels=np.arange(0,3.5, 0.2)

    f, (ax1 ax2) = plt.subplots(1,2, figsize=(14,7))
    X,Y = np.meshgrid(alon, alat)
    m = rp.get_basemap(ax=ax1)

    #im1=m.pcolormesh(X,Y,np.sum(tzompa,0), vmax=5.)
    im1=m.contourf(X,Y,np.sum(tzompa,0), levels=levels)

    divider = make_axes_locatable(ax1)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    f.colorbar(im1, cax=cax, orientation='vertical')
    ax1.set_title('Tzompa')

    X,Y = np.meshgrid(blon, blat)
    m = rp.get_basemap(ax=ax2)

    #im2=m.pcolormesh(X,Y,np.sum(ceds,0), vmax=5.)
    im2=m.contourf(X,Y,np.sum(ceds,0), levels=levels)

    divider = make_axes_locatable(ax2)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    f.colorbar(im2, cax=cax, orientation='vertical')
    ax2.set_title('CEDS')

    plt.savefig('plots/c2h6_emissions.png')
    plt.close()

# get annual totals
sae = rp.surface_area_earth(720,360)
sae=np.swapaxes(sae,0,1)
days=np.array([31,28,31,30,31,30,31,31,30,31,30,31])

new=[]
for n in range(12):
    x = ceds[n] * sae * ( 3600 * 24 * days[n]) * 1e-9
    new.append(x)
ceds = np.array(new)
print('CEDS:', np.round(np.nansum(ceds),2), 'Tg yr-1')


sae = rp.surface_area_earth(144,91)
sae=np.swapaxes(sae,0,1)
days=np.array([31,28,31,30,31,30,31,31,30,31,30,31])

new=[]
for n in range(12):
    x = tzompa[n] * sae * ( 3600 * 24 * days[n]) * 1e-9
    new.append(x)
tzompa = np.array(new)
print('Tzompa:', np.round(np.nansum(tzompa),2), 'Tg yr-1')


sys.exit()
