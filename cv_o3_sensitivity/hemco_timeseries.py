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
import gcpy.grid as grid
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
year='2016'
path='/users/mjr583/scratch/GC/GEOS-Chem/rundirs/fullchem_hourly_default/OutputDir/'
flist=find_file_list(path, ['HEMCO_diagnostics.%s' %year])
a = xr.open_mfdataset(flist)

bname='2xethane_EU'
path='/users/mjr583/scratch/GC/GEOS-Chem/rundirs/%s/OutputDir/' %bname
flist=find_file_list(path, ['HEMCO_diagnostics.%s' %year])
b = xr.open_mfdataset(flist)

cname='2xethane_US'
path='/users/mjr583/scratch/GC/GEOS-Chem/rundirs/%s/OutputDir/' %cname
flist=find_file_list(path, ['HEMCO_diagnostics.%s' %year])
c = xr.open_mfdataset(flist)

dname='2xethane_Asia'
path='/users/mjr583/scratch/GC/GEOS-Chem/rundirs/%s/OutputDir/' %dname
flist=find_file_list(path, ['HEMCO_diagnostics.%s' %year])
d = xr.open_mfdataset(flist)

region=True
if region:
    
    #minlon,maxlon,minlat,maxlat=(-20.,40.,36.,72.)  # EU
    #minlon,maxlon,minlat,maxlat=(-160.,-60.,15.,72.)  # US
    minlon,maxlon,minlat,maxlat=(70.,165.,0.,50.)  # Asia

    a = a.\
            where(a.lon>= minlon, drop=True).\
            where(a.lon<= maxlon, drop=True).\
            where(a.lat>= minlat, drop=True).\
            where(a.lat<= maxlat, drop=True)

    b = b.\
            where(b.lon>= minlon, drop=True).\
            where(b.lon<= maxlon, drop=True).\
            where(b.lat>= minlat, drop=True).\
            where(b.lat<= maxlat, drop=True)
    c = c.\
            where(c.lon>= minlon, drop=True).\
            where(c.lon<= maxlon, drop=True).\
            where(c.lat>= minlat, drop=True).\
            where(c.lat<= maxlat, drop=True)

    d = d.\
            where(d.lon>= minlon, drop=True).\
            where(d.lon<= maxlon, drop=True).\
            where(d.lat>= minlat, drop=True).\
            where(d.lat<= maxlat, drop=True)

a=a.sum(dim='lat', keep_attrs=True)
b=b.sum(dim='lat', keep_attrs=True)
c=c.sum(dim='lat', keep_attrs=True)
d=d.sum(dim='lat', keep_attrs=True)

a=a.sum(dim='lon', keep_attrs=True)
b=b.sum(dim='lon', keep_attrs=True)
c=c.sum(dim='lon', keep_attrs=True)
d=d.sum(dim='lon', keep_attrs=True)

a=a.sum(dim='lev', keep_attrs=True)
b=b.sum(dim='lev', keep_attrs=True)
c=c.sum(dim='lev', keep_attrs=True)
d=d.sum(dim='lev', keep_attrs=True)

a1d = a.EmisC2H6_Anthro
b1d = b.EmisC2H6_Anthro
c1d = c.EmisC2H6_Anthro
d1d = d.EmisC2H6_Anthro

print(a1d.values)
print(b1d.values.mean())
print(c1d.values.mean())
print(d1d.values.mean())

a1d.plot.step("-", label='Ref')
b1d.plot.step("--", label=bname)
c1d.plot.step(":",label=cname)
d1d.plot.step("-.", label=dname)

#plt.ylim()
plt.legend()
plt.savefig('plots/Asia_anthro_ethane_emissions_%s.png' %year)
plt.close()
