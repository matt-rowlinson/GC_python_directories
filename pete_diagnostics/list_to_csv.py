#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#SBATCH --job-name=process_to_csv
#SBATCH --ntasks=1
#SBATCH --mem=13gb
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

for infile in glob.glob('/users/mjr583/scratch/GC/13.1.2_diagnostics/rundirs/test/OutputDir/GEOSChem.ConcAfterChem*20190701*.nc4' ):
        fh=Dataset(infile)
        keys=fh.variables.keys()#[11:]
        continue

variables = { 
                'CH4' : {
                        'Collection' : 'ConcAfterChem',
                        'varname' : 'concAfterChem_LVOC'},
               'CO' : {
                        'Collection' : 'ConcAfterChem',
                        'varname' : 'concAfterChem_CO'},
               'O3' : {
                        'Collection' : 'ConcAfterChem',
                        'varname' : 'concAfterChem_O3'},
               'C2H6' : {
                        'Collection' : 'ConcAfterChem',
                        'varname' : 'concAfterChem_C2H6'},
               'OH' : {
                        'Collection' : 'ConcAfterChem',
                        'varname' : 'concAfterChem_OH'},
               'TOLU' : {
                        'Collection' : 'ConcAfterChem',
                        'varname' : 'concAfterChem_TOLU'},
               'BENZ' : {
                        'Collection' : 'ConcAfterChem',
                        'varname' : 'concAfterChem_BENZ'},
                        }

def basic_map_plot(var, lon, lat, varname):
    f,ax=plt.subplots()
    m=rp.get_basemap(ax=ax)
    X,Y=np.meshgrid(lon,lat)
    from matplotlib.colors import LogNorm 
    ax.pcolormesh(X,Y,var)#, norm=LogNorm())
    plt.savefig('plots/%s.png' %varname)

#for n,v in enumerate(variables):

for n, key in enumerate(keys):
    if 'concAfterChem' not in key:
        continue
    #collection=variables[v]['Collection']
    #varname=variables[v]['varname']
    varname=key
    print(varname)

    var=[] ; times=[] ; air=[] 
    for i,infile in enumerate(sorted(glob.glob('/users/mjr583/scratch/GC/13.1.2_diagnostics/rundirs/test/OutputDir/GEOSChem.ConcAfterChem*20190701*.nc4'))):
        metfile = infile.replace( 'ConcAfterChem', 'StateMet')
        air.append( Dataset(metfile).variables['Met_AIRVOL'][:] )
        
        fh=Dataset(infile)
        var.append( fh.variables['%s' %(varname) ][:,:,:,:] )
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
    airvol=np.concatenate(air)

    var = var * airvol[:]

    var = np.mean(var,0)[0]

    f,ax=plt.subplots()
    m=rp.get_basemap(ax=ax)
    X,Y=np.meshgrid(lon,lat)
    im = ax.pcolormesh(X,Y,var)

    cax = f.add_axes([0.27, 0.1, 0.5, 0.05])
    cbar = f.colorbar( im, cax=cax, orientation='horizontal')
    cbar.ax.set_xlabel('molec')
    ax.set_title(longname)
    if var.max() == 0.:
        varname += '_ZERO'
    plt.savefig('plots/%s.png' %varname)
    plt.close()
    

