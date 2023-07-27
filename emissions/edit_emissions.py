#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#SBATCH --job-name=timeseries
#SBATCH --ntasks=1
#SBATCH --mem=10gb
#SBATCH --partition=interactive
#SBATCH --time=00:10:00
#SBATCH --output=LOGS/timeseries.log
import sys
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib
import glob
import numpy as np
from netCDF4 import Dataset
matplotlib.use('agg')
sys.path.append('/users/mjr583/python_lib')
import GC_tools as GC
import RowPy as rp
from CVAO_dict import CVAO_dict as d
import CVAO_tools as CV
plt.style.use('seaborn-darkgrid')

inpath = '/mnt/lustre/groups/chem-acm-2018/earth0_data/GEOS/HEMCO/GFED4/v2020-02/2017/'
outpath = '/users/mjr583/scratch/GC/12.9.3/rundirs/irma_025x03125_noAfrBB/Emissions/'

for infile in sorted(glob.glob('%s/GFED4_gen*' %inpath)):
    outfile=infile[-27:]
    emissions_outfile = Dataset('%s/%s' %(outpath, outfile),'w')#,format='NETCDF4_CLASSIC')

    print(infile)
    print(outfile)

    ## Read in all variables from GFED emissions file 
    fh = Dataset(infile,'r')
    lats = fh.variables['lat']
    lons = fh.variables['lon']
    T = fh.variables['time']
    
    time = emissions_outfile.createDimension('time',len(T))
    lat  =  emissions_outfile.createDimension('lat',len(lats))
    lon  = emissions_outfile.createDimension('lon',len(lons))

    time = emissions_outfile.createVariable('time',np.int64, ('time'))
    lat  =  emissions_outfile.createVariable('lat', float, ('lat'), fill_value=np.nan)
    lon  = emissions_outfile.createVariable('lon',float, ('lon'), fill_value=np.nan)

    time.standard_name=T.standard_name
    time.axis=T.axis
    time.units=T.units
    time.calendar=T.calendar
    time[:]=T[:]
    
    lat.standard_name=lats.standard_name
    lat.long_name=lats.long_name
    lat.units=lats.units
    lat.axis=lats.axis
    lat[:]=lats[:]

    lon.standard_name=lons.standard_name
    lon.long_name=lons.long_name
    lon.units=lons.units
    lon.axis=lons.axis
    lon[:]=lons[:]

    keys=fh.variables
    total_ems=np.zeros((len(time), len(lats),len(lons)))
    new_total_ems=np.zeros((len(time), len(lats),len(lons)))
    for key in keys:
        if 'DM' in key:
            em = fh.variables[key] 
            zero_ems = em[:]
            zero_ems[:,200:480,650:920] = np.nan

            new_em  = emissions_outfile.createVariable(key,float, ('time','lat','lon'), fill_value=np.nan)
            new_em.units=em.units
            new_em[:]=zero_ems[:]

            total_ems+=em[:]
            new_total_ems+=new_em[:]

    fltr=np.where(total_ems==0.)
    total_ems[fltr]=np.nan
    fltr=np.where(new_total_ems==0.)
    new_total_ems[fltr]=np.nan
    
    f,(ax1,ax2) = plt.subplots(1,2, figsize=(10,4))
    X,Y=np.meshgrid(lons, lats)
    m1=rp.get_basemap(ax=ax1)
    m2=rp.get_basemap(ax=ax2)

    im = m1.pcolormesh(X,Y, total_ems[0,:,:], zorder=100, norm=matplotlib.colors.LogNorm(), cmap='hot_r')
    i2 = m2.pcolormesh(X,Y, new_total_ems[0,:,:], zorder=100, norm=matplotlib.colors.LogNorm(), cmap='hot_r')

    f.subplots_adjust(bottom=0.1)
    cbar_ax = f.add_axes([0.13, 0.1, 0.75, 0.07])
    ax1.set_title('GFED4')
    ax2.set_title('GFED4 No African fire emissions')
    f.colorbar(im, cax=cbar_ax, orientation='horizontal')
    plt.savefig('./%s.png' %outfile)
    emissions_outfile.close()
