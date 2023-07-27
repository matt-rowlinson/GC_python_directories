#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#SBATCH --job-name=strmfunction
#SBATCH --ntasks=1
#SBATCH --mem=10gb
#SBATCH --partition=interactive
#SBATCH --time=00:30:00
#SBATCH --output=LOGS/strmfunc.log
import cartopy.crs as ccrs
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from cartopy.util import add_cyclic_point
import matplotlib as mpl
import matplotlib.pyplot as plt
from netCDF4 import Dataset

mpl.rcParams['mathtext.default'] = 'regular'
import sys
sys.path.append('/users/mjr583/python_lib')
import cartopy_GC_tools as GC
import pandas as pd
import numpy as np
from matplotlib.colors import LogNorm

rundir='2018_4x5' ; version='12.9.3' ; jobid=0
month_lengths=[0,31,59,90,120,151,181,212,243,273,304,334,365]
variable='O3'
var, lat, lon, lev, time = GC.get_gc_var(rundir, variable, version)

hold=[]
for x in range(12):
    print(x)
    in_days=pd.date_range(time[month_lengths[x]], time[month_lengths[x+1]-1], freq='D')
    tp, mlat,mlon=GC.get_all_gc_input(in_days, var='PRECLSC', filetype='A1')
    tp=np.sum(tp, 0)

    #cosweights = np.cos(np.deg2rad(mlat))
    #area = np.reshape((np.tile(cosweights,1152)), (721,1152), order='f')
    
    #tp = tp * area
    #tp = tp * 86400

    mlat[0]=-90.
    mlat[-1]=90.
    fltr = mlon>=0.
    lons = np.concatenate((mlon[fltr],(mlon[~fltr] + 360.)))
    tp = np.concatenate((tp[:,fltr], tp[:,~fltr]), axis=-1)
    lats=mlat[::-1]
    tp=tp[::-1,:]

    X,Y=np.meshgrid(lons, lats)
    months=['01','02','03','04','05','06','07','08','09','10','11','12']

    # Plot total precip
    ax1 = plt.axes(projection=ccrs.PlateCarree(central_longitude=180))
    Max=np.round(np.nanmax(tp),2)
    levels=np.arange(.1,Max, .01)
    tot = ax1.contourf(X,Y, tp, cmap=plt.cm.Blues)
    #tot = ax1.contourf(X,Y, tp, cmap=plt.cm.Blues, levels=levels, norm=LogNorm(), extend='min')
    #tot = ax1.pcolormesh(X,Y,tp, cmap=plt.cm.Blues, vmin=.1)
    tot.cmap.set_under('w')
    ax1.coastlines()
    ax1.gridlines()
    ax1.set_xticks([0, 60, 120, 180, 240, 300, 359.99], crs=ccrs.PlateCarree())
    ax1.set_yticks([-90, -60, -30, 0, 30, 60, 90], crs=ccrs.PlateCarree())
    lon_formatter = LongitudeFormatter(zero_direction_label=True,
                                               number_format='.0f')
    lat_formatter = LatitudeFormatter()
    ax1.xaxis.set_major_formatter(lon_formatter)
    ax1.yaxis.set_major_formatter(lat_formatter)
    plt.colorbar(tot,orientation='horizontal')
    plt.title('Total Precipitation (kg $m^{-2}$ $s^{-1}$)', fontsize=16)
    plt.savefig('plots/gcmet_total_precipitation_%s.png' %months[x])
    plt.close()

    hold.append(tp)

tp=np.array(hold)
print(tp.shape)
outfile= Dataset('/users/mjr583/scratch/GC/12.9.3/2018_4x5/metfiles/total_precipitation.nc','w',format='NETCDF4' )

# Write dimensions (lon, lat, time)
month    = outfile.createDimension('months', len(tp))
lat     = outfile.createDimension('lat',len(lats))
lon     = outfile.createDimension('lon',len(lons))

# Write variables
mons    = outfile.createVariable('months',str, ('months'))
lat     = outfile.createVariable('lat',np.float32, ('lat'))
lon     = outfile.createVariable('lon',np.float32, ('lon'))

mons[:]     = np.array(months[:len(tp)])
lat[:]      = lats
lon[:]      = lons

precipitation = outfile.createVariable('PRECTOT',np.float,('months','lat','lon'))
precipitation.units = 'kg m-2 s-1 '
precipitation.longname = 'Total monthly precipiation'
precipitation[:] = tp

outfile.close()
