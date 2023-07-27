#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#SBATCH --job-name=strmfunction
#SBATCH --ntasks=1
#SBATCH --mem=10gb
#SBATCH --partition=interactive
#SBATCH --time=00:30:00
#SBATCH --output=LOGS/strmfunc.log
"""
Compute streamfunction and velocity potential from the long-term-mean
flow.

This example uses the standard interface.

Additional requirements for this example:

* netCDF4 (http://unidata.github.io/netcdf4-python/)
* matplotlib (http://matplotlib.org/)
* cartopy (http://scitools.org.uk/cartopy/)

"""
import cartopy.crs as ccrs
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from cartopy.util import add_cyclic_point
import matplotlib as mpl
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from windspharm.standard import VectorWind
from windspharm.tools import prep_data, recover_data, order_latdim
from windspharm.examples import example_data_path
mpl.rcParams['mathtext.default'] = 'regular'
import sys
sys.path.append('/users/mjr583/python_lib')
import cartopy_GC_tools as GC
import pandas as pd
import numpy as np

## Read in total precipitation values
infile='/users/mjr583/scratch/GC/12.9.3/2018_4x5/metfiles/total_precipitation.nc'
fh=Dataset(infile)
tp=fh.variables['PRECTOT'][:]
lat=fh.variables['lat'][:]
lon=fh.variables['lon'][:]
fltr=np.where(tp <= 5e4)
tp[fltr]=np.nan
print(tp.shape)

rundir='2018_4x5' ; version='12.9.3' ; jobid=0
variable='O3'
var, lat, lon, lev, time = GC.get_gc_var(rundir, variable, version)
in_months = pd.date_range(time[0], time[-1], freq='MS')

U, mlat, mlon = GC.get_all_gc_input(in_months, var='U', filetype='A3dyn')
V, mlat, mlon = GC.get_all_gc_input(in_months, var='V', filetype='A3dyn')
hpa='700-300hPa'
#U=U[:,17,:,:] # 22=510hPA, 27=314hPa, 17=695hPa
#V=V[:,17,:,:]
U=np.mean(U[:,17:28,:,:],1) # 22=510hPA, 27=314hPa, 17=695hPa
V=np.mean(V[:,17:28,:,:],1)

mlat[0]=-90.
mlat[-1]=90.

fltr = mlon>=0.
lons = np.concatenate((mlon[fltr],(mlon[~fltr] + 360.)))
U = np.concatenate((U[:,:,fltr],U[:,:,~fltr]), axis=-1)
V = np.concatenate((V[:,:,fltr],V[:,:,~fltr]), axis=-1)

lats=mlat[::-1]
uwnd = U[:,::-1,:]
vwnd = V[:,::-1,:]

# The standard interface requires that latitude and longitude be the leading
# dimensions of the input wind components, and that wind components must be
# either 2D or 3D arrays. The data read in is 3D and has latitude and
# longitude as the last dimensions. The bundled tools can make the process of
# re-shaping the data a lot easier to manage.
uwnd, uwnd_info = prep_data(uwnd, 'tyx')
vwnd, vwnd_info = prep_data(vwnd, 'tyx')

# It is also required that the latitude dimension is north-to-south. Again the
# bundled tools make this easy.
lats, uwnd, vwnd = order_latdim(lats, uwnd, vwnd)

# Create a VectorWind instance to handle the computation of streamfunction and
# velocity potential.
w = VectorWind(uwnd, vwnd)

# Compute the streamfunction and velocity potential. Also use the bundled
# tools to re-shape the outputs to the 4D shape of the wind components as they
# were read off files.
sf, vp = w.sfvp()
sf = recover_data(sf, uwnd_info)
vp = recover_data(vp, uwnd_info)

## Read in total precipitation values
infile='/users/mjr583/scratch/GC/12.9.3/2018_4x5/metfiles/total_precipitation.nc'
fh=Dataset(infile)
tp=fh.variables['PRECTOT'][:]
lat=fh.variables['lat'][:]
lon=fh.variables['lon'][:]
fltr=np.where(tp <= 5e4)
tp[fltr]=np.nan
print(tp.shape)


# Cycle through months and plot ITCZ approximation for each
mons=['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
months=['01','02','03','04','05','06','07','08','09','10','11','12']
days_in_month = [31,28,21,30,31,30,31,31,30,31,30,31]
for i in range(12):
    print(months[i])
    precip=tp[i] * days_in_month[i] * 86400
    print(precip.shape)

    sf_dec, lons_c = add_cyclic_point(sf[i], lons)
    vp_dec, lons_c = add_cyclic_point(vp[i], lons)
    preci_, lons_c = add_cyclic_point(precip, lons)
    # Plot streamfunction.
    ax1 = plt.axes(projection=ccrs.PlateCarree(central_longitude=180))
    clevs = np.arange(-120,130,10)#[-80,-70, -60,-50, -40, -30, -20, -10, 0, 10, 20, 30, 40, 50, 60, 70, 80]
    sf_fill = ax1.contourf(lons_c, lats, sf_dec * 1e-06, 
                                   transform=ccrs.PlateCarree(), cmap=plt.cm.PiYG,
                                                          extend='both', zorder=2)

    levels=[1e-1, 1e1, 1e3,5e3,1e4,5e4, 1e5, 5e5]
    colors=['#ffffd9','#edf8b1','#c7e9b4','#7fcdbb','#41b6c4','#1d91c0','#225ea8','#253494','#081d58']
    PRECTOT = ax1.contourf(lons_c,lats, preci_, colors=colors,levels=levels,extend='min',zorder=4)#norm=LogNorm(), extend='min')

    #PRECTOT = ax1.contourf(lons_c, lats, preci_, levels=np.arange(0.1, preci_.max(), 0.1), cmap='Blues', zorder=3)
    itcz = ax1.contour(lons_c, lats, sf_dec * 1e-06, levels=[ 0.],  transform=ccrs.PlateCarree(), 
                                                          linestyle='--', colors='yellow', zorder=5)
    ax1.coastlines()
    ax1.gridlines()
    ax1.set_xticks([0, 60, 120, 180, 240, 300, 359.99], crs=ccrs.PlateCarree())
    ax1.set_yticks([-90, -60, -30, 0, 30, 60, 90], crs=ccrs.PlateCarree())
    lon_formatter = LongitudeFormatter(zero_direction_label=True,
                                               number_format='.0f')
    lat_formatter = LatitudeFormatter()
    ax1.xaxis.set_major_formatter(lon_formatter)
    ax1.yaxis.set_major_formatter(lat_formatter)
    cbar=plt.colorbar(sf_fill, orientation='horizontal')
    cbar.set_label('Streamfunction ($10^6$m$^2$s$^{-1}$)')
    plt.title('ITCZ approximation - %s 2018' %mons[i], fontsize=16)
    plt.savefig('plots/gcmet_ITCZ_%s_%s.png' %(hpa, months[i]))
    plt.close()
