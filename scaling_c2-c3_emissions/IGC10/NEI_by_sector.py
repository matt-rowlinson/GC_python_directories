#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#SBATCH --job-name=timeseries
#SBATCH --ntasks=1
#SBATCH --mem=10gb
#SBATCH --partition=interactive
#SBATCH --time=00:10:00
#SBATCH --output=Logs/timeseries_%A.log
import sys
sys.path.append('/users/mjr583/python_lib')
import RowPy as rp
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib
import read as R
import re
import glob
from matplotlib.colors import LogNorm
import numpy as np
from netCDF4 import Dataset
matplotlib.use('agg')
plt.style.use('seaborn-darkgrid')
#sae = np.swapaxes(rp.surface_area_earth(3600,1800 , verbose=False),0,1)

import pickle
f = open('/users/mjr583/NEI/nei_sae.pkl', 'rb')
sae = pickle.load(f)
f.close()

def read_emissions(species, sectors=False, ignore_sector=False, sae=sae):
    non_voc=['CO','NO','NO2','SO2','SO4','BC','OC','NH3','HNO2']
    inpath = '/mnt/lustre/groups/chem-acm-2018/earth0_data/GEOS/HEMCO/NEI2011/v2017-02-MM/*/NEI11_0.1x0.1_2011*.nc'
    em=0

    for infile in sorted(glob.glob('%s' %inpath)):
        fh = Dataset(infile)
        keys=list(fh.variables.keys())

        for key in keys:
            try:
                nei11_name = fh.variables[key].NEI11_name
            except:
                continue
            if nei11_name in non_voc:
                continue
            #print(nei11_name)
            if 'lev' in keys:
                #print('Aqui')
                em +=  np.sum( fh.variables[key][:] , 1 )
            else:
                em +=  fh.variables[key][:]
    lat = fh.variables['lat'][:]
    lon = fh.variables['lon'][:]

    nei_lats = np.arange( -89.95 , 90.05 , 0.1 )
    nei_lons = np.arange( -179.95 , 180.05 , 0.1 )
   
    lat_fltr = np.where(( nei_lats >= lat.min() ) & ( nei_lats <= lat.max() ))
    lon_fltr = np.where(( nei_lons >= lon.min()-0.1 ) & ( nei_lons <= lon.max() ))
    
    sae = sae[lat_fltr[0][0]:lat_fltr[0][-1]+1, lon_fltr[0][0]:lon_fltr[0][-1]+1]

    em = np.sum(em, 0)
    em = em * sae * (3600*24)*365/12 * 1e-6
    #ann_em.append(em)

    time=fh.variables['time']
    t = re.findall(r'\b\d+\b', time.units)
    t = [ int(x) for x in t ] 
    time = rp.days_to_datetime(time, start_year=t[0], start_month=t[1], start_day=t[2], start_hour=t[3])

    lat = fh.variables['lat'][:]
    lon = fh.variables['lon'][:]
    if not sectors and not ignore_sector:
        return em, lon, lat, time
    em = sector_select(fh, sectors=sectors, ignore_sector=ignore_sector)
    return em, lon, lat, time

def sector_select(fh, sectors=False, ignore_sector=False):
    if sectors=='shp':
        sectors='shipping'
    if ignore_sector=='shp':
        ignore_sector='shipping'
    if not sectors and not ignore_sector:
        return fh.variables['sumallsectors'][:]
    
    if sectors:
        if type(sectors) == str:
            return fh.variables[sectors][:]
        ems=[]
        for sector in sectors:
            ems.append( fh.variables[sector][:] )
        em = np.sum( ems , axis=0 )
        return em
    
    if ignore_sector:
        exclude = fh.variables[ignore_sector][:]
        return fh.variables['sumallsectors'][:] - exclude
 

def get_region(em, lon, lat, lonlat=(-180, 180, -90, 90), only_land=False):
    from mpl_toolkits.basemap import Basemap
    m = Basemap(projection='cyl', llcrnrlat=-90,urcrnrlat=90,llcrnrlon=-180, urcrnrlon=180.,\
            resolution='c',area_thresh=1000.)
    print( np.nansum( em ) )
    
    for ilon in range(len(lon)):
        if lonlat[0] < lon[ilon] < lonlat[1]:
            for ilat in range(len(lat)):
                if lonlat[2] < lat[ilat] < lonlat[3] :
                    if only_land:
                        if m.is_land(lon[ilon],lat[ilat]) == False:
                            em[ ilat, ilon] = np.nan
                        else:
                            pass
                    else:
                        pass                        
                else:
                    em[ ilat, ilon] = np.nan
        else:
            em[ :, ilon] = np.nan
    print( np.nansum( em ) )

    return em

def select_years(em, time, start_year, end_year):
    if start_year:
        for n, t in enumerate(time):
            if float(t.year) == float(start_year):
                em = em[n:]
                time=time[n:]
                break
    if end_year:
        for n, t in enumerate(time):
            if float(t.year) == float(end_year):
                em = em[:n+1]
                time=time[:n+1]
                break
    return em, time

def plot_map(Z,X,Y, year='', only_land=False):
    y = str(year.year)
    f,(ax) =plt.subplots(1,1)
    m=rp.get_basemap(lllat=Y.min()-5, urlat=Y.max()+5, lllon=X.min()-5, urlon=X.max()+5, ax=ax)
    X,Y=np.meshgrid(X,Y)
    m.pcolormesh(X,Y, Z, cmap='viridis', norm=LogNorm())
    
    f.text(0.16, .685, f'NEI {y}: %s Gg' %np.round(np.nansum(Z),3),
                                    fontsize=14, bbox=dict(facecolor='w', alpha=1.) )

    plt.savefig('plots/NEI_land_%s' %y)
    plt.close()
    return


def get(species, region=False, lonlat=False, start_year=False, end_year=False, sectors=False,
                                         ignore_sector=False, only_land=False, make_map=False):
    em, lon, lat, time = read_emissions(species, sectors=sectors, ignore_sector=ignore_sector)
    if start_year or end_year:
        em, time = select_years(em, time, start_year=start_year, end_year=end_year)
    
    if region:
        em = get_region(em, lon, lat, lonlat=lonlat, only_land=only_land)
    
    if make_map:
        for n, y in enumerate(time):
            plot_map(em, lon, lat, year=y)

    return em, time, lon, lat

if __name__ == "__main__":

    e1, T, ln, lt = get(species='NMVOC', region='All_Eur', lonlat=(-139.95, -50.05, 20.05, 59.95), sectors=False, 
                                                    ignore_sector=False, only_land=True, make_map=True)

    print( np.nansum(e1) )
