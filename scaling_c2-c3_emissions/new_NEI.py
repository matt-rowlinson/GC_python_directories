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
import re
from netCDF4 import Dataset
matplotlib.use('agg')
sys.path.append('/users/mjr583/python_lib')
import GC_tools as GC
import RowPy as rp
from CVAO_dict import CVAO_dict as d
import CVAO_tools as CV
#sae = np.swapaxes(rp.surface_area_earth(3600,1800 , verbose=False),0,1)

import pickle
f = open('/users/mjr583/NEI/nei_sae.pkl', 'rb')
sae = pickle.load(f)
print(sae.shape)
f.close()

non_voc=['lev', 'CO','NO','NO2','NOX','SO2','SO4','CH4','BC','OC','NH3','HNO2', 'CH4_INV','CO2_INV','N2O_INV']

def read_emissions(years , sae=sae ):
    ann_sum=[] ; ann_em=[]
    for y in years:
        print(y)
        inpath = '/mnt/lustre/users/mjr583/GC/emissions/NEI/*/*.ncf'
        em=0
        for infile in sorted(glob.glob('%s' %inpath)):
            
            if "_3D_" in infile:
                continue                
            if "_merge_0" not in infile:
                continue

            fh = Dataset(infile)
            keys=list(fh.variables.keys())
            for key in keys[3:]:
                if key in non_voc:
                    continue
                if key != "VOC_INV":
                    continue

                if 'lev' in keys:
                    #print(key)
                    em +=  np.sum( fh.variables[key][:] , 1 )
                else:
                    em +=  fh.variables[key][:]
            #sys.exit()
        lat = fh.variables['lat'][:]
        lon = fh.variables['lon'][:]

        nei_lats = np.arange( -89.95 , 90.05 , 0.1 )
        nei_lons = np.arange( -179.95 , 180.05 , 0.1 )
       
        lat_fltr = np.where(( nei_lats >= lat.min() ) & ( nei_lats <= lat.max() ))
        lon_fltr = np.where(( nei_lons >= lon.min()-0.1 ) & ( nei_lons <= lon.max() ))
        
        sae = sae[lat_fltr[0][0]:lat_fltr[0][-1]+1, lon_fltr[0][0]:lon_fltr[0][-1]+1]
        em = np.sum(em, 0)
        em = em * sae * (3600*24)*365/12 * 1e-6

        ann_em.append(em)
    emissions=np.array(ann_em)
    em=emissions
    return emissions, lon, lat


def read_species(years , sae=sae, species=False):
    ann_sum=[] ; ann_em=[]
    for y in years:
        inpath = '/mnt/lustre/users/mjr583/GC/emissions/NEI/*/*.ncf'
        em=0
        for infile in sorted(glob.glob('%s' %inpath)):
            
            if "_3D_" in infile:
                continue                
            if "_merge_0" not in infile:
                continue

            fh = Dataset(infile)
            keys=list(fh.variables.keys())
            if species in keys:
                if 'lev' in keys:
                    em +=  np.sum( fh.variables[species][:] , 1 )
                else:
                    em +=  fh.variables[species][:] 
        lat = fh.variables['lat'][:]
        lon = fh.variables['lon'][:]

        nei_lats = np.arange( -89.95 , 90.05 , 0.1 )
        nei_lons = np.arange( -179.95 , 180.05 , 0.1 )
       
        lat_fltr = np.where(( nei_lats >= lat.min() ) & ( nei_lats <= lat.max() ))
        lon_fltr = np.where(( nei_lons >= lon.min()-0.1 ) & ( nei_lons <= lon.max() ))
        
        sae = sae[lat_fltr[0][0]:lat_fltr[0][-1]+1, lon_fltr[0][0]:lon_fltr[0][-1]+1]

        em = np.sum(em, 0)
        em = em * sae * (3600*24)*365/12 * 1e-6
        ann_em.append(em)
    emissions=np.array(ann_em)
    return emissions, lon, lat




def get_region(em, lon, lat, is_land=True, lonlat=(-180, 180, -90, 90)):
    from mpl_toolkits.basemap import Basemap
    m = Basemap(projection='cyl', llcrnrlat=-90,urcrnrlat=90,llcrnrlon=-180, urcrnrlon=180.,\
            resolution='c',area_thresh=1000.)

    for ilon in range(len(lon)):
        if lonlat[0] < lon[ilon] < lonlat[1]:
            for ilat in range(len(lat)):
                if lonlat[2] < lat[ilat] < lonlat[3] :
                    if is_land:
                        if m.is_land(lon[ilon],lat[ilat]) == True:
                            pass
                        else:
                            #pass
                            em[:, ilat, ilon] = np.nan
                    else:
                        pass
                else:
                    em[:, ilat, ilon] = np.nan
        else:
            em[:, :, ilon] = np.nan
    return em

def year_list(start_year, end_year=None):
    end_year = start_year if end_year is None else end_year
    x = range(int(start_year), int(end_year)+1, 1)
    years = [ str(xx) for xx in x ]
    return years

def get(start_year, end_year=None, is_land=True, region=False, lonlat=False):
    years = year_list(start_year, end_year)
    em, lon, lat = read_emissions(years )
    em = get_region(em, lon, lat, is_land=is_land, lonlat=lonlat)
    for n, y in enumerate(years):
        print(np.nansum(em[n]))

    return em, years

def get_species(start_year, end_year=None, is_land=True, region=False, lonlat=False, species=False):
    end_year = start_year if end_year is None else end_year
    years = year_list(start_year, end_year)
    em, lon, lat = read_species(years, species=species)
    em = get_region(em, lon, lat, is_land=is_land, lonlat=lonlat)
    for n, y in enumerate(years):
        print(np.nansum(em[n]))

    return em, years

def main():
    e, years = get(start_year=2017, end_year=2017, lonlat=(-9.55, 18.6, 35.95, 61.65))

    print( e )

if __name__ == "__main__":
    main()

    #get(start_year=2012, end_year=2017, region='UK', lonlat=(-8.55, 2.75, 49.15, 61.65)) 
    #get(start_year=2012, end_year=2017, region='Iberia', lonlat=(-9.55, 3.38, 35.95, 43.28))
    #get(start_year=2012, end_year=2017, region='W_Europe', lonlat=(-9.55, 18.6, 35.95, 61.65))
    #get(start_year=2017, end_year=2017, region='North_America', is_land=False, lonlat=(-139.95, -50.05, 20.05, 59.95) )
    #sys.exit()
    #get_species(start_year=2017, end_year=2017, region='North_America', is_land=False, lonlat=(-139.95, -50.05, 20.05, 59.95), 
    #                            species='ETHA')
