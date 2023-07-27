#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#SBATCH --job-name=timeseries
#SBATCH --ntasks=1
#SBATCH --mem=10gb
#SBATCH --partition=interactive
#SBATCH --time=00:10:00
#SBATCH --output=LOGS/timeseries.log
import os
import glob
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import sys
sys.path.append('/users/mjr583/python_lib')
import RowPy as rp

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
    
    sae = rp.surface_area_earth(144,91)
    sae=np.swapaxes(sae,0,1)
    days=np.array([31,28,31,30,31,30,31,31,30,31,30,31])

    new=[]
    for n in range(12):
        x = tzompa[n] * sae * ( 3600 * 24 * days[n]) * 1e-9
        new.append(x)
    tzompa = np.array(new)
    print('Tzompa:', np.round(np.nansum(tzompa),2), 'Tg yr-1')
    
    return tzompa, alon, alat


def xiao():
    path='/mnt/lustre/groups/chem-acm-2018/earth0_data/GEOS/ExtData/HEMCO//XIAO/v2014-09/'
    flist=find_file_list(path, ['1x1.nc'])
    print(flist)

    hold=np.zeros((1,181,360))
    for n, infile in enumerate(flist):
        fh=Dataset(infile)
        total = fh.variables['C3H8'][:]
        print(total.shape)
        
        
        hold[n] = total

        alat=fh.variables['lat'][:]
        alon=fh.variables['lon'][:]
    xiao=np.sum(hold,0)
    print(xiao.shape)
    factor = 44.1 / (3*12.011)
    xiao = xiao * factor
    
    sae = rp.surface_area_earth(360,181)
    sae=np.swapaxes(sae,0,1)
    days=np.array([31,28,31,30,31,30,31,31,30,31,30,31])

    new=[]
    xiao = xiao * sae *  3600 * 24 * 365 * 1e-9
    print('Xiao:', np.round(np.nansum(xiao),2), 'Tg yr-1')
    
    return xiao, alon, alat

def bf():
    path='/mnt/lustre/groups/chem-acm-2018/earth0_data/GEOS/ExtData/HEMCO/BIOFUEL/v2019-08/'
    flist=find_file_list(path, ['2x25'])

    hold=np.zeros((12,91,144))
    for n, infile in enumerate(flist):
        print(infile)
        fh=Dataset(infile)
        total = fh.variables['BIOFUEL_C3H8'][:]
        
        
        hold[n] = total

        alat=fh.variables['lat'][:]
        alon=fh.variables['lon'][:]
    bf=hold
    factor = 44.1 / (3*12.011)
    bf = bf * factor
    
    sae = rp.surface_area_earth(144,91)
    sae=np.swapaxes(sae,0,1)
    days=np.array([31,28,31,30,31,30,31,31,30,31,30,31])

    new=[]
    for n in range(12):
        x = bf[n] * sae * ( 3600 * 24 * days[n]) * 1e-9
        new.append(x)
    bf = np.array(new)
    print('Biofuel:', np.round(np.nansum(bf),2), 'Tg yr-1')
    
    return bf, alon, alat

def plot_map(Z,X,Y, year='2017',lonlat=False, out='plot.png'):
    y = str(year)
    f,(ax) =plt.subplots(1,1)
    if lonlat:
        m=rp.get_basemap(lllon=lonlat[0], urlon=lonlat[1], lllat=lonlat[2], urlat=lonlat[3], ax=ax)
    else:
        m=rp.get_basemap(lllat=Y.min(), urlat=Y.max(), lllon=X.min(), urlon=X.max(), ax=ax)
    X,Y=np.meshgrid(X,Y)
    m.pcolormesh(X,Y, Z, cmap='viridis', norm=LogNorm())
    
    f.text(0.16, .685, f'NEI {y}: %s Gg' %np.round(np.nansum(Z),3),
                                    fontsize=14, bbox=dict(facecolor='w', alpha=1.) )

    plt.savefig('plots/%s' %out)
    plt.close()
    return


def get_region(em, lon, lat, lonlat=(-180, 180, -90, 90), only_land=False):
    from mpl_toolkits.basemap import Basemap
    m = Basemap(projection='cyl', llcrnrlat=-90,urcrnrlat=90,llcrnrlon=-180, urcrnrlon=180.,\
            resolution='c',area_thresh=1000.)
    
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

    return em

def get( inv, lonlat=False, only_land=False):
    if inv=='tzompa':
        em, lon, lat = tzompa()
    elif inv=='xiao':
        em, lon, lat = xiao()
    elif inv=='bf':
        em, lon, lat = bf()

    if len(em.shape) >= 3:
        em = np.nansum(em,0)
    em =  get_region(em, lon, lat, lonlat=lonlat, only_land=only_land)
    em = em * 1e3
    return em, lon, lat



if __name__ == "__main__":
    
    em, ln, lt = get('xiao', lonlat=(-139.95, -50.05, 20.05, 59.95), only_land=True)
    print( np.nansum(em))
 
    plot_map( em , ln, lt, out='xiao_ propane.png')

