#!/usr/bin/env python3
import os
import xarray as xr
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('agg')
from netCDF4 import Dataset
import sys
sys.path.append('/users/mjr583/python_lib')
import GC_tools as GC
import RowPy as rp
import matplotlib.colors as mcolors
from matplotlib.colors import LogNorm
plt.style.use("seaborn-darkgrid")
import cartopy.crs as ccrs

rgb_WhGrYlRd = np.genfromtxt('/users/mjr583/python_lib/colormaps/WhGrYlRd.txt',
                                     delimiter=' ')
WhGrYlRd = mcolors.ListedColormap(rgb_WhGrYlRd/255.0)
cmap=WhGrYlRd

regions={
    'Global' :              { 'minlon' : -180, 'maxlon' : 180, 'minlat' : -90., 'maxlat' : 90. },
    'East Asia' :           { 'minlon' : 70, 'maxlon' : 165, 'minlat' : 0., 'maxlat' : 50. },
    'North America'   :     { 'minlon' : -160, 'maxlon' : -60, 'minlat' : 15., 'maxlat' : 72. },
    'Europe'   :            { 'minlon' : -30, 'maxlon' : 90, 'minlat' : 30., 'maxlat' : 90. },
    'NH'   :                { 'minlon' : -180, 'maxlon' : 180, 'minlat' : 0., 'maxlat' : 90. },
    'Southern Hemisphere':  { 'minlon' : -180, 'maxlon' : 180, 'minlat' : -90., 'maxlat' : 0. },
        }

def find_file_list(path, substrs):
    file_list =[]
    print( path, substrs )
    for root, directory, files in os.walk(path):
        for f in files:
            for s in substrs:
                if s in f:
                    file_list.append(os.path.join(root, f))
    file_list.sort()
    return file_list

def read_hemco_diag(rundir, year, variable='', version='13.4.0'):
    path=f'/mnt/lustre/users/mjr583/GC/{version}/rundirs/{rundir}/OutputDir/'
    flist=find_file_list(path, [f'HEMCO_diagnostics.{year}'])
    ds = xr.open_mfdataset(flist, combine='by_coords')

    lon  = ds['lon'][:]
    lat  = ds['lat'][:]
    lev  = ds['lev'][:]
    area = ds['AREA'][:]
    ds = ds[variable]
    return ds, lon, lat, lev, area

def total_by_region( ds, ds_ ):
    #ds = emission_to_tg( ds, area, sumems=False)#
    for region in regions:
        minlon,maxlon,minlat,maxlat=(regions[region]['minlon'],regions[region]['maxlon'], 
                regions[region]['minlat'], regions[region]['maxlat'])
        em = ds.\
                where(ds.lon>= minlon, drop=True).\
                where(ds.lon<= maxlon, drop=True).\
                where(ds.lat>= minlat, drop=True).\
                where(ds.lat<= maxlat, drop=True)
        reg_area = area.\
                where(area.lon>= minlon, drop=True).\
                where(area.lon<= maxlon, drop=True).\
                where(area.lat>= minlat, drop=True).\
                where(area.lat<= maxlat, drop=True)

        em = emission_to_tg( em, reg_area, sumems=False)#
        em = np.nansum( em )
        ds_.append( em )
    return ds

def emission_to_tg(ems, area, unit="kg/m2/s", sumems=False):
    hold=[]
    hold = ems.values * area.values                                       # m-2 to gridbox
    hold = hold * 86400 * np.mean([31,28,31,30,31,30,31,31,30,31,30,31])  # s-1 to yr-1
    hold = hold * 1e-9                                                    # Kg to Tg
    if sumems:
        em = np.sum( np.sum( hold , 1 ) , 1 )
        return em
    return hold

def plot(s, ds0, ds1, lon, lat, sname="HEMCO"):
    X, Y = np.meshgrid( lon, lat )
    #f, (ax1,ax2,ax3) = plt.subplots(1,3, figsize=(12,4))
    f=plt.figure(figsize=(11,4))
    ax1 = f.add_subplot( 131,projection=ccrs.Robinson(), aspect='auto')
    ax1.coastlines(lw=.3)
    ax2 = f.add_subplot( 132,projection=ccrs.Robinson(), aspect='auto')
    ax2.coastlines(lw=.3)
    ax3 = f.add_subplot( 133,projection=ccrs.Robinson(), aspect='auto')
    ax3.coastlines(lw=.3)

    #m1 = rp.get_basemap(ax=ax1)
    #m2 = rp.get_basemap(ax=ax2)
    #m3 = rp.get_basemap(ax=ax3)
    diff = - ( 100 - ( np.squeeze(ds1) / np.squeeze(ds0) * 100 ) )
    #print( diff )
    #print( type( diff ) )
    #sys.exit()
    #diff = np.squeeze(ds0) - np.squeeze(ds1)
    MAX=np.nanmax(ds1)
    norm=LogNorm( vmax=MAX )
    #im1=m1.pcolormesh( X, Y, np.squeeze(ds0),cmap=cmap, norm=norm, label="Control")
    #im2=m2.pcolormesh( X, Y, np.squeeze(ds1),cmap=cmap, norm=norm, label="Scaled")
    #im3=m3.pcolormesh( X, Y, diff,cmap='bwr', center=0)#, vmin=-1000, vmax=1000, label="Difference")
    im1=np.squeeze(ds0).plot.imshow( x='lon',y='lat', ax=ax1, transform=ccrs.PlateCarree(), vmax=MAX, cmap=cmap, add_colorbar=False)
    im2=np.squeeze(ds1).plot.imshow( x='lon',y='lat', ax=ax2, transform=ccrs.PlateCarree(), vmax=MAX, cmap=cmap, add_colorbar=False)
    im3=diff.plot.imshow( x='lon',y='lat', ax=ax3, transform=ccrs.PlateCarree(), center=0., cmap='bwr', add_colorbar=False)
    
    ax1.set_title('CEDS')
    ax2.set_title('Scaled CEDS')
    ax3.set_title('Scaled - CEDS')
    plt.subplots_adjust(top=1.5, bottom=0.1)
    cbar_ax = f.add_axes([0.046, 0.15, 0.602, 0.05])
    f.colorbar(im1, cax=cbar_ax, orientation='horizontal')
    cbar_ax.set_xlabel(f'{s} (Tg)')
    cbar_ax = f.add_axes([0.675,0.15,0.3, 0.05])
    f.colorbar(im3, cax=cbar_ax, orientation='horizontal')
    cbar_ax.set_xlabel('%')
    plt.tight_layout()
    plt.savefig(f'plots/{sname}.png')
    plt.close()

def main():
    global area, regions
    rundirs=['ceds_only','all_scaled']
    species=['C2H6','C3H8','BENZ','EOH','CH2O','MEK','ALK4','TOLU','XYLE','PRPE']
    df_=[]
    for s in species:
        version='13.4.0' 
        ds0, lon, lat, lev, area = read_hemco_diag(rundirs[0], year='20150101', variable=f'Emis{s}_Anthro')
        ds1, lon, lat, lev, area = read_hemco_diag(rundirs[1], year='20150101', variable=f'Emis{s}_Anthro')
        
        ds0 = ds0.sum(dim='lev', keep_attrs=True)
        ds1 = ds1.sum(dim='lev', keep_attrs=True)
         
        ds0_=[]
        ds0 = total_by_region( ds0, ds0_ )

        ds1_=[]
        ds1 = total_by_region(ds1, ds1_)
        
        df = pd.DataFrame( {f'CEDS {s}':ds0_, f'Scaled {s}': ds1_}, index=regions)
        df = df.transpose()
        plot(s, ds0, ds1, lon, lat, sname=f"TEST.13.4.0_scale_{s}")
        df_.append(df)
    #df = pd.concat(df_, axis=0)
    #print( df )#.transpose() )
    #df.to_csv('ceds-and-scaled_emissions_vocs.csv')
if __name__ == "__main__":
    main()
