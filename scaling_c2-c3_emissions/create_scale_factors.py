#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#SBATCH --job-name=timeseries
#SBATCH --ntasks=1
#SBATCH --mem=10gb
#SBATCH --partition=interactive
#SBATCH --time=00:10:00
#SBATCH --output=LOGS/timeseries.log
import sys
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
from netCDF4 import Dataset
matplotlib.use('agg')
sys.path.append('/users/mjr583/python_lib')
import RowPy as rp
from matplotlib.colors import LogNorm
import copy

region_dict={
    False : False,
    'Asia' : { 'minlon' : 70, 'maxlon' : 165, 'minlat' : 0., 'maxlat' : 50. },
    'NEI'   : { 'minlon' : -160, 'maxlon' : -60, 'minlat' : 15., 'maxlat' : 72. },
    'EMEP'   : { 'minlon' : -30, 'maxlon' : 90, 'minlat' : 30., 'maxlat' : 90. },
        }

species_SFs={
        'C2H6' : {'Asia' : 2.50, 'NEI' : 6.21, 'EMEP' : 2.22 },
        'C3H8' : {'Asia' : 3.50, 'NEI' : 3.36, 'EMEP' : 1.95} 
            }

def regional(em, lat ,lon, region, sf):
    reg_em = copy.deepcopy(em)
    d = region_dict
    for ilon in range(len(lon)):
        for ilat in range(len(lat)):
            if d[region]['minlon'] < lon[ilon] <  d[region]['maxlon']:
                if  d[region]['minlat'] < lat[ilat] <  d[region]['maxlat']:
                    reg_em[ilat, ilon] = sf
                else:
                    pass                    
            else:
                pass
    return reg_em

def plot_sf(grid, lons, lats, species='C2H6'):
    X, Y = np.meshgrid( lons, lats )
    f,ax=plt.subplots(1,1)
    m = rp.get_basemap(ax=ax)
    m.pcolormesh( X, Y, grid)

    plt.savefig(f'{species}_scale_factors.png')
    plt.close()
    return

def save_to_ncdf(grid, lons, lats, species='C2H6'):
    sf_outfile = Dataset(f'./CEDS_{species}-scale_factors.nc', 'w')

    lat  =  sf_outfile.createDimension('lat',len(lats))
    lon  = sf_outfile.createDimension('lon',len(lons))
    lat  =  sf_outfile.createVariable('lat', float, ('lat'))
    lon  = sf_outfile.createVariable('lon',float, ('lon'))

    lat.standard_name='latitude'
    lat.long_name='Latitude'
    lat.units='degrees_north'
    lat.axis='Y'
    lat[:]=lats

    lon.standard_name='longitude'
    lon.long_name='Longitude'
    lon.units='degrees_east'
    lon.axis='X'
    lon[:]=lons

    sf  = sf_outfile.createVariable(f'SF_{species}',float, ('lat','lon'), fill_value=np.nan)
    sf.standard_name=f'Global {species} scale factors'
    sf.long_name=sf.standard_name
    sf.units='unitless'
    sf[:]=grid
    sf_outfile.close()
    return

if __name__=="__main__":
    species='C2H6'
    nlat=180
    nlon=360
    lats = np.arange(-89.95,90.00, 180/nlat)
    lons = np.arange(-180, 180, 360/nlon)
    
    lowest=min( species_SFs[species], key=species_SFs[species].get)

    grid = np.zeros((nlat, nlon))
    grid = grid + species_SFs[species][lowest]

    grid = regional( grid, lats, lons, region='EMEP', sf=species_SFs[species]['EMEP'])
    grid = regional( grid, lats, lons, region='NEI', sf=species_SFs[species]['NEI'])
    grid = regional( grid, lats, lons, region='Asia', sf=species_SFs[species]['Asia'])
    
    plot_sf(grid, lons, lats, species=species)
    save_to_ncdf(grid, lons, lats, species=species)
