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
import matplotlib
import numpy as np
from netCDF4 import Dataset
matplotlib.use('agg')
sys.path.append('/users/mjr583/python_lib')
import GC_tools as GC
import RowPy as rp
from CVAO_dict import CVAO_dict as d
import CVAO_tools as CV
from matplotlib.colors import LogNorm

csv = '/mnt/lustre/users/mjr583/GC/emissions/Geo-CH4_emission_grid_files/Gridded_geoCH4_csv_ESSD/Total_geoCH4_output_2018.csv'
d = pd.read_csv(csv)
d.columns=['lon','lat','emission']
d=d.replace(-9999.0, np.nan)

em = d['emission']
em = np.reshape(em.values, (180, 360), order='C')

em = em * 1e3 # Tonnes to kg
em = em / (3600 * 24 * 365) # yr-1 to s-1
sae = np.swapaxes(rp.surface_area_earth(360,180),0,1)
em = em / sae # gridcell-1 to m-2

lons = np.arange(-179.5, 180.5, 1)
lats = np.arange(-89.5, 90.5, 1)
X,Y = np.meshgrid(lons,lats)

fig, ax = plt.subplots(figsize=(9,6))
cax = fig.add_axes([0.15, 0.1, 0.7, 0.05])

m = rp.get_basemap(ax=ax)
im = m.pcolormesh(X,Y, em, norm=LogNorm())
cbar = fig.colorbar(im, cax=cax, orientation='horizontal')
cbar.set_label('kg m-2 s-1')

x = em * sae * ( 3600 * 24 * 365) * 1e-9
plt.text(.1, 1.1, 'Global Geological $CH_4$ emission = %s Tg yr$^{-1}$' %np.round(np.nansum(x),2) , 
            fontsize=16, transform=ax.transAxes)
plt.savefig('plots/GeoCH4_gridded_ESSD_Etiope2019.png')
plt.close()

print('Global Geological CH4 emissions:', np.round(np.nansum(x),2), 'Tg yr-1')

## Scale to C2H6 emissions and plot
em = em * ( 3 / np.nansum(x))
fig, ax = plt.subplots(figsize=(9,6))
cax = fig.add_axes([0.15, 0.1, 0.7, 0.05])

m = rp.get_basemap(ax=ax)
im = m.pcolormesh(X,Y, em, norm=LogNorm())
cbar = fig.colorbar(im, cax=cax, orientation='horizontal')
cbar.set_label('kg m-2 s-1')

x = em * sae * ( 3600 * 24 * 365) * 1e-9
plt.text(.1, 1.1, 'Global Geological $C_2H_6$ emission = %s Tg yr$^{-1}$' %np.round(np.nansum(x),2), 
            fontsize=16, transform=ax.transAxes)
plt.savefig('plots/Scaled_C2H6_Geo_emissions.png')
plt.close()

print('Scaled Global Geological C2H6 emissions:', np.round(np.nansum(x),2), 'Tg yr-1')


monthly_em = np.zeros((12,180,360))
days=[31,28,31,30,31,30,31,31,30,31,30,31]
for n in range(12):
    monthly_em[n] = em #/ 365 * days[n]

## Write to a netcdf suitable for HEMCO
mons=['01','02','03','04','05','06','07','08','09','10','11','12']
for i,ii in enumerate(mons):
    #print( np.nansum(monthly_em[i]) *  )

    emissions_outfile = Dataset('C2H6-em-Geo_%s.nc' %ii , 'w')

    time = emissions_outfile.createDimension('time',None)
    lat  =  emissions_outfile.createDimension('lat',180)
    lon  = emissions_outfile.createDimension('lon',360)

    time = emissions_outfile.createVariable('time',float, ('time'))
    lat  =  emissions_outfile.createVariable('lat', float, ('lat'))
    lon  = emissions_outfile.createVariable('lon',float, ('lon'))

    import time as T
    import datetime
    new_time=[]
    s = ["01/01/2010","01/02/2010","01/03/2010","01/04/2010","01/05/2010","01/06/2010","01/07/2010","01/08/2010",
            "01/09/2010","01/10/2010","01/11/2010","01/12/2010"]
    dt = datetime.datetime.strptime(s[i], '%d/%m/%Y')
    stamp = ( dt - datetime.datetime(1950,1,1)).days
    
    time.standard_name='time'
    time.long_name='Time'
    time.axis='T'
    time.units='days since 1950-01-01 00:00:00'
    time.calendar='standard'
    time[:]=stamp

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

    new_em  = emissions_outfile.createVariable('Geo_C2H6',float, ('lat','lon'), fill_value=np.nan)
    new_em.standard_name='Geological_C2H6'
    new_em.long_name='Natural C2H6 emissions from geological sources'
    new_em.units='kg m-2 s-1'
    print(monthly_em.shape)
    new_em[:]=monthly_em[i]

    emissions_outfile.description = 'Data from Etiope'
    import datetime
    emissions_outfile.date_created =  str(datetime.datetime.now())
    emissions_outfile.created_by = "Matthew Rowlinson (matthew.rowlinson@york.ac.uk)"
    emissions_outfile.grid = "1.0 x 1.0 degree latitude x longitude grid"
    emissions_outfile.description  = "File contains C2H6 natural geological emissions produced for use in \
    GEOS-Chem chemical transport model.C2H6 geological emission distribution is based on the CH4 geological \
    emissions grid from Etiope et al. 2019. Global emissions are scaled to 3 Tg yr-1, in accordance with \
    Dalsoren et al. 2018. File has no monthly or interannual variation, based on 'present-day' estimate annual emissions."
    emissions_outfile.references = "Dalsøren, S.B., Myhre, G., Hodnebrog,O. et al. Discrepancy between simulated and observed ethane and propane levels explained by underestimated fossil emissions. Nature Geosci 11, 178-184 (2018). https://doi.org/10.1038/s41561-018-0073-0 ; Etiope, Giuseppe & Ciotoli, Giancarlo & Schwietzke, Stefan & Schoell, Martin. (2019). Gridded maps of geological methane emissions and their isotopic signature. Earth System Science Data. 11. 1-22. 10.5194/essd-11-1-2019. "

    emissions_outfile.close()

