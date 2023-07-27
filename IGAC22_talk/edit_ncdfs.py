#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#SBATCH --job-name=edit_ncdf
#SBATCH --ntasks=1
#SBATCH --mem=13gb
#SBATCH --partition=interactive
#SBATCH --time=01:25:00
#SBATCH --output=Logs/edit_ncdf.log
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
plt.style.use('seaborn-darkgrid')

years=['2015']#','2015','2010','2005','2000','1995','1990','1985','1980']
to_remove=False
replace=True; err=''
for year in years:
    print( year )
    indir=f'j25-SS_run_2015'
    outdir='j25-SS_run_SMALL_2015'
    inpath  = f'/mnt/lustre/users/mjr583/GC/13.1.2/rundirs/{indir}/OutputDir/'
    outpath = f'/mnt/lustre/users/mjr583/GC/13.1.2/rundirs/{outdir}/OutputDir/'
    want=['O3','HNO2','NO','NO2']#'NIT','NITs','NO','NO2','CO','SO2','SO4','SO4s','NH3','HNO3','C2H6','C3H8','BENZ']
    print( inpath )
    for infile in sorted(glob.glob(f'{inpath}/GEOSChem.SpeciesConc.{year}*')):
        print( infile )
        outfile = infile.replace( inpath,'' )
        # Check if it already exists
        import os.path
        if os.path.isfile(f'{outpath}/{outfile}') == True:
            if replace:
                pass
            else:
                continue

        print( outfile )        
        trimmed_outfile = Dataset(f'{outpath}/{outfile}','w')#,format='NETCDF4_CLASSIC')
        ## Read in variables from original file 
        fh   = Dataset(infile,'r')
        keys = fh.variables.keys()
        lats = fh.variables['lat']
        lons = fh.variables['lon']
        levs = fh.variables['lev']
        hyams= fh.variables['hyam']
        T    = fh.variables['time']
        try:
            if len( T[:] ) < 28:
                to_remove=True
                err='Not enough time dimensions'
        except:
            to_remove=True

        ## Build new ncdf
        time = trimmed_outfile.createDimension('time',len(T))
        lat  = trimmed_outfile.createDimension('lat',len(lats))
        lon  = trimmed_outfile.createDimension('lon',len(lons))
        lev  = trimmed_outfile.createDimension('lev',len(levs[:35]))

        time = trimmed_outfile.createVariable('time',np.int64, ('time'))
        lat  = trimmed_outfile.createVariable('lat', float, ('lat'), fill_value=np.nan)
        lon  = trimmed_outfile.createVariable('lon',float, ('lon'), fill_value=np.nan)
        lev  = trimmed_outfile.createVariable('lev',float, ('lev'), fill_value=np.nan)
        hyam = trimmed_outfile.createVariable('hyam',float, ('lev'), fill_value=np.nan)

        time.axis=T.axis
        time.units=T.units
        time.calendar=T.calendar
        time.long_name=T.long_name
        time[:]=T[:]
        
        lat.long_name=lats.long_name
        lat.units=lats.units
        lat.axis=lats.axis
        lat[:]=lats[:]

        lon.long_name=lons.long_name
        lon.units=lons.units
        lon.axis=lons.axis
        lon[:]=lons[:]

        lev.long_name=levs.long_name
        lev.units=levs.units
        lev.axis=levs.axis
        lev[:]=levs[:35]

        hyam.long_name=hyams.long_name
        hyam.units=hyams.units
        hyam[:]=hyams[0]

        for species in keys:
            if species.replace('SpeciesConc_','') in want:
                v = trimmed_outfile.createVariable( species ,float, ('time','lev','lat','lon'), fill_value=np.nan)
                v.long_name=fh.variables[species].long_name
                v.units=fh.variables[species].units
                v.averaging_method = fh.variables[species].averaging_method
                try:
                    v[:]=fh.variables[species][:,:35,:,:]
                except:
                    to_remove=True
        ## Add OH in manually
        
        try:
            ofh = Dataset( infile.replace('SpeciesConc','ConcAfterChem') )
            v = trimmed_outfile.createVariable( 'SpeciesConc_OH', float, ('time','lev','lat','lon'), fill_value=np.nan)
            v.long_name='OH from ConcAfterChem'
            v.units=ofh.variables['OHconcAfterChem'].units
            v[:] = ofh.variables['OHconcAfterChem'][:,:35,:,:]
        except:
            to_remove=True
            err='problem with ConcAfterChem'
        
        trimmed_outfile.close()
        if to_remove:
            print( 'Incomplete file, deleting' )
            os.remove(f'{outpath}/{outfile}')
        to_remove=False
