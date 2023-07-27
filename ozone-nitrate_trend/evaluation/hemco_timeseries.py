#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#SBATCH --job-name=hemco_timeseries
#SBATCH --ntasks=1
#SBATCH --mem=18gb
#SBATCH --partition=interactive
#SBATCH --time=01:30:00
#SBATCH --output=Logs/hemco_timeseries_%A.log
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('agg')
import os
import numpy as np
import sys
sys.path.append('/users/mjr583/python_lib/')
import RowPy as rp


from itertools import chain
def longlist2array(longlist):
    flat = np.fromiter(chain.from_iterable(longlist), np.array(longlist[0][0]).dtype, -1) # Without intermediate list:)
    return flat.reshape((len(longlist), -1))

def find_file_list(path, substrs):
    file_list =[]
    for root, directory, files in os.walk(path):
        for f in files:
            for s in substrs:
                if s in f:
                    file_list.append(os.path.join(root, f))
    file_list.sort()
    return file_list

days=np.array([31,28,31,30,31,30,31,31,30,31,30,31])
regions={
    False : False,
    'Asia' : { 'minlon' : 70, 'maxlon' : 165, 'minlat' : 0., 'maxlat' : 50. },
    'US'   : { 'minlon' : -160, 'maxlon' : -60, 'minlat' : 15., 'maxlat' : 72. },
    'EU'   : { 'minlon' : -20, 'maxlon' : 40, 'minlat' : 36., 'maxlat' : 72. },
    'NH'   : { 'minlon' : -180, 'maxlon' : 180, 'minlat' : 0., 'maxlat' : 90. },
    'SH'   : { 'minlon' : -180, 'maxlon' : 180, 'minlat' : -90., 'maxlat' : 0. },
        }

emissions = ['HcoLightningFlashRate_Total','EmisNH3_Total', 'EmisHNO3_Ship', 'EmisNO_Total', 'EmisCO_Total','EmisC2H6_Total', 'EmisSO2_Total']
rundirs=['nitrate_photol_control']
for region in regions:
    for rundir in rundirs:
        path=f'/users/mjr583/scratch/GC/13.1.2/rundirs/{rundir}/OutputDir/'
        flist=find_file_list(path, ['HEMCO_diagnostics.201'])
        a = xr.open_mfdataset(flist)
        print( "Got files" )

        #if region:
        #    out=region
        #    minlon,maxlon,minlat,maxlat=(regions[region]['minlon'],regions[region]['maxlon'], 
        #                                 regions[region]['minlat'], regions[region]['maxlat'])
        #    a = a.\
        #            where(a.lon>= minlon, drop=True).\
        #            where(a.lon<= maxlon, drop=True).\
        #            where(a.lat>= minlat, drop=True).\
        #            where(a.lat<= maxlat, drop=True)
        #else:
        #    out='global'
        
        AREA=a.AREA
        
        EMSET=a.sum(dim='lev', keep_attrs=True)
        time = EMSET['time']
        for emis in emissions:
            em = EMSET[emis]
            em=np.array(em)
            new_em=[]
            for n in range(len(em)):
                x = em[n] * AREA[n] * ( 3600 * 24 * days.mean()) * 1e-9
                new_em.append(x.values)
            em = np.array(new_em)
            
            a = np.sum( np.sum(em, 1), 1)
            plt.plot( time,  a )
            plt.legend()
            plt.savefig(f'plots/{emis}_anthro_ethane_emissions.png')
            plt.close()
    sys.exit()
