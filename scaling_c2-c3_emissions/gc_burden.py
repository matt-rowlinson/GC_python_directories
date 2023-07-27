#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#SBATCH --job-name=gc_burde
#SBATCH --ntasks=1
#SBATCH --mem=2Gb
#SBATCH --partition=interactive
#SBATCH --time=00:05:00
#SBATCH --output=Logs/gc_burden_%A.log
import sys
import pandas as pd
import numpy as np
import glob
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib
matplotlib.use('agg')
sys.path.append('/users/mjr583/python_lib')
import GC_tools as GC
import RowPy as rp
from CVAO_dict import CVAO_dict as d
import CVAO_tools as CV
plt.style.use('seaborn-darkgrid')

### Set constants ###
zboltz=1.3807e-23       #Boltzmann constant
Rd=287.05               #gas constant
avc=6.022e23            #Avogadro's constant
mm_da=avc*zboltz/Rd     #molar mass of dry air (approx 28 g/mol)

mm = {
        'O3'   : 0.0479982,
        'moh'  : 0.03204,
        'CO'   : 0.02801,
        'CH4'  : 0.01604,
        'C2H6' : 0.03007
        }


def get_tropopause_ps(rundir, version, year):
    for i,infile in enumerate(sorted(glob.glob('/users/mjr583/scratch/GC/%s/rundirs/%s/OutputDir/GEOSChem.StateMet.*%s*.nc4' % ('13.1.2', rundir, year)))):
        print(infile)
        fh=Dataset(infile)
        Tlev=fh.variables['Met_TropLev'][:]
        ad=fh.variables['Met_AD'][:]
        return Tlev, ad

def trop_only(var, lon, lat, lev):
    for ilon in range(len(lon)):
        #print(ilon)
        for ilat in range(len(lat)):
            for ilev in range(1,73):            
                if (ilev > T[ilat,ilon]):
                    var[:,ilev-1,ilat,ilon] = np.nan
    return var

def calc_burden(VAR, variable):
    mm=mm[variable] ; burden_=[]
    print( mm )
    var=np.zeros(72) * np.nan
    for v in range(len(VAR)):
        for i in range(72):
            var[i] = np.nansum( VAR[v,i,:,:] * AD[i,:,:] * mm/mm_da) * 1e-9
            
        burd = np.round(np.nansum(var),2)
        burden_.append(burd)
    print(f'Tropospheric {variable} burden is {np.mean(burden)} Tg')
    return burden_


def plot():
    plt.plot(time, burden_, label=labels[r])
    plt.ylabel(variable+' (Tg)')
    plt.legend()
    plt.savefig('plots/ethane_burden.png')
    plt.close()


def main():

    variable='O3'
    rundirs=['full_nitrate_photol_control']
    labels =['Base']
    year='2017'

    for r,rundir in enumerate(rundirs):
        var, lat,lon,lev,time = GC.get_gc_var(rundir=rundir, variable=variable, version='13.1.2',year='2017')
        VAR= var / float(d[variable]['scale'])

    T, AD = get_tropopause_ps(rundir='nitrate_photol_control', version='13.1.2', year='2016')
    T = T[0]
    AD= AD[0]
    print( VAR.shape )

    VAR = trop_only(VAR, lon, lat, lev)
    burden = calc_burden(VAR, variable)

    #plot(time, burden_, label=labels[r])

if __name__=="__main__":
    main()
