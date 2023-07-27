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
import read as R
plt.style.use('seaborn-darkgrid')

NEI, tlon, tlot = R.read_em(source='NEI', variable='ETHA')

import pandas as pd
years = range(2000, 2020, 1)
years = [ str(xx) for xx in years ]
x = pd.to_datetime( years, format="%Y" )

fh = Dataset('/mnt/lustre/groups/chem-acm-2018/earth0_data/GEOS/HEMCO/CEDS/v2021-06/2017/C2H6-em-anthro_CMIP_CEDS_2017.nc')
lat = fh.variables['lat'][:]
lon = fh.variables['lon'][:]

v21_voc = np.load('temp_files/v21_voc.npy')
v20_voc = np.load('temp_files/v20_voc.npy')
v18_voc = np.load('temp_files/v18_voc.npy')

v21_ethane = np.load('temp_files/v21_ethane.npy')
v20_ethane = np.load('temp_files/v20_ethane.npy')
v18_ethane = np.load('temp_files/v18_ethane.npy')

regions = ['Global','Europe','North America', 'South America','Asia','Africa','Oceania','UK']
regions = ['North America','UK']
emep=False
for region in regions:
    print(region)
    if region == 'UK':
        ## Read and plot Luke's EMEP scale factors
        emep = pd.read_csv('temp_files/EMEP_NMVOC_ScaleFactors_allYears.csv', index_col=1)
        emep = np.array(emep.iloc[4][11:])
    elif region == 'North America':
        nei = np.zeros((len(v21_voc))) + 0.05
        
    v21_voc_region=[] ; v21_ethane_region=[]
    for i in range(len(v21_voc)):
        v21_voc_region.append( R.regional_totals( v21_voc[i,:,:], lat, lon, region=region) )
        v21_ethane_region.append( R.regional_totals( v21_ethane[i,:,:], lat, lon, region=region) )
    v20_voc_region=[] ; v20_ethane_region=[]
    for i in range(len(v20_voc)):
        v20_voc_region.append( R.regional_totals( v20_voc[i,:,:], lat, lon, region=region) )
        v20_ethane_region.append( R.regional_totals( v20_ethane[i,:,:], lat, lon, region=region) )
    v18_voc_region=[] ; v18_ethane_region=[]
    for i in range(len(v18_voc)):
        v18_voc_region.append( R.regional_totals( v18_voc[i,:,:], lat, lon, region=region) )
        v18_ethane_region.append( R.regional_totals( v18_ethane[i,:,:], lat, lon, region=region) )
    
    v18_voc_region = np.array(v18_voc_region)
    v20_voc_region = np.array(v20_voc_region)
    v21_voc_region = np.array(v21_voc_region)

    v18_ethane_region = np.array(v18_ethane_region)
    v20_ethane_region = np.array(v20_ethane_region)
    v21_ethane_region = np.array(v21_ethane_region)

    fig,ax = plt.subplots(figsize=(9,5))
    p1, = ax.plot(  x[:len(v18_ethane_region)], v18_ethane_region / v18_voc_region, ls='-', label='CEDS v2018 (to 2014)' )
    p2, = ax.plot(  x[:len(v20_ethane_region)], v20_ethane_region / v20_voc_region, ls='--', label='CEDS v2020 (to 2017)' )
    p3, = ax.plot(  x, v21_ethane_region / v21_voc_region, ls=':', label='CEDS v2021 (to 2019)' )
    if region == 'UK':
        ax.plot( x[:len(emep)] , emep, ls='--', label='EMEP Scale Factors')
    if region == 'North America':
        ax.plot( x[:len(nei)] , nei, ls='--', label='NEI Scale Factors')

    ax.legend()
    ax.set_ylabel( 'C2H6 / NMVOCs Scale Factor'  )
    plt.savefig( 'plots/%s_ratio.png'  %region.replace(" ", "_") )
    plt.close()


