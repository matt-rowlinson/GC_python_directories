#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#SBATCH --job-name=timeseries
#SBATCH --ntasks=1
#SBATCH --mem=8gb
#SBATCH --partition=interactive
#SBATCH --time=02:00:00
#SBATCH --output=Logs/bar_plots.log
import sys
import new_NEI as nei_species
import NEI_by_sector as nei
sys.path.append('/users/mjr583/EMEP_emissions')
import EMEP_by_sector as emep
import CEDS_by_sector as ceds
import read_tzompa_xiao as tzompa_xiao
import read_ceds
sys.path.append('/users/mjr583/python_lib')
import RowPy as rp
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib
import read as R
import re
from matplotlib.colors import LogNorm
import numpy as np
from netCDF4 import Dataset
matplotlib.use('agg')
plt.style.use('seaborn-darkgrid')

def EMEP_scale_factors( total, year ):
    sf = pd.read_csv("EMEP_NMVOC_ScaleFactors_allYears.csv", index_col=1)
    sf = sf[year]
    ef = pd.DataFrame(sf * total )
    ef.loc['Total',:]= ef.sum(axis=0)

    return ef, sf

def nei_data():
    return


def main():
    emep_sectors=False
    ceds_sectors=False
    ignore_sector=False
    start_year=2000
    end_year=2000
    only_land=False

    years=range(2017, 2018, 1)
    years=[ str(xx) for xx in years]
    
    for n, year in enumerate(years):
        print(year)
        start_year=year ; end_year=year
         
        vocs_ceds, voc_names = read_ceds.read_by_species(years=[year], lonlat=(-29.95, 89.95, 30.05, 89.95), is_land=only_land)

        vocs_ceds = np.array(vocs_ceds)
        vocs_ceds = np.nansum(np.nansum(np.nansum( vocs_ceds , 1), 1), 1)
        CEDS = pd.DataFrame( vocs_ceds, index = voc_names)
        CEDS_sf = CEDS / CEDS.sum(axis=0)
        
        EMEP_sf = pd.read_csv("EMEP_NMVOC_ScaleFactors_allYears.csv",index_col=1)
        EMEP_sf = EMEP_sf[year]
        print( EMEP_sf )
        print( type( EMEP_sf ))

        NEI_sf = pd.read_csv("NEI_NMVOC_ScaleFactors_2017.csv",index_col=0)['nei_spec']
        print( NEI_sf )
        print( type( NEI_sf) )
        print( NEI_sf.shape)

        ### Match up VOCs from seperate inventories ###
        EMEP_sf = EMEP_sf.transpose()
        NEI_sf  = NEI_sf.transpose()
        CEDS_sf = CEDS_sf.transpose()
        
        CEDS_sf['ALK4'] = CEDS_sf['butanes'] + CEDS_sf['hexanes_plus_higher_alkanes'] + CEDS_sf['pentanes'] 

        CEDS_sf = CEDS_sf[['methanal','alcohols','ethane','propane','propene','ketones','ALK4','benzene','toluene','xylene']]
        EMEP_sf = EMEP_sf[['CH2O','EOH','C2H6','C3H8','PRPE','MEK','ALK4','BENZ','TOLU','XYLE']]
        NEI_sf  = NEI_sf[['methanal','alcohols','ethane','propane','propene','ketones','ALK4','benzene','toluene','xylene']]
        
        CEDS_sf['alcohols'] = CEDS_sf['alcohols'] * 0.375

        CEDS_sf = CEDS_sf.transpose()
        EMEP_sf = EMEP_sf.transpose()
        NEI_sf  = NEI_sf.transpose()

        df = pd.DataFrame( {"CEDS" : np.array(np.squeeze(CEDS_sf.values))  ,
                            "EMEP" : np.array(np.squeeze(EMEP_sf.values))  
                                            }, index = EMEP_sf.index)

        df.index = ['CH2O','EOH','Ethane','Propane','Propene','Ketones','ALK4','Benzene','Toluene','Xylene']
        
        c = ['#66c2a5','#fc8d62','#8da0cb']
        #df=df.sort_index()
        print( df )
        f,ax=plt.subplots(figsize=(9,3.75))
        df.plot.bar(rot=0, ax=ax, color=c)
        plt.ylabel('NMVOC Speciation Factors')
        plt.title(year)
        plt.savefig(f'plots/EMEP_ScaleFactors_VOCs_{year}_%s.png' %str(only_land) )
        plt.close()

        print( df.EMEP / df.CEDS )


if __name__ == "__main__":
    main()
