#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#SBATCH --job-name=timeseries
#SBATCH --ntasks=1
#SBATCH --mem=4gb
#SBATCH --partition=interactive
#SBATCH --time=00:20:00
#SBATCH --output=Logs/bar_plots.log
import sys
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
sys.path.append('/users/mjr583/NEI')
import NEI_by_sector as nei
import new_NEI as nei_species
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

    return ef

def EMEP_factors():
    plots_=[]
    emep_sectors=False
    ceds_sectors=False
    ignore_sector=False
    start_year=2017
    end_year=2017
    only_land=False

    years=range(2017, 2018, 1)
    years=[ str(xx) for xx in years]
    
    for n, year in enumerate(years):
        start_year=year ; end_year=year
        vocs_ceds, voc_names = read_ceds.read_by_species(years=[year], lonlat=(-29.95, 89.95, 30.05, 89.95), is_land=only_land)
        # Gives CEDS VOC emissions with dimensions (VOC, Months, lat, lon)
        vocs_ceds = np.array(vocs_ceds)
        vocs_ceds = np.nansum(np.nansum(np.nansum( vocs_ceds , 1), 1), 1)
        
        cf = pd.DataFrame( vocs_ceds, index = voc_names)
        cf.loc['Total',:]= cf.sum(axis=0)

        e1, time = emep.get_emissions( species='NMVOC', region='All_Eur', lonlat=(-29.95, 89.95, 30.05, 89.95), only_land=only_land, 
                                        start_year=start_year, end_year=end_year, sectors=emep_sectors, ignore_sector=ignore_sector,
                                        make_map=True) 
        NMVOC_EMEP= np.nansum(e1)
        ef = EMEP_scale_factors(NMVOC_EMEP, year)
        cf = cf.transpose()
        ef = ef.transpose()
        cf['ALK4'] = cf['butanes'] + cf['hexanes_plus_higher_alkanes'] + cf['pentanes'] 
        cf = cf[['ethane','propane','toluene','benzene','alcohols','xylene','ALK4','propene','trimethylbenzenes','ketones','methanal','Total']]
        cf['alcohols'] = cf['alcohols'] * 0.375
        ef['EOH'] = ef['EOH'] * 0.375
        ef = ef[['C2H6','C3H8','TOLU','BENZ','EOH','XYLE','ALK4', 'PRPE', 'TMB','MEK','CH2O','Total']]
        cf = cf.transpose()
        ef=ef.transpose()
        
        df = pd.DataFrame( {"CEDS VOCs" : np.array(np.squeeze(cf.values))  ,
                            "EMEP NMVOC (with NAEI speciation)" : np.array(np.squeeze(ef.values)) } ,
                            index = ef.index )
        df = df.transpose()
        df = df[['Total','CH2O','EOH','C2H6','C3H8','PRPE','MEK','ALK4','BENZ','TOLU','XYLE']]
        df=df.transpose()
        df.index = ['Total','CH2O','EOH','Ethane','Propane','Propene','Ketones','ALK4','Benzene','Toluene','Xylene']
        plots_.append( df )
    return df

def NEI_factors():
    emep_sectors=False
    ceds_sectors=False
    ignore_sector=False
    start_year=2000
    end_year=2000
    only_land=False

    years=range(2017, 2018, 1)
    years=[ str(xx) for xx in years]
    
    for n, year in enumerate(years):
        start_year=year ; end_year=year
        vocs_ceds, voc_names = read_ceds.read_by_species(years=[year], lonlat=(-139.95, -50.05, 20.05, 59.95), is_land=only_land)
        # Gives CEDS VOC emissions with dimensions (VOC, Months, lat, lon)
        ceds = vocs_ceds
        vocs_ceds = np.array(vocs_ceds)
        vocs_ceds = np.nansum(np.nansum(np.nansum( vocs_ceds , 1), 1), 1)

        cf = pd.DataFrame( vocs_ceds, index = voc_names)
        cf.loc['Total',:]= cf.sum(axis=0)

        e1, T, Nlon, Nlat = nei.get(species='NMVOC', region='All_Eur', lonlat=(-139.95, -50.05, 20.05, 59.95), sectors=False, 
                                                    ignore_sector=False, only_land=only_land, make_map=True)
        

        PRPA, t = nei_species.get_species(start_year=2017, end_year=2017, region='North_America', is_land=only_land,
                                          lonlat=(-139.95, -50.05, 20.05, 59.95), species='PRPA')

        ETHA, t = nei_species.get_species(start_year=2017, end_year=2017, region='North_America', is_land=only_land, 
                                          lonlat=(-139.95, -50.05, 20.05, 59.95), species='ETHA')
        RCHO, t = nei_species.get_species(start_year=2017, end_year=2017, region='North_America', is_land=only_land, 
                                          lonlat=(-139.95, -50.05, 20.05, 59.95), species='ALDX')
        BENZ, t = nei_species.get_species(start_year=2017, end_year=2017, region='North_America', is_land=only_land, 
                                          lonlat=(-139.95, -50.05, 20.05, 59.95), species='BENZ')
        EOH, t = nei_species.get_species(start_year=2017, end_year=2017, region='North_America', is_land=only_land, 
                                          lonlat=(-139.95, -50.05, 20.05, 59.95), species='ETOH')
        CH2O, t = nei_species.get_species(start_year=2017, end_year=2017, region='North_America', is_land=only_land, 
                                          lonlat=(-139.95, -50.05, 20.05, 59.95), species='FORM')
        MEK, t = nei_species.get_species(start_year=2017, end_year=2017, region='North_America', is_land=only_land, 
                                          lonlat=(-139.95, -50.05, 20.05, 59.95), species='KET')
        ALK4, t = nei_species.get_species(start_year=2017, end_year=2017, region='North_America', is_land=only_land, 
                                          lonlat=(-139.95, -50.05, 20.05, 59.95), species='PAR')

        TOLU, t = nei_species.get_species(start_year=2017, end_year=2017, region='North_America', is_land=only_land, 
                                          lonlat=(-139.95, -50.05, 20.05, 59.95), species='TOL')
        XYLE, t = nei_species.get_species(start_year=2017, end_year=2017, region='North_America', is_land=only_land, 
                                          lonlat=(-139.95, -50.05, 20.05, 59.95), species='XYLMN')

        OLE, t = nei_species.get_species(start_year=2017, end_year=2017, region='North_America', is_land=only_land, 
                                          lonlat=(-139.95, -50.05, 20.05, 59.95), species='OLE')
        IOLE, t = nei_species.get_species(start_year=2017, end_year=2017, region='North_America', is_land=only_land, 
                                          lonlat=(-139.95, -50.05, 20.05, 59.95), species='IOLE')
        PRPE = OLE + IOLE
        
        NMVOC_NEI= np.nansum(e1)
        nf = pd.DataFrame( {'Total' : NMVOC_NEI } , index=T )
        nf.loc['ethane',:]= np.nansum( ETHA )
        nf.loc['propane',:]= np.nansum( PRPA )
        nf.loc['benzene',:]= np.nansum( BENZ )
        nf.loc['EOH',:]= np.nansum( EOH )
        nf.loc['methanal',:]= np.nansum( CH2O )
        nf.loc['ketones',:]= np.nansum( MEK )
        nf.loc['ALK4',:]= np.nansum( ALK4 )
        nf.loc['toluene',:]= np.nansum( TOLU )
        nf.loc['xylene',:]= np.nansum( XYLE )
        nf.loc['propene',:]= np.nansum( PRPE )

        cf = cf.transpose()
        nf = nf.transpose()
        cf['ALK4'] = cf['butanes'] + cf['hexanes_plus_higher_alkanes'] + cf['pentanes']

        cf = cf[['methanal','alcohols','ethane','propane','propene','ketones','ALK4','benzene','toluene','xylene','Total']]
        nf  = nf[['methanal','EOH','ethane','propane','propene','ketones','ALK4','benzene','toluene','xylene']]
        cf['alcohols'] = cf['alcohols'] * 0.375
        nf['Total']= np.nansum(e1)
        cf = cf.transpose()
        nf = nf.transpose()
        
        df = pd.DataFrame( {"CEDS" : np.array(np.squeeze(cf.values))  ,
                            "NEI" : np.array(np.squeeze(nf.values)) } , 
                             index = cf.index)
        df=df.transpose()
        df = df[['Total','methanal','alcohols','ethane','propane','propene','ketones','ALK4','benzene','toluene','xylene']]
        df=df.transpose()
        df.index = ['Total','CH2O','EOH','Ethane','Propane','Propene','Ketones','ALK4','Benzene','Toluene','Xylene']
        return df

def main():
    df0 = EMEP_factors()
    print( df0 )

    sys.exit()

    df1 = NEI_factors()

    sf0=df0.transpose()
    sf0 = sf0.iloc[:,1:].div(sf0.Total, axis=0)
    sf0=sf0.transpose()

    sf1=df1.transpose()
    sf1 = sf1.iloc[:,1:].div(sf1.Total, axis=0)
    sf1=sf1.transpose()


    c = ['#66c2a5','#fc8d62']#,'#8da0cb']
    f,(ax1,ax2)=plt.subplots(2,1, figsize=(9,7.5))
    df0.plot.bar(rot=0, ax=ax1, width=.4, color=c)
    c = ['#66c2a5','#fc8d62']#,'#8da0cb']
    df1.plot.bar(rot=0, ax=ax2, width=.4, color=c)
    ax1.set_ylabel('VOC (Gg $yr^{-1}$)')
    ax2.set_ylabel('VOC (Gg $yr^{-1}$)')
    ax1.legend()
    ax2.legend()
    ax1.set_yscale('log')
    ax2.set_yscale('log')
    plt.tight_layout()
    plt.savefig(f'paper_figures/figure_2.png' )
    plt.close()

    c = ['#66c2a5','#fc8d62','#8da0cb']
    f,(ax1,ax2)=plt.subplots(2,1, figsize=(9,7.5))
    df0.plot.bar(rot=0, ax=ax1, width=.4, color=c)
    c = ['#66c2a5','#8da0cb']
    df1.plot.bar(rot=0, ax=ax2, width=.4, color=c)
    ax1.set_ylabel('VOC (Gg $yr^{-1}$)')
    ax1.set_ylabel('VOC (Gg $yr^{-1}$)')
    ax1.legend()
    ax2.legend()
    plt.tight_layout()
    plt.savefig(f'paper_figures/figure_2_nonlog.png' )
    plt.close()

    c = ['#66c2a5','#fc8d62','#8da0cb']
    f,(ax1,ax2)=plt.subplots(2,1, figsize=(9,7.5))
    sf0.plot.bar(rot=0, ax=ax1, width=.4, color=c)
    c = ['#66c2a5','#8da0cb']
    sf1.plot.bar(rot=0, ax=ax2, width=.4, color=c)
    ax1.set_ylabel('VOC scale factor')
    ax1.set_ylabel('VOC scale factor')
    ax1.legend()
    ax2.legend()
    plt.tight_layout()
    plt.savefig(f'paper_figures/figure_2b.png' )
    plt.close()



if __name__ == "__main__":
    main()

