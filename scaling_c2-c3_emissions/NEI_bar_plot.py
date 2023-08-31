#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#SBATCH --job-name=barplot
#SBATCH --ntasks=1
#SBATCH --mem=8gb
#SBATCH --partition=interactive
#SBATCH --time=02:00:00
#SBATCH --output=Logs/bar_plots_%A.log
import sys
sys.path.append('/users/mjr583/NEI')
import NEI_by_sector as nei
import new_NEI as nei_species
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

def bar_plot( ems, width=.4, labels='', suff=''):
    f,ax=plt.subplots(figsize=(14,7))
    
    for n, em in enumerate(ems):
        Y = np.array(np.squeeze(em.values))
        p1 =  ax.bar( range(len(ems[0].index)) + (n*width) , Y, width=.4, label=labels[n])
    ax.set_xticks(ems[0].index)
    ax.set_ylabel('VOC (Gg)')
    plt.legend() 
    plt.savefig( 'plots/bars_%s.png'  %suff)
    plt.close()
    return

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


if __name__ == "__main__":
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
        print("Get CEDS")
        vocs_ceds, voc_names, lon, lat = read_ceds.read_by_species(years=[year], lonlat=(-139.95, -50.05, 20.05, 59.95), is_land=only_land)
        # Gives CEDS VOC emissions with dimensions (VOC, Months, lat, lon)
        ceds = vocs_ceds
        vocs_ceds = np.array(vocs_ceds)
        vocs_ceds = np.nansum(np.nansum(np.nansum( vocs_ceds , 1), 1), 1)
        print(voc_names)

        cf = pd.DataFrame( vocs_ceds, index = voc_names)
        cf.loc['Total',:]= cf.sum(axis=0)

        print( "Get NEI")
        e1, T, Nlon, Nlat = nei.get(species='NMVOC', region='All_Eur', lonlat=(-139.95, -50.05, 20.05, 59.95), sectors=False, 
                                                    ignore_sector=False, only_land=only_land, make_map=True)
        
        print( "Make maps of CEDS and NEI NMVOCs")
        plot_map(np.nansum(np.nansum(ceds,0),0), lon, lat, lonlat=(-139.95, -50.05, 20.05, 59.95),
                                 out='CEDS_NMVOCs_%s.png' %str(only_land) )
        plot_map( e1, Nlon, Nlat, out='NEI_NMVOCs_%s.png' %str(only_land) )

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
        c = ['#66c2a5','#8da0cb']

        f,ax=plt.subplots(figsize=(12,5))
        df.plot.bar(rot=0, ax=ax, color=c)
        plt.ylabel('VOC (Gg $yr^{-1}$)')
        ax.set_yscale('log')
        plt.tight_layout()
        plt.savefig(f'bar_plots/LOGGY_NEI_VOCs_barplot_{year}_%s.png' %str(only_land))
        plt.close()
        print( "Ethane ratio: ", df["NEI"].values[0] / df["CEDS"].values[0] )
        print( "Propane ratio: ", df["NEI"].values[1] / df["CEDS"].values[1] )
