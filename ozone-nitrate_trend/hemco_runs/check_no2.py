
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('agg')
from netCDF4 import Dataset
import sys
sys.path.append('/users/mjr583/python_lib')
import GC_tools as GC
plt.style.use("seaborn-darkgrid")

def emission_to_tg(ems, area, unit="kg/m2/s"):
    hold=[] 
    for em in ems:
        x = em * area           ## kg/m2/s to kg/s per gridbox
        x = x * 86400 * 30      ## kg/s/gridbox to kg/gridbox/month
        em =  x * 1e-9          ## kg to Tg
        hold.append( em )
    em = np.sum( np.sum( np.array(hold) , 1 ) , 1 )
    return em

def annual_means( em , times=False):
    try:
        em = np.mean(em.reshape(-1, 12), axis=1)
    except:
        for i in range(1,12):
            try:
                em = np.sum(em[:-i].reshape(-1, 12), axis=1)
                break
            except:
                continue
    if type(times) != bool:
        times = pd.to_datetime(times)
        times = times.year.drop_duplicates().tolist()
        return em, times
    return em


if __name__ == "__main__":
    sectors = ['Total']
    rundir='control'
    version='13.1.2'

    for sector in sectors:
        NO, times, lat, lon, lev, area, var_list, unit = GC.HEMCO_Diag_read(rundir, version, variable=f'EmisNO_{sector}', standalone=False)
        if len( NO.shape) == 4:      # If vertical levels then sum over the vertical
            NO = np.nansum( NO, 1 ) 
        NO = emission_to_tg( NO, area )
        NO, time = annual_means( NO, times=times )

        NO2_shp, times, lat, lon, lev, area, var_list, unit = GC.HEMCO_Diag_read(rundir, version, variable=f'EmisNO2_Ship', standalone=False)
        if len( NO2_shp.shape) == 4:  
            NO2_shp = np.nansum( NO2_shp, 1 ) 
        NO2_shp = emission_to_tg( NO2_shp, area )
        NO2_shp = annual_means( NO2_shp )

        NO2, times, lat, lon, lev, area, var_list, unit = GC.HEMCO_Diag_read(rundir, version, variable=f'EmisNO2_Anthro', standalone=False)
        if len( NO2.shape) == 4:  
            NO2 = np.nansum( NO2, 1 ) 
        NO2 = emission_to_tg( NO2, area )
        NO2 = annual_means( NO2 )

        NO2 = NO2 + NO2_shp
        
        print(NO.mean())
        print(NO2.mean())
        
        print( NO2 / NO * 100 )

        f, ax = plt.subplots( figsize=(12,4) )
        ax.plot( time[:len(NO)], NO, label=f"{sector} NO" )
        ax.plot( time[:len(NO)], NO2, label=f"{sector} NO2" )
        plt.ylabel("Tg yr-1")
        plt.legend()
        plt.savefig(f'plots/NO2_emissions.png')
        plt.close()
