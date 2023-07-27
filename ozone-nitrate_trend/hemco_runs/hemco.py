
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('agg')
from datetime import datetime
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
        em = np.sum(em.reshape(-1, 12), axis=1)
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
    sectors = ['Total','Anthro']
    rundir='hemco_standalone'
    version='13.1.2'

    for sector in sectors:
        NO, times, lat, lon, lev, area, var_list, unit = GC.HEMCO_Diag_read(rundir, version, variable=f'EmisNO_{sector}', standalone=True)
        if len( NO.shape) == 4:      # If vertical levels then sum over the vertical
            NO = np.nansum( NO, 1 ) 
        NO = emission_to_tg( NO, area )
        NO, time = annual_means( NO, times=times )

        SO2, times, lat, lon, lev, area, var_list, unit = GC.HEMCO_Diag_read(rundir, version, variable=f'EmisSO2_{sector}', standalone=True)
        if len( SO2.shape) == 4:  
            SO2 = np.nansum( SO2, 1 )
        SO2 = emission_to_tg( SO2, area )

        df = pd.DataFrame( {'SO2' : SO2}, index=times )
        SO2 = annual_means( SO2 )
        
        f, ax = plt.subplots( figsize=(12,4) )
        ax.plot( time[:len(NO)], NO, label=f"{sector} NO" )
        ax.plot( time[:len(NO)], SO2, label=f"{sector} SO2" )
        
        a = list(SO2 - NO).index(max(list(SO2 - NO)))
        p = plt.axvline(time[a], color='#2ca02c', alpha=0.5)

        b = list(SO2 - NO).index(min(list(SO2 - NO)))

        print(len(time))
        print( a, time[a])
        print( b, time[b])
        
        c = list(SO2 - NO).index( min(list(SO2 - NO) , key=abs) )
        print(c, time[c])

        q = plt.axvline(time[b], color='red', alpha=0.5)

        r = plt.axvline(time[c], color='orange', alpha=0.5)
 

        plt.ylabel("Tg yr-1")
        plt.legend()
        plt.savefig(f'plots/{sector}_emissions_NO.png')
        plt.close()
